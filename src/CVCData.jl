module CVCData

using Distributions, Random, SnpArrays, LinearAlgebra
using DataFrames, CSV, StatsBase, Roots
using DelimitedFiles: readdlm, writedlm
using JLD2
using TextWrap, Glob

export save_obj, CVCDataset, show

include("GenoRaw.jl")
include("Tools.jl")
include("Geno.jl")
include("GenoArch.jl")
include("GenoPart.jl")
include("Covariate.jl")

struct CVCDataset
    N::Int
    M::Int
    M_k::Vector{Int}
    C::Vector{Int}
    K::Int
    nldbins::Int
    nmafbins::Int
    h2::Real
    h2hat::Real
    cr::Vector{AbstractFloat}
    maflb::Real
    mafub::Real
    cvr::Real
    cvrhat::Real
    a::Real
    b::Real
    part_genodir::AbstractString
    covdir_parent::AbstractString
    covdir::Vector
    phenodir_parent::Vector{AbstractString}
    phenodir::Vector{AbstractString}
    ηw::Vector
    vc_array::AbstractArray
    rep::Vector{Int}
    datasetobjdir::AbstractString
    datasetobjpath::Vector
    normalize::Bool
    phig::Real
    phie::Real
    conditional_ldbin::Bool
    conditional_mafbin::Bool
    rng::Random.AbstractRNG
    ldtype::Symbol
    rho::AbstractFloat
end

function CVCDataset(
    ga::GenoArch,
    gp::GenoPart,
    covpath::AbstractString,
    covdir_parent::AbstractString,
    datasetobjdir::AbstractString
    
    )
    N = get_N(gp)
    M = get_M(gp)
    K = get_K(gp)
    covdir = Vector{AbstractString}(undef, 1)
    C = Vector{Int}(undef, 1)
    if isempty(covpath) || isempty(covdir_parent)
        covdir[] = ""
        T = Float64
        ηw = zeros(T, N)
        C[] = 0
    else
        mkpath(covdir_parent)
        @assert isfile(covpath) "Covariate file does not exist."
        covdir[], w = align_cov(covpath, gp.part_genodir[], covdir_parent, gp.ldtype, ga.rho)
        T = eltype(w)
        ηw = zeros(T, size(w, 1))
        update_fixed_effect_component!(ηw, w, ga.rng)
        C[] = size(w, 2)
    end
    nldbins = get_nldbins(gp)
    nmafbins = get_nmafbins(gp)
    h2 = get_h2(ga)
    h2hat = get_h2hat(ga)
    maflb = get_maflb(ga)
    mafub = get_mafub(ga)
    cvr = get_cvr(ga)
    cvrhat = round(get_cvrhat(ga), sigdigits=2)
    a = get_a(ga)
    b = get_b(ga)
    normalize = get_normalize(ga)
    phig = get_phig(ga)
    phie = get_phie(ga)
    part_genodir = get_part_genodir(gp)
    phenodir_parent = Vector{AbstractString}(undef, 1)
    phenodir = similar(phenodir_parent)
    phenodir_parent[] = ""
    phenodir[] = ""
    rep = Vector{Int}(undef, 1)
    M_k = get_m.(part_genodir, 1:K)
    vc_array = Array{Vector{T}, 1}(undef, K)
    vcmat = get_vcmat(ga)
    set_vcarray!(vc_array, vcmat, part_genodir)
    datasetobjpath = Vector{AbstractString}(undef, 1)
    conditional_ldbin = get_conditional_ldbin(gp)
    conditional_mafbin = get_conditional_mafbin(gp)
    cr = Vector{AbstractFloat}(undef, 1)
    cr[] = NaN
    CVCDataset(N, M, M_k, C, K, nldbins, nmafbins, h2, h2hat, cr, maflb, mafub, cvr, cvrhat, 
                a, b, part_genodir, covdir_parent, covdir, phenodir_parent, phenodir, ηw, 
                vc_array, rep, datasetobjdir, datasetobjpath, normalize, 
                phig, phie, conditional_ldbin, conditional_mafbin, ga.rng, gp.ldtype, ga.rho)
end

function set_vcarray!(
    vcarray::AbstractArray, 
    vcmat::DataFrame, 
    part_genodir::AbstractString
    )
    for k in eachindex(vcarray)
        bimk = DataFrame(readdlm("$(part_genodir)/G$k.bim", header = false), :auto)
        rename!(bimk, :x1 => :chr, :x2 => :snpid)
        vcarray[k] = leftjoin(bimk[!, [:chr, :snpid]], vcmat, on = [:chr, :snpid])[:,:vc]
    end
end

# TODO: create a function that subset and aligns covariate matrix and BED files by f.eid
function align_cov(
    covpath::AbstractString,
    part_genodir::AbstractString,
    covdir_parent::AbstractString,
    ldtype::Symbol,
    rho::Real=NaN
)
    w, header = readdlm("$(covpath)", header=true)
    fampath = get_fampath(part_genodir);
    k = get_k(part_genodir)
    fam = readdlm(fampath);
    feid = vec(string.(Int.(fam[:,1])));
    wdf = DataFrame(w, vec(header));
    wdf."f.eid" = string.(Int.(wdf."f.eid"))
    wdf = wdf[in.(wdf."f.eid", Ref(feid)), :]
    n = size(wdf, 1)
    c = size(wdf, 2) - 1
    covdirname = "sim_cov_N_$(n)_C_$(c)_K_$(k)_ldtype_$(ldtype)_rho_$(rho)"
    covdir = joinpath(covdir_parent, covdirname)
    mkpath(covdir)
    CSV.write("$(covdir)/w.txt", wdf[:,Not(:"f.eid")]; delim = "\t", header=false)
    w = Matrix(wdf[:,Not(:"f.eid")])
    insertcols!(wdf, 1, :FID => feid); rename!(wdf, :"f.eid" => :IID);
    CSV.write("$(covdir)/waug.txt", wdf; delim = "\t", header=false)
    covdir, w
end

function update_fixed_effect_component!(
    ηw::Vector{T},
    w::Matrix{T},
    rng::Random.AbstractRNG
    ) where T
    function simulate_alpha(
        rng::Random.AbstractRNG, 
        C::Int
        )
        Random.seed!(rng)
        [-C / 2 ÷ 2 + 0.5 * i for i in 1:C]
    end
    alpha = convert(Vector{T}, simulate_alpha(rng, size(w,2)))
    mul!(ηw, w, alpha)
end

function save_obj(s::CVCDataset)
    mkpath(s.datasetobjdir)
    set_datasetobjpath!(s)
    JLD2.save_object(s.datasetobjpath[], s)
end

function set_datasetobjpath!(s::CVCDataset)
    s.datasetobjpath[] = "$(s.datasetobjdir)/$(get_objname(s)).jld2"
end

function get_objname(s::CVCDataset)
    if isnan(s.cr[])
        "dataset_N_$(s.N)_M_$(s.M)_ldtype_$(s.ldtype)_rho_$(s.rho)_\
        C_$(s.C[1])_K_$(s.K)_nldbins_$(s.nldbins)_nmafbins_$(s.nmafbins)_\
        h2_$(s.h2)_cr_NA_maflb_$(s.maflb)_mafub_$(s.mafub)_cvr_$(s.cvr)_a_$(s.a)_b_$(s.b)"
    else
        "dataset_N_$(s.N)_M_$(s.M)_ldtype_$(s.ldtype)_rho_$(s.rho)_\
        C_$(s.C[1])_K_$(s.K)_nldbins_$(s.nldbins)_nmafbins_$(s.nmafbins)_\
        h2_$(s.h2)_cr_$(s.cr[])_maflb_$(s.maflb)_mafub_$(s.mafub)_cvr_$(s.cvr)_a_$(s.a)_b_$(s.b)"
    end
end

function Base.show(io::IO, s::CVCDataset)
    Mks = ["M$(k) = $(Mk)" for (k, Mk) in enumerate(s.M_k)]
    output = "N=$(s.N), M=$(s.M), $(join(Mks, ", ")), \
    C=$(s.C[1]), K=$(s.K), nldbins=$(s.nldbins), nmafbins=$(s.nmafbins), \
    ldtype=$(s.ldtype), rho = $(s.rho), h2=$(s.h2), h2hat=$(s.h2hat), cr=$(s.cr), \
    maflb=$(s.maflb), mafub=$(s.mafub), cvr=$(s.cvr), cvrhat=$(s.cvrhat), \
    a=$(s.a), b=$(s.b), normalize=$(s.normalize), \
    conditional_ldbin=$(s.conditional_ldbin), conditional_mafbin=$(s.conditional_mafbin) \
    phig=$(s.phig), phie=$(s.phie)"
    println("Showing simulation setting:")
    println_wrapped(output, width=50)
    println(io)
end

include("Simulate.jl")
include("Pheno.jl")


end
