module CVCData

using Distributions, Random, SnpArrays, LinearAlgebra
using DataFrames, CSV, StatsBase, Roots
using DelimitedFiles: readdlm, writedlm
using Serialization
using TextWrap, Glob, Printf

export save_obj, CVCDataset, show, align!

include("GenoRaw.jl")
include("Tools.jl")
include("Geno.jl")
include("GenoPart.jl")
include("GenoArch.jl")
include("Covariate.jl")

struct CVCDataset
    N::Int
    M::Int
    M_k::Vector{Int}
    C::Int
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
    covdir::AbstractString
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
    covdir::AbstractString,
    datasetobjdir::AbstractString
    )
    part_genodir = get_part_genodir(gp)
    @assert check_alignment(covdir, part_genodir) 
    "Subject IDs (rows) of covariate and genotype files are not aligned."
    N = get_N(gp)
    M = get_M(gp)
    K = get_K(gp)
    nldbins = get_nldbins(gp)
    nmafbins = get_nmafbins(gp)
    conditional_ldbin = get_conditional_ldbin(gp)
    conditional_mafbin = get_conditional_mafbin(gp)
    M_k = get_m.(part_genodir, 1:K)
    w = readdlm("$covdir/w.txt", '\t', header = false)
    ηw = zeros(eltype(w), size(w, 1))
    update_fixed_effect_component!(ηw, w, ga.rng)
    C = size(w, 2)
    vc_array = Array{Vector{eltype(w)}, 1}(undef, K)
    vcmat = get_vcmat(ga)
    set_vcarray!(vc_array, vcmat, part_genodir)
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
    phenodir_parent = Vector{AbstractString}(undef, 1)
    phenodir = similar(phenodir_parent)
    phenodir_parent[] = ""
    phenodir[] = ""
    rep = Vector{Int}(undef, 1); rep[] = zero(Int);
    datasetobjpath = Vector{AbstractString}(undef, 1)
    cr = Vector{AbstractFloat}(undef, 1)
    cr[] = NaN
    dataset = CVCDataset(N, M, M_k, C, K, nldbins, nmafbins, h2, h2hat, cr, maflb, mafub, cvr, cvrhat, 
                a, b, part_genodir, covdir, phenodir_parent, phenodir, ηw, 
                vc_array, rep, datasetobjdir, datasetobjpath, normalize, 
                phig, phie, conditional_ldbin, conditional_mafbin, ga.rng, gp.ldtype, ga.rho)
    save_obj(dataset)
    return dataset
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
    open("$(s.datasetobjpath[])", "w") do io
        Serialization.serialize(io, s)
    end
end

function set_datasetobjpath!(s::CVCDataset)
    s.datasetobjpath[] = "$(s.datasetobjdir)/$(get_objname(s)).jls"
end

function get_objname(s::CVCDataset)
    "dataset_N_$(s.N)_M_$(s.M)_ldtype_$(s.ldtype)_rho_$(s.rho)_\
    C_$(s.C[1])_K_$(s.K)_nldbins_$(s.nldbins)_nmafbins_$(s.nmafbins)_\
    h2_$(s.h2)_cr_$(s.cr[])_maflb_$(s.maflb)_mafub_$(s.mafub)_cvr_$(s.cvr)_a_$(s.a)_b_$(s.b)_rep_$(s.rep[])"
end

function Base.show(io::IO, s::CVCDataset)
    Mks = ["M$(k) = $(Mk)" for (k, Mk) in enumerate(s.M_k)]
    output = "N=$(s.N), M=$(s.M), $(join(Mks, ", ")), \
    C=$(s.C[1]), K=$(s.K), nldbins=$(s.nldbins), nmafbins=$(s.nmafbins), \
    ldtype=$(s.ldtype), rho = $(s.rho), h2=$(s.h2), h2hat=$(s.h2hat), cr=$(s.cr), \
    maflb=$(s.maflb), mafub=$(s.mafub), cvr=$(s.cvr), cvrhat=$(s.cvrhat), \
    a=$(s.a), b=$(s.b), normalize=$(s.normalize), \
    conditional_ldbin=$(s.conditional_ldbin), conditional_mafbin=$(s.conditional_mafbin) \
    phig=$(s.phig), phie=$(s.phie), rep=$(s.rep)"
    println("Showing simulation setting:")
    println_wrapped(output, width=50)
    println(io)
end

function align!(
    gp::GenoPart,
    covdir::AbstractString;
    subjectidspath::Union{AbstractString, Nothing} = nothing,
    phenopath::Union{AbstractString, Nothing} = nothing
    )
    covpath = "$covdir/wdf.txt"
    wdf = CSV.read(covpath, DataFrame)
    id_cov = string.(Int.(wdf."f.eid"))
    part_genodir = gp.part_genodir[]
    ldtype = gp.ldtype
    rho = gp.rho
    plink_binpath = gp.plink_binpath
    fampath = get_fampath(part_genodir);
    fam = readdlm(fampath);
    id_geno = vec(string.(Int.(fam[:,1])));
    @assert id_geno != 0 "No individuals in the genotype file."    
    
    # Initialize with intersection of covariate and genotype IDs
    id = intersect(id_cov, id_geno)
    
    # Add subject ID intersection if provided
    id_subject = nothing
    if !isnothing(subjectidspath) && isfile(subjectidspath)
        subject_data = readdlm(subjectidspath)
        id_subject = vec(string.(Int.(subject_data[:,1])))
        id = intersect(id, id_subject)
    end
    
    # Add phenotype ID intersection if provided
    id_pheno = nothing
    pheno_df = nothing
    pheno_basename = nothing
    pheno_dirname = nothing
    if !isnothing(phenopath) && isfile(phenopath)
        pheno_df = CSV.read(phenopath, DataFrame)
        pheno_basename = basename(phenopath)
        pheno_dirname = basename(dirname(phenopath))
        # Assume first column contains IDs - adjust if needed
        id_pheno = string.(Int.(pheno_df[:,1]))
        id = intersect(id, id_pheno)
    end
    
    n = length(id)
    @assert n > 0 "No individuals remain after intersection of all provided ID sets."
    
    # Check if ACTUAL alignment is needed - only true if any id set differs from intersection
    need_alignment = (id != id_cov) || (id != id_geno) || 
                     (!isnothing(id_subject) && id != id_subject) || 
                     (!isnothing(id_pheno) && id != id_pheno)
    
    # If no alignment is needed and no additional files, just return original directories
    if !need_alignment && isnothing(phenopath) && isnothing(subjectidspath)
        return covdir, part_genodir
    end
    
    K = get_k(part_genodir)
    covdirparent = splitdir(covdir)[1]
    partgenodirparent = splitdir(part_genodir)[1]
    
    # Create directory name components based on what's being aligned
    c = size(wdf, 2) - 1
    covdirname = "aligned_cov_N_$(n)_C_$(c)_K_$(K)_rho_$(rho)_ldtype_$(ldtype)"
    out_covdir = joinpath(covdirparent, covdirname)
    
    genodirname = "aligned_partgeno_N_$(n)_C_$(c)_K_$(K)_rho_$(rho)_ldtype_$(ldtype)"
    out_partgenodir = joinpath(partgenodirparent, genodirname)
    
    # Only create output directories if alignment is needed
    if need_alignment
        mkpath(out_covdir)
        mkpath(out_partgenodir)
        
        # Create the aligned ID file
        CSV.write("$out_partgenodir/id.txt", DataFrame(FID=id, IID=id); delim = "\t", header=false)
        
        # Filter covariate data
        mask = in.(string.(wdf."f.eid"), Ref(id))
        w = wdf[mask, Not(:"f.eid")]
        wdf = wdf[mask, :]
        
        # Write covariate data
        CSV.write("$out_covdir/wdf.txt", wdf; delim = "\t", header=true)
        CSV.write("$out_covdir/w.txt", w; delim = "\t", header=false)
        
        # Create filtered genotype files
        for k in 1:K
            plink_cmd = `$plink_binpath/plink2 \
            --bfile $part_genodir/G$(k) \
            --keep-fam $out_partgenodir/id.txt \
            --make-bed \
            --out $out_partgenodir/G$(k)`
            run(plink_cmd)
        end
    else
        # No alignment needed, but we still need directories for additional files
        mkpath(out_covdir)
        mkpath(out_partgenodir)
        CSV.write("$out_partgenodir/id.txt", DataFrame(FID=id, IID=id); delim = "\t", header=false)
        
        # Create symbolic links for covariate files
        for file in ["wdf.txt", "w.txt"]
            src = joinpath(covdir, file)
            dst = joinpath(out_covdir, file)
            isfile(dst) && rm(dst)
            isfile(src) && symlink(src, dst)
        end

        # Create symbolic links for genotype files
        for k in 1:K
            for ext in [".bed", ".bim", ".fam"]
            src = joinpath(part_genodir, "G$k$ext")
            dst = joinpath(out_partgenodir, "G$k$ext")
            isfile(dst) && rm(dst)
            isfile(src) && symlink(src, dst)
            end
        end
    end
    
    # Handle phenotype data if provided
    out_phenodir = nothing
    if !isnothing(phenopath)
        # Get phenotype parent directory automatically
        phenodir_parent = dirname(phenopath)
        
        # Create a simple aligned phenotype directory name
        out_phenodir = joinpath(phenodir_parent, "aligned_$(pheno_dirname)")
        mkpath(out_phenodir)
        
        if need_alignment && !isnothing(pheno_df)
            # Filter phenotype data
            pheno_mask = in.(string.(Int.(pheno_df[:,1])), Ref(id))
            pheno_filtered = pheno_df[pheno_mask, :]
            CSV.write(joinpath(out_phenodir, pheno_basename), pheno_filtered; delim = "\t")
        else
            # Just create a symbolic link
            dst = joinpath(out_phenodir, pheno_basename)
            isfile(dst) && rm(dst)
            symlink(phenopath, dst)
        end
    end
    
    gp.part_genodir[] = out_partgenodir
    
    # Return appropriate directories
    if !isnothing(out_phenodir)
        return out_covdir, out_partgenodir, out_phenodir
    else
        return out_covdir, out_partgenodir
    end
end

function check_alignment(
    covdir::AbstractString,
    partgenodir::AbstractString
    )
    wdf = CSV.read("$covdir/wdf.txt", DataFrame);
    id_cov = string.(Int.(wdf."f.eid"));
    fampath = get_fampath(partgenodir);
    id_geno = vec(string.(Int.(readdlm(fampath)[:,1])));
    return id_cov == id_geno
end
include("Simulate.jl")
include("Pheno.jl")


end
