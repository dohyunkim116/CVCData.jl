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
    subjectidspath::Union{AbstractString, Nothing}=nothing
    )
    # Progress tracking
    function status(msg)
        println(">>> ALIGN: $msg")
        flush(stdout)
    end
    
    status("Starting alignment process")
    
    # Add thread detection at the beginning
    nthreads = haskey(ENV, "NSLOTS") ? parse(Int, ENV["NSLOTS"]) : Threads.nthreads()
    status("Detected $nthreads threads for PLINK operations")
    
    # Read covariate IDs
    status("Reading covariate data from $covdir/wdf.txt")
    covpath = "$covdir/wdf.txt"
    wdf = CSV.read(covpath, DataFrame)
    id_cov = string.(Int.(wdf."f.eid"))
    status("Found $(length(id_cov)) individuals in covariate file")

    # Read genotype IDs
    part_genodir = gp.part_genodir[]
    fampath = get_fampath(part_genodir)
    status("Reading genotype data from $fampath")
    fam = readdlm(fampath)
    id_geno = vec(string.(Int.(fam[:,1])))
    status("Found $(length(id_geno)) individuals in genotype file")
    @assert !isempty(id_geno) "No individuals in the genotype file."

    # Intersection of cov & geno IDs
    status("Computing intersection of IDs")
    id = intersect(id_cov, id_geno)
    status("$(length(id)) individuals in common between covariates and genotypes")

    # Optional subject ID filtering
    if !isnothing(subjectidspath) && isfile(subjectidspath)
        status("Applying additional subject ID filtering from $subjectidspath")
        subj = readdlm(subjectidspath)
        id_subj = vec(string.(Int.(subj[:,1])))
        status("Found $(length(id_subj)) IDs in subject filter file")
        id = intersect(id, id_subj)
        status("$(length(id)) individuals remain after all filtering")
    end

    @assert !isempty(id) "No individuals remain after intersection."

    # Determine if alignment is needed via element‐wise checks
    status("Checking if alignment is needed")
    cov_ok  = length(id) == length(id_cov)  && all(id .== id_cov)
    geno_ok = length(id) == length(id_geno) && all(id .== id_geno)
    # if subjectidspath was used, check that too
    subj_ok = true
    if !isnothing(subjectidspath) && isfile(subjectidspath)
        id_subj = vec(string.(Int.(readdlm(subjectidspath)[:,1])))
        subj_ok = length(id) == length(id_subj) && all(id .== id_subj)
    end
    need_alignment = !(cov_ok && geno_ok && subj_ok)
    
    if need_alignment
        status("Alignment needed: covariate_aligned=$(cov_ok), genotype_aligned=$(geno_ok), subject_filter_aligned=$(subj_ok)")
    else
        status("No alignment needed - all files already have matching IDs")
    end

    # Prepare output paths
    n = length(id)
    c = size(wdf,2)-1
    K = get_k(part_genodir)
    covdir_parent = dirname(covdir)
    geno_parent   = dirname(part_genodir)

    status("Preparing output paths for $n individuals and $c covariates")
    aligned_cov_out  = joinpath(covdir_parent,
        "aligned_cov_N_$(n)_C_$(c)_K_$(K)_rho_$(gp.rho)_ldtype_$(gp.ldtype)")
    aligned_partgeno_out = joinpath(geno_parent,
        "aligned_partgeno_N_$(n)_C_$(c)_K_$(K)_rho_$(gp.rho)_ldtype_$(gp.ldtype)")
    
    status("Output paths:\n  - Covariates: $aligned_cov_out\n  - Genotypes: $aligned_partgeno_out")

    if need_alignment
        status("Starting alignment process")
        mkpath(aligned_cov_out);  mkpath(aligned_partgeno_out)
        
        # write aligned IDs and filter
        status("Filtering and writing covariate files")
        mask = in.(string.(wdf."f.eid"), Ref(id))
        wdf = wdf[mask, :]
        id_to_row = Dict(string(eid) => i for (i, eid) in enumerate(wdf."f.eid"))
        ordered_indices = [id_to_row[i] for i in id]
        wdf = wdf[ordered_indices, :]
        w = wdf[:, Not(:"f.eid")]        
        CSV.write("$aligned_cov_out/wdf.txt", wdf; delim="\t", header=true)
        CSV.write("$aligned_cov_out/w.txt", w; delim="\t", header=false)

        status("Writing aligned ID file for PLINK filtering")
        CSV.write("$aligned_partgeno_out/id.txt", DataFrame(FID=id, IID=id); delim="\t", header=false)

        # PLINK‐filter genotype with threading
        status("Starting PLINK filtering for $K partitions using $nthreads threads")
        for k in 1:K
            status("  - Filtering partition $k/$K with PLINK")
            cmd = `$(gp.plink_binpath)/plink2 \
                   --threads $nthreads \
                   --bfile $part_genodir/G$(k) \
                   --keep-fam $aligned_partgeno_out/id.txt \
                   --indiv-sort f \
                   --make-bed \
                   --out $aligned_partgeno_out/G$(k)`
            run(cmd)
            status("  - Completed partition $k/$K")
        end
    else
        # no alignment → symlink originals
        status("Creating symbolic links (no filtering required)")
        mkpath(aligned_cov_out);  mkpath(aligned_partgeno_out)
        
        status("Writing ID file and creating covariate symlinks")
        CSV.write("$aligned_partgeno_out/id.txt", DataFrame(FID=id, IID=id);
                  delim="\t", header=false)
        for f in ["wdf.txt","w.txt"]
            src, dst = joinpath(covdir,f), joinpath(aligned_cov_out,f)
            isfile(dst) && rm(dst);  isfile(src) && symlink(src,dst)
        end
        
        status("Creating genotype file symlinks")
        for k in 1:K, ext in [".bed",".bim",".fam"]
            src, dst = joinpath(part_genodir,"G$k$ext"), joinpath(aligned_partgeno_out,"G$k$ext")
            isfile(dst) && rm(dst);  isfile(src) && symlink(src,dst)
        end
    end
    status("Verifying alignment between covariate and genotype files")
    if !check_alignment(aligned_cov_out, aligned_partgeno_out)
        status("Alignment verification failed, cleaning up")
        rm(aligned_cov_out, recursive=true)
        rm(aligned_partgeno_out, recursive=true)
        error("Alignment check failed: covariate and genotype IDs do not match after alignment.")
    end
    status("Alignment verification successful")

    # Create a deep copy of the original GenoPart object
    status("Creating and updating aligned GenoPart object")
    aligned_gp = deepcopy(gp)
    
    # Update the aligned GenoPart with the new sample size and directory
    status("Updating sample count to $n and directory path")
    update_N!(aligned_gp, n)  # Update sample size to the aligned count
    aligned_gp.part_genodir[] = aligned_partgeno_out
    
    # Save the aligned GenoPart object
    aligned_objpath = joinpath(aligned_gp.part_genoobj_dir, 
        "aligned_$(basename(get_objname(aligned_gp))).jls")
    status("Saving aligned GenoPart to $aligned_objpath")
    open(aligned_objpath, "w") do io
        serialize(io, aligned_gp)
    end
    
    status("Alignment complete")
    
    # Return both directories and the aligned GenoPart object
    return aligned_cov_out, aligned_partgeno_out, aligned_gp
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
