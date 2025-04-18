export Geno, update_ldmafmat!, get_objpath

struct Geno
    plink_binpath::AbstractString
    gcta_binpath::AbstractString
    ldak_binpath::AbstractString
    raw_genodir::AbstractString
    raw_genobasename::AbstractString
    qced_genodir_parent::AbstractString
    qced_genodir::Vector{AbstractString}
    qced_genoobj_dir::AbstractString
    chrs::Union{Vector{Int}, UnitRange{Int}, Int}
    geno::AbstractFloat  # maximum missing genotype rate allowed
    maf::AbstractFloat # minimum MAF allowed
    ldmafmat::Vector{DataFrame}
    N::Vector{Int}
    M::Vector{Int}
    ldtype::Vector{Symbol}
    rho::AbstractFloat
    imputed::Bool
    ldcomputed::Vector{Bool}
    objpath::Vector{AbstractString}
end

function Geno(
    gs::Vector{Geno},
    qced_genodir_parent::AbstractString,
    qced_genoobj_dir::AbstractString
)
    @assert check_merge_condition(gs)
    ldmafmat = sizehint!(DataFrame[], 1)
    ldmafmat_array = [g.ldmafmat[1] for g in gs]
    push!(ldmafmat, reduce(vcat, ldmafmat_array))
    plink_binpath = gs[1].plink_binpath
    gcta_binpath = gs[1].gcta_binpath
    ldak_binpath = gs[1].ldak_binpath
    raw_genodir = gs[1].raw_genodir
    raw_genobasename = gs[1].raw_genobasename
    qced_genodir = Vector{AbstractString}(undef, 1)
    chrs = reduce(vcat, [g.chrs for g in gs]) |> sort |> unique
    geno = gs[1].geno
    maf = gs[1].maf
    N = gs[1].N
    M = sum([g.M for g in gs])
    ldtype = gs[1].ldtype
    rho = gs[1].rho
    imputed = gs[1].imputed
    ldcomputed = gs[1].ldcomputed
    objpath = Vector{AbstractString}(undef, 1)
    g = Geno(plink_binpath, gcta_binpath, ldak_binpath, raw_genodir, raw_genobasename,
        qced_genodir_parent, qced_genodir, qced_genoobj_dir, chrs, geno, maf, ldmafmat, 
        N, M, ldtype, rho, imputed, ldcomputed, objpath)
    update_qced_genodir_merged!(g, gs)
    g.ldcomputed[1] ? _update_ldmafmat_rowid!(g) : nothing
    update_objpath!(g)
    save_obj(g)
    for g in gs
        rm(g.objpath[1])
    end
    g
end

function Geno(
    plink_binpath::AbstractString, 
    gcta_binpath::AbstractString, 
    ldak_binpath::AbstractString, 
    raw_genodir::AbstractString, 
    raw_genobasename::AbstractString,
    qced_genodir_parent::AbstractString, 
    qced_genoobj_dir::AbstractString,
    chrs::Union{Vector{Int}, UnitRange{Int}, Int},
    geno::AbstractFloat, 
    maf::AbstractFloat,
    rho::AbstractFloat=NaN,
    imputed::Bool=false;
    subjectids_keep_path::AbstractString="",
    sub_sample_size::Int=get_n(raw_genodir;imputed=imputed),
    snpidspath_template::AbstractString="" 
)
    @assert geno >= 0.0 && geno <= 1.0 && maf >= 0.0 && maf <= 0.5 
        "geno should be in [0,1], and maf should be in [0,0.5]"
    chrs = sort(unique(chrs))
    qced_genodir = Vector{AbstractString}(undef, 1)
    ldmafmat = Vector{DataFrame}(undef, 1)
    N = Vector{Int}(undef, 1)
    M = Vector{Int}(undef, 1)
    ldcomputed = Vector{Bool}(undef, 1)
    ldcomputed[1] = false
    ldtype = Vector{Symbol}(undef, 1)
    ldtype[1] = :none
    objpath = Vector{AbstractString}(undef, 1)
    g = Geno(plink_binpath, gcta_binpath, ldak_binpath, raw_genodir, raw_genobasename,
        qced_genodir_parent, qced_genodir, qced_genoobj_dir, chrs, geno, maf, ldmafmat, 
        N, M, ldtype, rho, imputed, ldcomputed, objpath)
    qc_snps!(
        g,
        subjectids_keep_path,
        sub_sample_size,
        snpidspath_template
    );
    update_qced_genodir!(g);
    update_objpath!(g)
    save_obj(g)
    g
end

function qc_snps!(
    g::Geno,
    subjectids_keep_path::AbstractString,
    sub_sample_size::Int,
    snpidspath_template::AbstractString
)
    mkpath(g.qced_genodir_parent)
    g.qced_genodir[1] = mktempdir(g.qced_genodir_parent; cleanup=false)
    subjectids_path = get_subjectids_path(
        g.raw_genodir, 
        g.qced_genodir[1], 
        subjectids_keep_path, 
        sub_sample_size;
        imputed=g.imputed
    )
    qc_snps!(g, subjectids_path, snpidspath_template)
    g.N[1] = get_n(g.qced_genodir[1])
    g.M[1] = sum(get_m.(g.qced_genodir[1], g.chrs))
end

function qc_snps!(
    g::Geno,
    subjectids_path::AbstractString,
    snpidspath_template::AbstractString
)
    for chr in g.chrs
        snpids_path = get_snpids_path(snpidspath_template, chr, g.imputed)
        _qc_snps(
            g.plink_binpath, 
            g.raw_genodir, 
            g.raw_genobasename, 
            g.qced_genodir[1],
            subjectids_path,
            snpids_path,
            chr;
            geno=g.geno, 
            maf=g.maf,
            imputed=g.imputed
        )
    end
end

function _qc_snps(
    plink_binpath::AbstractString,
    genodir::AbstractString,
    genobasename::AbstractString,
    qced_genodir::AbstractString,
    subjectids_path::AbstractString,
    snpids_path::AbstractString,
    chr::Int;
    geno::AbstractFloat=0.01,
    maf::AbstractFloat=0.01,
    imputed::Bool=false
)
    bedfile_basename = "$(qced_genodir)/G$chr"
    
    # Common filters applied to both imputed and non-imputed data
    common_filters = [
        "--geno $geno",
        "--maf $maf minor",
        "--hwe midp 1e-7",
        "--autosome",
        "--snps-only just-acgt"
    ]
    
    if imputed
        # Step 1: Import data from BGEN format
        import_cmd_args = [
            "--bgen $(genodir)/$(genobasename)$chr.bgen 'ref-first'",
            "--sample $(genodir)/$(genobasename)$chr.sample",
            "--hard-call-threshold 0.1",
            "--import-dosage-certainty 0.9",
            "--mach-r2-filter 0.3",
            "--rm-dup 'exclude-mismatch' 'list'",
            "--keep-fam $subjectids_path"
        ]
        
        # Add SNP extraction if provided
        !isempty(snpids_path) && push!(import_cmd_args, "--extract $snpids_path")
        
        # Create BED files
        push!(import_cmd_args, "--make-bed", "--out $bedfile_basename")
        run(`$plink_binpath/plink2 $(split(join(import_cmd_args, " ")))`)
        
        # Step 2: Apply common filters
        filter_cmd_args = ["--bfile $bedfile_basename"]
        append!(filter_cmd_args, common_filters)
        push!(filter_cmd_args, "--make-bed", "--out $bedfile_basename")
        run(`$plink_binpath/plink2 $(split(join(filter_cmd_args, " ")))`)
        
        # Clean up temporary files
        rm.(filter(f -> occursin(r".bim~|.fam~|.bed~", f), 
            readdir(qced_genodir, join=true)), force=true)
    else
        # For non-imputed data: Single command with all operations
        cmd_args = [
            "--bfile $(genodir)/$(genobasename)$chr",
            "--keep-fam $subjectids_path"
        ]
        
        # Add SNP extraction if provided
        !isempty(snpids_path) && push!(cmd_args, "--extract $snpids_path")
        
        # Add common filters and output options
        append!(cmd_args, common_filters)
        push!(cmd_args, "--make-bed", "--out $bedfile_basename")
        
        run(`$plink_binpath/plink2 $(split(join(cmd_args, " ")))`)
    end
end

function get_subjectids_path(
    raw_genodir::AbstractString,
    qced_genodir::AbstractString,
    subjectids_keep_path::AbstractString,
    sub_sample_size::Int;
    imputed::Bool=false
    )
    if isempty(subjectids_keep_path)
        subjectids_keep = []
    else 
        subjectids_keep = vec(string.(Int.(readdlm(subjectids_keep_path))))
    end
    fampath = get_fampath(raw_genodir; imputed=imputed)
    if !imputed
        subjectids = string.(Int.(vec(readdlm(fampath)[:, 1])))
    else
        subjectids = string.(Int.(vec(readdlm(fampath)[: ,1])[2:end]))
    end
    subjectids = isempty(subjectids_keep) ? subjectids : intersect(subjectids, subjectids_keep)
    n = length(subjectids)
    if sub_sample_size > n
        sub_sample_size = n
    end
    @assert sub_sample_size <= n && sub_sample_size > 0 
        "Subsample size must be greater than zero and less than the size of subsetted samples."
    rowindex = sample(1:n, sub_sample_size, replace=false)
    subjectids_path = "$(qced_genodir)/subjectids.txt"
    open(subjectids_path, "w+") do io
        writedlm(io, subjectids[rowindex])
    end
    subjectids_path
end

function get_snpids_path(snpidspath_template::AbstractString, chr::Int=0, imputed::Bool=false)
    # Return empty string if snpidspath_template is empty
    isempty(snpidspath_template) && return ""
    
    # For chromosome-specific requests
    if chr > 0
        return replace(snpidspath_template, "{chr}" => chr)
    else
        # If requesting a global file but we're using per-chromosome files,
        # just return empty string for both imputed and non-imputed data
        return occursin("{chr}", snpidspath_template) ? "" : snpidspath_template
    end
end

function update_qced_genodir!(g::Geno)
    if g.imputed
        newdir = "$(g.qced_genodir_parent)/qced_imputed_geno_chrs_$(g.chrs[1])_to_$(g.chrs[end])\
        _N_$(g.N[1])_M_$(g.M[1])_ldtype_$(g.ldtype[1])"
    else
        newdir = "$(g.qced_genodir_parent)/qced_geno_chrs_$(g.chrs[1])_to_$(g.chrs[end])\
        _N_$(g.N[1])_M_$(g.M[1])_ldtype_$(g.ldtype[1])"
    end
    mv("$(g.qced_genodir[1])", newdir, force=true)
    g.qced_genodir[1] = newdir
end

function update_ldmafmat!(
    g::Geno;
    ldtype::Symbol=:ldak,
    ldwind::Int=2000,
    isannot::Bool=false
)
    if !isannot
        @assert ldtype == :ldsc || ldtype == :ldak
        lddata_array = Vector{Matrix{Any}}(undef, length(g.chrs));
        header = [];
        for (i, chr) in enumerate(g.chrs)
            lddata_array[i], header = 
                _compute_ld(
                    g.qced_genodir[1], 
                    "G", 
                    chr,
                    g.qced_genodir[1],
                    g.gcta_binpath,
                    g.ldak_binpath; 
                    ldtype=ldtype,
                    ldwind=ldwind
                )
        end
        chr = reduce(vcat, [repeat([chr], size(lddata_array[i], 1)) for (i, chr) in enumerate(g.chrs)])
        data = reduce(vcat, lddata_array)
        g.ldmafmat[1] = DataFrames.DataFrame(data, vec(header))
        if ldtype == :ldsc
            rename!(g.ldmafmat[1], :ldscore => :ldsc, :MAF => :maf, :SNP => :snpid)
            g.ldmafmat[1][!, :ldsc] = convert(Vector{Float64}, g.ldmafmat[1][!, :ldsc])
            g.ldmafmat[1][!, :maf] = convert(Vector{Float64}, g.ldmafmat[1][!, :maf])
            g.ldmafmat[1][!, :snpid] = convert(Vector{String}, g.ldmafmat[1][!, :snpid])
            g.ldmafmat[1][!, :chr] = chr
            select!(g.ldmafmat[1], :chr, :snpid, :maf, :ldsc)
        else
            rename!(g.ldmafmat[1], :Weight => :ldak, :Predictor => :snpid)
            g.ldmafmat[1][!, :ldak] = convert(Vector{Float64}, g.ldmafmat[1][!, :ldak])
            g.ldmafmat[1][!, :snpid] = convert(Vector{String}, g.ldmafmat[1][!, :snpid])
            g.ldmafmat[1][!, :chr] = chr
            select!(g.ldmafmat[1], :chr, :snpid, :ldak)
        end
        if ldtype == :ldak
            _update_ldmafmat_maf!(g)
        end
        _update_ldmafmat_rowid!(g)
        g.ldcomputed[1] = true
        g.ldtype[1] = ldtype
        update_qced_genodir!(g)
        rm(g.objpath[1])
        update_objpath!(g)
        save_obj(g)
    else
        lddata_array = Vector{Matrix{Any}}(undef, length(g.chrs));
        header = [];
        for (i, chr) in enumerate(g.chrs)
            lddata_array[i], header = 
                load_bim(g.qced_genodir[1], "G", chr)
        end
        chr = reduce(vcat, [repeat([chr], size(lddata_array[i], 1)) for (i, chr) in enumerate(g.chrs)])
        data = reduce(vcat, lddata_array)
        g.ldmafmat[1] = DataFrames.DataFrame(data, vec(header))
        g.ldmafmat[1][!, :chr] = chr
        g.ldmafmat[1][!, :snpid] = convert(Vector{String}, g.ldmafmat[1][!, :snpid])
        g.ldmafmat[1][!, :cm] = convert(Vector{Float64}, g.ldmafmat[1][!, :cm])
        g.ldmafmat[1][!, :bp] = convert(Vector{Int}, g.ldmafmat[1][!, :bp])
        g.ldmafmat[1][!, :a1] = convert(Vector{String}, g.ldmafmat[1][!, :a1])
        g.ldmafmat[1][!, :a2] = convert(Vector{String}, g.ldmafmat[1][!, :a2])
        save_obj(g)
    end
    
end

function load_bim(
    genodir::AbstractString,
    genobasename::AbstractString,
    chr::Int
)
    data = readdlm("$(genodir)/$(genobasename)$chr.bim", header=false)
    header = [:chr, :snpid, :cm, :bp, :a1, :a2]
    data, header
end

function _update_ldmafmat_maf!(g::Geno)
    function construct_sa_array!(
        sa_array::Array{SnpArray,1}, 
        genodir::AbstractString,
        genobasename::AbstractString,
        parts::Union{Vector{Int},Vector{String},UnitRange{Int},Int}
        )
        for (i, k) in enumerate(parts)
            sa_array[i] = SnpArray("$(genodir)/$(genobasename)$(k).bed")
        end
    end
    sa_array = Vector{SnpArray}(undef, length(g.chrs));
    construct_sa_array!(sa_array, g.qced_genodir[1], "G", g.chrs)
    g.ldmafmat[1][!, :maf] = reduce(vcat, maf.(sa_array))
end

function _compute_ld(
    genodir::AbstractString,
    genobasename::AbstractString,
    chr::Int,
    tempdir_path::AbstractString,
    gcta_binpath::AbstractString,
    ldak_binpath::AbstractString;
    ldtype::Symbol=:ldsc,
    ldwind::Int=2000
    )
    @assert ldtype == :ldsc || ldtype == :ldak
    if ldtype == :ldsc
        gctacmd = `$gcta_binpath/gcta-1.94.1 \
                    --bfile $genodir/$genobasename$chr \
                    --ld-score \
                    --ld-wind $ldwind \
                    --out $tempdir_path/$genobasename$chr`
        run(gctacmd)
        data, header = readdlm("$tempdir_path/$(genobasename)$(chr).score.ld", header=true);
    else
        ldakcmd1 = `$ldak_binpath/ldak5.2  \
                    --cut-weights $tempdir_path/sections$chr \
                    --bfile $genodir/$genobasename$chr`
        ldakcmd2 = `$ldak_binpath/ldak5.2 \
                    --calc-weights-all $tempdir_path/sections$chr \
                    --bfile $genodir/$genobasename$chr`
        run(ldakcmd1); run(ldakcmd2);
        data, header = readdlm("$tempdir_path/sections$(chr)/weights.all", header=true);
    end
    data, header
end

_update_ldmafmat_rowid!(g::Geno) = g.ldmafmat[1][!, :rowid] = 1:size(g.ldmafmat[1], 1);

function save_obj(g::Geno)
    mkpath(g.qced_genoobj_dir)
    #JLD2.save_object(g.objpath[1], g)
    open("$(g.objpath[])", "w") do io
        Serialization.serialize(io, g)
    end
end

function update_objpath!(g::Geno)
    g.objpath[] = "$(g.qced_genoobj_dir)/$(get_objname(g)).jls"
    g.objpath[]
end

function get_objname(g::Geno)
    if g.imputed
        return "qced_imputed_geno_chrs_$(g.chrs[1])_to_$(g.chrs[end])_N_$(g.N[1])_M_$(g.M[1])_ldtype_$(g.ldtype[1])"
    else
        return "qced_geno_chrs_$(g.chrs[1])_to_$(g.chrs[end])_N_$(g.N[1])_M_$(g.M[1])_ldtype_$(g.ldtype[1])"
    end
end

get_N(g::Geno) = g.N[1]
get_M(g::Geno) = g.M[1]
get_qced_genodir(g::Geno) = g.qced_genodir[1]
get_ldmafmat(g::Geno) = deepcopy(g.ldmafmat[1])
get_ldtype(g::Geno) = g.ldtype[1]
get_ldcomputed(g::Geno) = g.ldcomputed[1]
get_plink_binpath(g::Geno) = g.plink_binpath
get_chrs(g::Geno) = g.chrs
get_objpath(g::Geno) = isassigned(g.objpath, 1) ? g.objpath[1] : nothing

function check_merge_condition(gs::Vector{Geno})
    raw_genodirs = [g.raw_genodir for g in gs]
    raw_genobasenames = [g.raw_genobasename for g in gs]
    genos = [g.geno for g in gs]
    mafs = [g.maf for g in gs]
    ldtypes = [g.ldtype[1] for g in gs]
    imputeds = [g.imputed for g in gs]
    ldcomputeds = [g.ldcomputed[1] for g in gs]
    Ns = [g.N[1] for g in gs]
    @assert length(unique(raw_genodirs)) == 1 "Raw genotype file directories are not the same."
    @assert length(unique(raw_genobasenames)) == 1 "Raw genotype file base names are not the same."
    @assert length(unique(genos)) == 1 "Maximum missing genotype rates are not the same."
    @assert length(unique(mafs)) == 1 "Minimum MAFs are not the same."
    @assert length(unique(ldtypes)) == 1 "ldtypes are not the same."
    @assert length(unique(imputeds)) == 1 "Imputation statuses are not the same."
    @assert length(unique(ldcomputeds)) == 1 "LD computation statuses are not the same."
    @assert length(unique(Ns)) == 1 "The number of subjects are not the same."
    true
end

function update_qced_genodir_merged!(g_merged::Geno, gs::Vector{Geno})
    mkpath(g_merged.qced_genodir_parent)
    if g_merged.imputed
        merged_genodir = "$(g_merged.qced_genodir_parent)/qced_imputed_geno_chrs_$(g_merged.chrs[1])_to_$(g_merged.chrs[end])\
        _N_$(g_merged.N[1])_M_$(g_merged.M[1])_ldtype_$(g_merged.ldtype[1])"
    else
        merged_genodir = "$(g_merged.qced_genodir_parent)/qced_geno_chrs_$(g_merged.chrs[1])_to_$(g_merged.chrs[end])\
        _N_$(g_merged.N[1])_M_$(g_merged.M[1])_ldtype_$(g_merged.ldtype[1])"
    end
    mkpath(merged_genodir)
    for g in gs
        for file_path in readdir(g.qced_genodir[1], join=true)
            mv(file_path, "$(merged_genodir)/$(basename(file_path))", force = true)
        end
        rm(g.qced_genodir[1])
    end
    g_merged.qced_genodir[1] = merged_genodir
end