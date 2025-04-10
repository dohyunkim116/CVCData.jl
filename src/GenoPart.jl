export 
    GenoPart,
    get_part_genodir,
    get_ldmafmat,
    partition_by_ldmaf!,
    merge_across_chr!,
    _merge_across_chr!,
    check_partition,
    convert_bin_names!
    
struct GenoPart
    N::Int
    M::Int
    K::Vector
    rho::AbstractFloat
    chrs::Union{Vector{Int}, UnitRange{Int}, Int}
    mafknots::Vector
    nmafbins::Vector
    nldbins::Int
    qced_genodir::AbstractString
    part_genodir_parent::AbstractString
    part_genodir::Vector
    ldmafmat::DataFrame
    ldtype::Symbol
    plink_binpath::AbstractString
    part_genoobj_dir::AbstractString
    conditional_ldbin::Bool
    conditional_mafbin::Bool
    nmafbins_per_interval::Int
    imputed::Bool
    isannot::Bool
end

function GenoPart(
    g::Geno,
    part_genodir_parent::AbstractString,
    part_genoobj_dir::AbstractString;
    mafknots::Vector=[],
    nmafbins::Int=4,
    nldbins::Int=4,
    conditional_ldbin::Bool=false,
    conditional_mafbin::Bool=false,
    nmafbins_per_interval::Int=9,
    isannot::Bool=false,
    annot_path::AbstractString=""
    )
    mkpath(part_genodir_parent)
    N = get_N(g)
    M = get_M(g)
    K = Vector{Int}(undef, 1)
    nmafbins = isempty(mafknots) ? nmafbins : length(mafknots) + 1
    ldmafmat = get_ldmafmat(g)
    ldtype = get_ldtype(g)
    qced_genodir = get_qced_genodir(g)
    chrs = get_chrs(g)
    part_genodir = mktempdir(part_genodir_parent; cleanup=false)
    plink_binpath = get_plink_binpath(g)
    rho = g.rho
    gp = GenoPart(N, M, K, rho, chrs, mafknots, [nmafbins], nldbins, qced_genodir, 
                part_genodir_parent, [part_genodir], ldmafmat, ldtype, plink_binpath, 
                part_genoobj_dir, conditional_ldbin, conditional_mafbin, nmafbins_per_interval,
                g.imputed, isannot)
    if !isannot
        update_ldmafmat_mafbin!(gp)
        update_ldmafmat_ldbin!(gp)
        update_ldmafmat_ldmafbin!(gp)
        update_K!(gp)
        update_nmafbins!(gp)
        update_nldbins!(gp)
    else
        merge_annot_to_ldmafmat!(gp, annot_path)
        update_ldmafmat_annotbin!(gp)
        update_K!(gp)
    end
    update_part_genodir!(gp)
    save_snpsets(gp)
    save_obj(gp)
    gp
end

function update_ldmafmat_mafbin!(gp::GenoPart)
    mafs = gp.ldmafmat[!, :maf];
    if isempty(gp.mafknots)
        gp.ldmafmat[!, :mafbin] = bin_maf(mafs, gp.nmafbins[1])
    else
        if gp.conditional_mafbin
            gp.ldmafmat[!, :mafbin] = bin_maf_conditional(mafs, gp.mafknots, gp.nmafbins_per_interval)
        else
            gp.ldmafmat[!, :mafbin] = bin_maf(mafs, gp.mafknots)    
        end
    end
end

bin_maf(mafs::Vector, nmafbins::Int)=bin_by_quantiles(mafs, nmafbins)

function bin_maf(mafs::Vector, mafknots::Vector)
    mafbins = Vector{Int}(undef, length(mafs));
    for i in eachindex(mafs)
        for j in 1:(length(mafknots) + 1) # j-th maf interval
            if j == 1
                if mafs[i] <= mafknots[j]
                    mafbins[i] = j
                end
            elseif j == length(mafknots) + 1
                if mafs[i] > mafknots[j - 1]
                    mafbins[i] = j
                end
            else
                if mafs[i] > mafknots[j - 1] && mafs[i] <= mafknots[j]
                    mafbins[i] = j
                end
            end
        end
    end
    mafbins
end

function bin_maf_conditional(mafs::Vector, mafknots::Vector, nmafbins_per_interval::Int)
    mafbins_by_knots = bin_maf(mafs, mafknots)
    mafbins = similar(mafbins_by_knots)
    for mafbin_by_knots in sort(unique(mafbins_by_knots))
        mafbins[mafbins_by_knots .== mafbin_by_knots] .= 
            bin_by_quantiles(mafs[mafbins_by_knots .== mafbin_by_knots], nmafbins_per_interval) .+ (mafbin_by_knots - 1) * nmafbins_per_interval;
    end
    mafbins
end

function update_ldmafmat_ldbin!(gp::GenoPart)
    lds = gp.ldmafmat[!, gp.ldtype];
    if gp.conditional_ldbin
        mafbins = gp.ldmafmat[!, :mafbin];
        ldbins = zeros(Int, size(gp.ldmafmat, 1));
        for mafbin in unique(mafbins)
            ldbins[mafbins .== mafbin] .= 
                bin_by_quantiles(lds[mafbins .== mafbin], gp.nldbins);
        end
        gp.ldmafmat[!, :ldbin] = ldbins;
    else
        gp.ldmafmat[!, :ldbin] = bin_by_quantiles(lds, gp.nldbins);
    end
end

function bin_by_quantiles(x::Vector, nbins::Int)
    function get_quantiles(x::Vector, n::Int)
        pvec = collect(0:(1/n):1)
        quantile(x, pvec);
    end
    quantiles = get_quantiles(x, nbins)
    bins = Vector{Int}(undef, length(x));
    for i in eachindex(x)
        for j in 1:length(quantiles)-1
            if j == 1
                if x[i] <= quantiles[j+1]
                    bins[i] = j
                end
            else
                if x[i] > quantiles[j] && x[i] <= quantiles[j+1]
                    bins[i] = j
                end
            end
        end
    end
    bins
end

function update_ldmafmat_ldmafbin!(gp::GenoPart)    
    ldbins = gp.ldmafmat[!, :ldbin];
    mafbins = gp.ldmafmat[!, :mafbin];
    gp.ldmafmat[!, :ldmafbin] = ["$(ldbins[i])" * "." * "$(mafbins[i])" for i in axes(gp.ldmafmat,1)]
end

function update_K!(gp::GenoPart)
    if !gp.isannot
        gp.K[1] = gp.ldmafmat[!,:ldmafbin] |> unique |> length
    else
        gp.K[1] = gp.ldmafmat[!,:annotbin] |> unique |> length
    end
end
function update_nmafbins!(gp::GenoPart)
    gp.nmafbins[1] = gp.ldmafmat[!, :mafbin] |> unique |> length
end
function update_nldbins!(gp::GenoPart)
    gp.nmafbins[1] = gp.ldmafmat[!, :ldbin] |> unique |> length
end

function update_part_genodir!(gp::GenoPart)
    chrstart = gp.chrs[1];
    chrend = gp.chrs[end];
    if !gp.isannot
        nldbins = get_nldbins(gp)
        nmafbins = get_nmafbins(gp)
        nmafbins_per_interval = get_nmafbins_per_interval(gp)
        if !gp.imputed
            newdir = 
                "$(gp.part_genodir_parent)/part_geno_chrs_$(chrstart)_to_$(chrend)\
                _N_$(gp.N)_M_$(gp.M)_K_$(gp.K[])_rho_$(gp.rho)_ldtype_$(gp.ldtype)_nldbins_$(nldbins)_nmafbins_$(nmafbins)\
                _nmafbins_per_interval_$(nmafbins_per_interval)"
        else
            newdir = 
                "$(gp.part_genodir_parent)/part_imputed_geno_chrs_$(chrstart)_to_$(chrend)\
                _N_$(gp.N)_M_$(gp.M)_K_$(gp.K[1])_rho_$(gp.rho)_ldtype_$(gp.ldtype)_nldbins_$(nldbins)_nmafbins_$(nmafbins)\
                _nmafbins_per_interval_$(nmafbins_per_interval)"
        end
        
    else
        if !gp.imputed
            newdir = 
                "$(gp.part_genodir_parent)/annot_part_geno_chrs_$(chrstart)_to_$(chrend)\
                _N_$(gp.N)_M_$(gp.M)_K_$(gp.K[1])_rho_$(gp.rho)"
        else
            newdir = 
                "$(gp.part_genodir_parent)/annot_part_imputed_geno_chrs_$(chrstart)_to_$(chrend)\
                _N_$(gp.N)_M_$(gp.M)_K_$(gp.K[1])_rho_$(gp.rho)"
        end
    end
    mv("$(gp.part_genodir[1])", newdir, force=true)
    gp.part_genodir[1] = newdir
end

function save_snpsets(gp::GenoPart)
    chrs_unique = unique(gp.ldmafmat[:,:chr]);
    part_genodir = get_part_genodir(gp);
    if !gp.isannot
        ldmafbins_unique = unique(gp.ldmafmat[:,:ldmafbin]);
        ldmat_array = Vector{DataFrame}(undef, length(chrs_unique) * length(ldmafbins_unique))
        i = 0
        # should be able to parallelize this; also you might be able to just use bin k to get the snpset
        for chr in chrs_unique, ldmafbin in ldmafbins_unique
            i += 1
            ldmat_array[i] = filter(row -> row.chr == chr && row.ldmafbin == ldmafbin, gp.ldmafmat)
            if size(ldmat_array[i], 1) > 0
                writedlm("$(part_genodir)/snpset_chr$(gp.ldtype)mafbin$(chr).$(ldmafbin).txt", ldmat_array[i][:,:snpid])
            end
        end
    else
        ldmat_array = Vector{DataFrame}(undef, length(chrs_unique) * gp.K[1])
        i = 0
        for chr in chrs_unique, k in 1:gp.K[1]
            i += 1
            ldmat_array[i] = filter(row -> row.chr == chr && row.annotbin == k, gp.ldmafmat)
            if size(ldmat_array[i], 1) > 0
                writedlm("$(part_genodir)/snpset_chr_$(chr)_annotbin_$(k).txt", ldmat_array[i][:,:snpid])
            end
        end
    end
end

function save_obj(gp::GenoPart)
    mkpath(gp.part_genoobj_dir)
    objpath = "$(gp.part_genoobj_dir)/$(get_objname(gp)).jld2"
    JLD2.save_object(objpath, gp)
end

function get_objname(gp::GenoPart)
    if !gp.isannot
        if gp.imputed
            return "part_imputed_geno_chrs_$(gp.chrs[1])_to_$(gp.chrs[end])_N_$(gp.N)_M_$(gp.M)_K_$(gp.K[1])_rho_$(gp.rho)_ldtype_$(gp.ldtype)"
        else
            return "part_geno_chrs_$(gp.chrs[1])_to_$(gp.chrs[end])_N_$(gp.N)_M_$(gp.M)_K_$(gp.K[1])_rho_$(gp.rho)_ldtype_$(gp.ldtype)"
        end
    else
        if gp.imputed
            return "annot_part_imputed_geno_chrs_$(gp.chrs[1])_to_$(gp.chrs[end])_N_$(gp.N)_M_$(gp.M)_rho_$(gp.rho)"
        else
            return "annot_part_geno_chrs_$(gp.chrs[1])_to_$(gp.chrs[end])_N_$(gp.N)_M_$(gp.M)_rho_$(gp.rho)"
        end
    end
end

function partition_by_ldmaf!(
    gp::GenoPart
)
    chrs_unique = unique(gp.ldmafmat[:, :chr])
    for chr in chrs_unique
        partition_by_ldmaf!(gp, chr)
    end
end

function partition_by_ldmaf!(gp::GenoPart, chr::Int)
    @assert !gp.isannot
    ldmafbins_unique = unique(gp.ldmafmat[:, :ldmafbin])
    for ldmafbin in ldmafbins_unique
        ldmafbin_path = "$(gp.part_genodir[1])/snpset_chr$(gp.ldtype)mafbin$(chr).$(ldmafbin).txt"
        bedbasename = "$(gp.part_genodir[1])/G_chr$(gp.ldtype)mafbin$(chr).$(ldmafbin)"
        if isfile(ldmafbin_path)
            cmd = `$(gp.plink_binpath)/plink2 \
                    --bfile $(gp.qced_genodir)/G$(chr) \
                    --extract $ldmafbin_path \
                    --make-bed \
                    --out $bedbasename`
            run(cmd)
        end
        rm(ldmafbin_path; force = true)
    end
end

function partition_by_annot!(
    gp::GenoPart
)
    @assert gp.isannot
    for chr in gp.chrs
        partition_by_annot!(gp, chr)
    end
end

function partition_by_annot!(
    gp::GenoPart,
    chr::Int
)
    @assert gp.isannot
    for k in gp.K[1]
        snpset_path = "$(gp.part_genodir[1])/snpset_chr_$(chr)_annotbin_$(k).txt"
        bedbasename = "$(gp.part_genodir[1])/G_chr_$(chr)_annotbin_$(k)"
        if isfile(snpset_path)
            cmd = `$(gp.plink_binpath)/plink2 \
                    --bfile $(gp.qced_genodir)/G$(chr) \
                    --extract $snpset_path \
                    --make-bed \
                    --out $bedbasename`
            run(cmd)
        end
        rm(snpset_path; force = true)
    end
end

function merge_across_chr!(
    gp::GenoPart;
    genobasename::AbstractString="G"
)
    if !gp.isannot
        ldmafbins_unique = unique(gp.ldmafmat[:, :ldmafbin])
        for ldmafbin in ldmafbins_unique
            _merge_across_chr!(gp; genobasename = genobasename, ldmafbin=ldmafbin)
        end
        check_partition(gp) ? convert_bin_names!(gp, genobasename) : error("Partititoning is incorrect.")
    else
        for annotbin in 1:gp.K[1]
            _merge_across_chr!(gp; genobasename = genobasename, annotbin=annotbin)
        end
        check_partition(gp)
    end
end

function _merge_across_chr!(
    gp::GenoPart;
    genobasename::AbstractString="G",
    ldmafbin::AbstractString="0",
    annotbin::Int=0    
)
    if !gp.isannot
        mergelistpath = "$(gp.part_genodir[1])/mergelist_$(gp.ldtype)mafbin$(ldmafbin).txt"
        rm(mergelistpath; force=true)
        outbedbasepath = "$(gp.part_genodir[1])/$(genobasename)$(ldmafbin)"
        chrs_unique = unique(gp.ldmafmat[:, :chr])
        # save mergelist
        for chr in chrs_unique
            bedbasepath = "$(gp.part_genodir[1])/G_chr$(gp.ldtype)mafbin$(chr).$(ldmafbin)"
            if ispath("$(bedbasepath).bed") 
                if ispath(mergelistpath)
                    open(mergelistpath, "a") do io
                        println(io, bedbasepath)
                    end
                else
                    open(mergelistpath, "w") do io
                        println(io, bedbasepath)
                    end
                end
            end
        end
    else
        mergelistpath = "$(gp.part_genodir[1])/mergelist_annotbin$(annotbin).txt"
        rm(mergelistpath; force=true)
        outbedbasepath = "$(gp.part_genodir[1])/$(genobasename)$(annotbin)"
        chrs_unique = unique(gp.ldmafmat[:, :chr])
        # save mergelist
        for chr in chrs_unique
            bedbasepath = "$(gp.part_genodir[1])/G_chr_$(chr)_annotbin_$(k)"
            if ispath("$(bedbasepath).bed") 
                if ispath(mergelistpath)
                    open(mergelistpath, "a") do io
                        println(io, bedbasepath)
                    end
                else
                    open(mergelistpath, "w") do io
                        println(io, bedbasepath)
                    end
                end
            end
        end
    end
    if countlines(mergelistpath) > 1
        cmd = `$(gp.plink_binpath)/plink2 \
            --pmerge-list $mergelistpath bfile \
            --make-bed \
            --out $outbedbasepath`; run(cmd);

    elseif countlines(mergelistpath) == 1
        cmd = `$(gp.plink_binpath)/plink2 \
            --bfile $(readlines(mergelistpath)[1]) \
            --make-bed \
            --out $outbedbasepath`; run(cmd);
    else
        error("No bed files found for ldmafbin: $ldmafbin");
    end
    rm(outbedbasepath * ".pgen"; force=true)
    rm(outbedbasepath * ".psam"; force=true)
    rm(outbedbasepath * ".pvar"; force=true)
    rm(outbedbasepath * ".log"; force=true)
    for bedbasepath in eachline(mergelistpath)
        rm(bedbasepath * ".bed"; force=true)
        rm(bedbasepath * ".bim"; force=true)
        rm(bedbasepath * ".fam"; force=true)
        rm(bedbasepath * ".log"; force=true)
    end
    rm(mergelistpath; force=true)
end

function clean_part_genodir(gp::GenoPart)
    regex_text = r".pgen|.psam|.pvar|.log|mergelist|G_chrldakmafbin|G_chrldscmafbin"
    rm.(filter(f -> occursin(regex_text, f), readdir(gp.part_genodir[1], join=true)))
end

function check_partition(gp::GenoPart)
    ldmafbins_unique = gp.ldmafmat[:,:ldmafbin] |> unique
    bim_array = Array{DataFrame}(undef, length(ldmafbins_unique))
    for (i, ldmafbin) in enumerate(ldmafbins_unique)
        bim_array[i] = DataFrame(readdlm("$(gp.part_genodir[1])/G$ldmafbin.bim", header = false), :auto)
        rename!(bim_array[i], :x1 => :chr, :x2 => :snpid)
        bim_array[i][!,:ldmafbin2] .= ldmafbin
        select!(bim_array[i], [:chr, :snpid, :ldmafbin2])
    end
    ldbim = reduce(vcat,bim_array)
    ldbim[!,:chr] = convert(Vector{Int}, ldbim[!,:chr])
    ldbim[!,:snpid] = convert(Vector{String}, ldbim[!,:snpid])
    ldbim[!,:ldmafbin2] = convert(Vector{String}, ldbim[!,:ldmafbin2])
    ldbim = leftjoin(ldbim, gp.ldmafmat[:, [:chr, :snpid, :ldmafbin]], on = [:chr, :snpid])
    all(ldbim[:,:ldmafbin] .== ldbim[:,:ldmafbin2])
end

function convert_bin_names!(
    gp::GenoPart,
    genobasename::String="G",
    )
    ldmafbins_unique = gp.ldmafmat[:,:ldmafbin] |> unique
    dict_ldmafbin_to_bin = Dict(zip(sort(ldmafbins_unique), collect(1:length(ldmafbins_unique))));
    part_genodir = get_part_genodir(gp)
    for g in split.(filter(contains(r"G\d+.\d+.bed"), readdir(part_genodir)), ".bed")
        g = g[1]
        ldmafbin = split(g, genobasename)[2]
        bin = dict_ldmafbin_to_bin[ldmafbin]
        mv("$(part_genodir)/$g.bed", "$(part_genodir)/$genobasename$bin.bed", force=true)
        mv("$(part_genodir)/$g.bim", "$(part_genodir)/$genobasename$bin.bim", force=true)
        mv("$(part_genodir)/$g.fam", "$(part_genodir)/$genobasename$bin.fam", force=true)
    end
    ldmafbin_map_tab = hcat(collect(keys(dict_ldmafbin_to_bin)), Int.(collect(values(dict_ldmafbin_to_bin))))
    open("$(part_genodir)/ldmafbin_bin_map", "w+") do stream
        writedlm(stream, ldmafbin_map_tab)
    end
    ldmafbin_map_tab = DataFrame(ldmafbin_map_tab, [:ldmafbin, :bin])
    leftjoin!(gp.ldmafmat, ldmafbin_map_tab, on = :ldmafbin)
end

# getters
get_N(gp::GenoPart) = gp.N
get_M(gp::GenoPart) = gp.M
get_K(gp::GenoPart) = gp.K[1]
get_part_genodir(gp::GenoPart) = gp.part_genodir[1]
get_nldbins(gp::GenoPart) = gp.nldbins
get_nmafbins(gp::GenoPart) = gp.nmafbins[1]
get_ldtype(gp::GenoPart) = gp.ldtype
get_ldmafmat(gp::GenoPart) = gp.ldmafmat
get_sla_array(gp::GenoPart) = gp.sla_array
get_conditional_ldbin(gp::GenoPart) = gp.conditional_ldbin
get_conditional_mafbin(gp::GenoPart) = gp.conditional_mafbin
get_nmafbins_per_interval(gp::GenoPart) = gp.conditional_mafbin ? gp.nmafbins_per_interval : zero(Int)

function merge_geno(
    gp::GenoPart,
    genobasename::AbstractString="G"
)   
    part_genodir = get_part_genodir(gp)
    bedfnames = filter(f -> occursin(r".bed", f), readdir(part_genodir, join=false));
    basenames = [split(bedfname, ".bed")[1] for bedfname in bedfnames];
    basepaths = part_genodir .* "/" .* basenames
    mergelistpath = "$(part_genodir)/mergelist.txt"
    rm(mergelistpath; force=true)
    open(mergelistpath, "w") do io
        writedlm(io, basepaths)
    end
    fampath = basepaths[1] * ".fam"
    indivsortdf = CSV.read(fampath, DataFrame; 
                            header = false, select = [1,2], types = String)
    indivsortpath = "$(part_genodir)/indivsort.txt"
    if(ispath(indivsortpath))
        rm(indivsortpath)
    end
    CSV.write(indivsortpath, indivsortdf, header=false, delim = "\t")
    cmd = `$(gp.plink_binpath)/plink \
        --merge-list $mergelistpath \
        --indiv-sort file $indivsortpath \
        --make-bed \
        --out $(part_genodir)/$genobasename`
    run(cmd);
    rm.(filter(f -> occursin(r".pgen|.psam|.pvar|.log|mergelist.txt|indivsort.txt|.nosex", f), 
        readdir(part_genodir, join=true)));
end

function generate_group_files(
    gp::GenoPart,
    genobasename::AbstractString="G",
    mS_per_group::String="0.0001,0.001,0.01"
    )
    ldmafbins_unique = unique(gp.ldmafmat[:,:ldmafbin]);
    dict_ldmafbin_to_bin = Dict(zip(sort(ldmafbins_unique), collect(1:length(ldmafbins_unique))));
    part_genodir = get_part_genodir(gp)
    bim = DataFrame(readdlm("$(part_genodir)/$genobasename.bim", header = false), :auto)
    rename!(bim, :x1 => :chr, :x2 => :snpid); select!(bim, [:chr, :snpid])
    bim = leftjoin(bim, gp.ldmafmat[:, [:chr, :snpid, :ldmafbin, :mafbin, :ldbin]], on = [:chr, :snpid])
    ldmafbins_unique = sort(unique(bim[:,:ldmafbin]))
    grouppath = "$(part_genodir)/G.group"
    mSpath = "$(part_genodir)/G.mS"
    group = bim[:,:ldmafbin] .|> x -> dict_ldmafbin_to_bin[x] - 1
    mS = repeat(mS_per_group * ";", length(ldmafbins_unique) - 1) * mS_per_group
    bin_map_path = "$(part_genodir)/bin_ldmafbin_map"
    bin_map_tab = hcat(collect(values(dict_ldmafbin_to_bin)), collect(keys(dict_ldmafbin_to_bin)))
    open(grouppath, "w+") do stream
        writedlm(stream, group)
    end
    open(mSpath, "w+") do stream
        write(stream, mS)
    end
    open(bin_map_path, "w+") do stream
        writedlm(stream, bin_map_tab)
    end
end