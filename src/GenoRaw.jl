export simulate_geno, partition_geno, simulate_bim_file, simulate_fam_file

function simulate_geno(
    rng::Random.AbstractRNG, 
    N::Int, 
    M::Int,
    K::Int, 
    rho::Real,
    mafs::Vector{<:AbstractFloat},
    genodir::AbstractString;
    save_unpartitioned::Bool = false,
    compress_partitioned::Bool = false,
    compress_unpartitioned::Bool = true,
    )
    @assert M == size(mafs, 1) "The length of MAFs should be equal to M"
    mkpath(genodir)
    open("$genodir/MAFs_true.txt", "w+") do stream
        writedlm(stream, mafs);
    end
    G = SnpArrays.simulate(rng, N, M, mafs, rho);
    partition_geno(G, K, genodir, mafs, compress=compress_partitioned);
    bim_array = Vector{Matrix{Any}}(undef, K)
    for k in 1:K
        bim_array[k] = readdlm("$genodir/G$k.bim", header=false)
    end
    if save_unpartitioned
        SnpArray("$genodir/G.bed", G)
        bim = reduce(vcat, bim_array)
        fam = readdlm("$genodir/G1.fam", header=false)
        open("$genodir/G.bim", "w+") do io
            writedlm(io, bim)
        end
        open("$genodir/G.fam", "w+") do io
            writedlm(io, fam)
        end
        if compress_unpartitioned
            compress_plink("$genodir/G")
            rm("$genodir/G.bed")
        end
    end
end

function partition_geno(
    G::SnpArray, 
    K::Int, 
    genodir::AbstractString, 
    MAFs::AbstractArray; 
    compress::Bool=false
    )
    N, M = size(G)
    function get_start_end(k::Int, K::Int, M::Int)
        s = 1 + (k - 1) * (M ÷ K)
        e = k == K ? M : k * (M ÷ K)
        s, e
    end
    for k in 1:K
        s, e = get_start_end(k, K, M)
        Mₖ = e - s + 1
        Gₖ = SnpArray(undef, N, Mₖ)
        copyto!(Gₖ.data, @view(G.data[:, s:e]))
        SnpArray("$genodir/G$k.bed", Gₖ) # save `geno` as a bed file
        simulate_bim_file(genodir, s, e, k)
        simulate_fam_file(genodir, N, k)
        if compress
            compress_plink("$genodir/G$k")
            rm("$datapath/G$k.bed")
        end
        open("$genodir/MAFs$(k)_true.txt", "w+") do io
            writedlm(io, MAFs[s:e])
        end
    end
end

function simulate_bim_file(genodir::AbstractString, s::Int, e::Int, k::Int)
    Mₖ = e - s + 1
    tab = DataFrame(chr=k, rsid="rs" .* string.(collect(s:e)), cm=repeat([0],Mₖ), bp=s:e, ALT=repeat(["A"],Mₖ), REF=repeat(["a"],Mₖ))
    open("$genodir/G$k.bim", "w+") do io
        writedlm(io, eachrow(tab))
    end
end

function simulate_fam_file(genodir::AbstractString, N::Int, k::Int)
    tab = DataFrame(FID=string.(1:N), IID=string.(1:N), IIDF=repeat([0],N), IIDM=repeat([0],N), SEX=repeat([0],N), PHENO=repeat([-9],N))
    open("$genodir/G$k.fam", "w+") do io
        writedlm(io, eachrow(tab))
    end
end