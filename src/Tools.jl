function get_fampath(genodir::AbstractString; imputed::Bool=false)
    if imputed
        famfilename = filter(contains(r".sample"), readdir(genodir))[1]
    else
        famfilename = filter(contains(r".fam"), readdir(genodir))[1]
    end
    joinpath(genodir, famfilename)
end

function get_n(genodir::AbstractString; imputed::Bool=false)
    fampath = get_fampath(genodir; imputed=imputed)
    n = open(fampath, "r") do io
        countlines(io)
    end
    imputed ? n - 1 : n
end

function get_m(genodir::AbstractString, k::Int)
    bimpath = joinpath(genodir, "G$(k).bim")
    M_k = open(bimpath, "r") do io
        countlines(io)
    end
    M_k
end

get_k(genodir::AbstractString)=filter(contains(r"G\d+.fam"), readdir(genodir)) |> length

function construct_sa_array!(sa_array::Array{SnpArray,1}, genopath::AbstractString)
    N = get_n(genopath)
    K = get_k(genopath)
    @assert K == length(sa_array)
    for k in eachindex(sa_array)
        Gpath = isfile("$(genopath)/G$k.bed.gz") ? "$(genopath)/G$k.bed.gz" : "$(genopath)/G$k.bed"
        sa_array[k] = SnpArray(Gpath, N)
    end
end

function construct_sla_array!(sla_array::Array{SnpLinAlg{T}, 1}, sa_array::Array{SnpArray,1}) where T
    @assert length(sla_array) == length(sa_array)
    for k in eachindex(sa_array)
        Xₖ = SnpLinAlg{T}(sa_array[k], model = ADDITIVE_MODEL, center = true, scale = true)
        standardize!(Xₖ, sa_array[k])
        sla_array[k] = Xₖ
    end
end

function construct_sla_array!(sla_array::Array{SnpLinAlg{T}, 1}, genopath) where T
    N = get_n(genopath)
    K = get_k(genopath)
    @assert K == length(sla_array)
    for k in eachindex(sla_array)
        Gpath = isfile("$(genopath)/G$k.bed.gz") ? "$(genopath)/G$k.bed.gz" : "$(genopath)/G$k.bed"
        Gₖ = SnpArray(Gpath, N)
        Xₖ = SnpLinAlg{T}(Gₖ, model = ADDITIVE_MODEL, center = true, scale = true)
        standardize!(Xₖ,Gₖ)
        sla_array[k] = Xₖ
    end
end

function standardize!(X::SnpLinAlg, G::SnpArray)
    n, p = size(G)
    @inbounds @views for j in 1:p
        #X.μ[j] = (X.μ[j] * G.columncounts[2,j] + G.columncounts[3,j] + 2G.columncounts[4,j]) / n
        X.σinv[j] = (abs2(X.μ[j]) * G.columncounts[1, j] + 
            abs2(1.0 - X.μ[j]) * G.columncounts[3, j] + 
            abs2(2.0 - X.μ[j]) * G.columncounts[4,j]) / n
        X.σinv[j] = X.σinv[j] > 0 ? inv(sqrt(X.σinv[j])) : one(eltype(X.σinv))
    end
end