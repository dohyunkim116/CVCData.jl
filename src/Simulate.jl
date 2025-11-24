function simulate_right_censored_data(
    m           :: Integer;
    cr          :: Real = 0.2,
    Y_dist      :: Distributions.UnivariateDistribution = Normal(),
    C_dist      :: Distributions.UnivariateDistribution = Normal(),
    Y_shift     = 0.0,
    Y_scale     = 1.0
)
    y = Vector{Float64}(undef, m)
    rand!(Y_dist, y)
    # standardize y
    y .= (y .- mean(Y_dist)) .* inv(sqrt(var(Y_dist)))
    # generate event times
    y .= y .* Y_scale .+ Y_shift
    
    if(cr == 0) return y, ones(Bool, m) end
    
    σy = std(y)

    c₀ = Vector{Float64}(undef, m)
    # generate censoring times
    rand!(C_dist, c₀)
    # standardize c
    c₀ .= (c₀ .- mean(C_dist)) .* inv(sqrt(var(C_dist)))
    # f(μ) calculates the proportion of censoring if using μ in C = σy * c₀ + μ
    f(μ) = begin
        cnt = 0
        @inbounds @simd for i in 1:m
            if μ + σy * c₀[i] ≤ y[i]
                cnt += 1
            end
            #cnt += μ + σy * c₀[i] ≤ y[i] ? 1 : 0
        end
        cnt / m - cr
    end
    μc = find_zero(f, (minimum(y) - 3σy, maximum(y) + 3σy), Bisection())
    c  = μc .+ σy .* c₀
    # T̃ are right-censored times
    ỹ  = copy(y)
    isobserved = Vector{Bool}(undef, m)
    @inbounds @simd for i in 1:m
        isobserved[i] = y[i] ≤ c[i]
        if !isobserved[i]; ỹ[i] = c[i]; end
    end
    ỹ, isobserved
end

function simulate_gumbel(
    m           :: Integer;
    cr          :: Real = 0.2,
    Y_shift     = 0.0,
    Y_scale     = 1.0
)
    y = Vector{Float64}(undef, m)
    rand!(Gumbel(), y)
    # standardize y
    y .= (y .- mean(Gumbel())) .* inv(sqrt(var(Gumbel())))
    # generate event times
    y .= y .* Y_scale .+ Y_shift
    
    if(cr == 0) return y, ones(Bool, m) end
    
    σy = std(y)

    c₀ = Vector{Float64}(undef, m)
    # generate censoring times
    rand!(Normal(), c₀)
    # standardize c
    c₀ .= (c₀ .- mean(Normal())) .* inv(sqrt(var(Normal())))
    # f(μ) calculates the proportion of censoring if using μ in C = σy * c₀ + μ
    f(μ) = begin
        cnt = 0
        @inbounds @simd for i in 1:m
            if μ + σy * c₀[i] ≤ y[i]
                cnt += 1
            end
            #cnt += μ + σy * c₀[i] ≤ y[i] ? 1 : 0
        end
        cnt / m - cr
    end
    μc = find_zero(f, (minimum(y) - 3σy, maximum(y) + 3σy), Bisection())
    c  = μc .+ σy .* c₀
    # T̃ are right-censored times
    ỹ  = copy(y)
    isobserved = Vector{Bool}(undef, m)
    @inbounds @simd for i in 1:m
        isobserved[i] = y[i] ≤ c[i]
        if !isobserved[i]; ỹ[i] = c[i]; end
    end
    ỹ, isobserved
end

function simulate_logistic(
    m           :: Integer;
    cr          :: Real = 0.2,
    Y_shift     = 0.0,
    Y_scale     = 1.0
)
    y = Vector{Float64}(undef, m)
    rand!(Logistic(), y)
    # standardize y
    y .= (y .- mean(Logistic())) .* inv(sqrt(var(Logistic())))
    # generate event times
    y .= y .* Y_scale .+ Y_shift
    
    if(cr == 0) return y, ones(Bool, m) end
    
    σy = std(y)

    c₀ = Vector{Float64}(undef, m)
    # generate censoring times
    rand!(Normal(), c₀)
    # standardize c
    c₀ .= (c₀ .- mean(Normal())) .* inv(sqrt(var(Normal())))
    # f(μ) calculates the proportion of censoring if using μ in C = σy * c₀ + μ
    f(μ) = begin
        cnt = 0
        @inbounds @simd for i in 1:m
            if μ + σy * c₀[i] ≤ y[i]
                cnt += 1
            end
            #cnt += μ + σy * c₀[i] ≤ y[i] ? 1 : 0
        end
        cnt / m - cr
    end
    μc = find_zero(f, (minimum(y) - 3σy, maximum(y) + 3σy), Bisection())
    c  = μc .+ σy .* c₀
    # T̃ are right-censored times
    ỹ  = copy(y)
    isobserved = Vector{Bool}(undef, m)
    @inbounds @simd for i in 1:m
        isobserved[i] = y[i] ≤ c[i]
        if !isobserved[i]; ỹ[i] = c[i]; end
    end
    ỹ, isobserved
end

function simulate_mixture_normal(
    m           :: Integer;
    cr          :: Real = 0.2,
    Y_shift     = 0.0,
    Y_scale     = 1.0
)
    y = Vector{Float64}(undef, m)
    d = MixtureModel(Normal, [(-3.0, 1.0), (-2.0, 1.2), (0.0, 1.0), (3.0, 2.5), (5.0, 1.5)], [0.1, 0.1, 0.3, 0.4, 0.1])
    rand!(d, y)
    # standardize y
    y .= (y .- mean(d)) .* inv(sqrt(var(d)))
    # generate event times
    y .= y .* Y_scale .+ Y_shift
    
    if(cr == 0) return y, ones(Bool, m) end
    
    σy = std(y)

    c₀ = Vector{Float64}(undef, m)
    # generate censoring times
    rand!(Normal(), c₀)
    # standardize c
    c₀ .= (c₀ .- mean(Normal())) .* inv(sqrt(var(Normal())))
    # f(μ) calculates the proportion of censoring if using μ in C = σy * c₀ + μ
    f(μ) = begin
        cnt = 0
        @inbounds @simd for i in 1:m
            if μ + σy * c₀[i] ≤ y[i]
                cnt += 1
            end
            #cnt += μ + σy * c₀[i] ≤ y[i] ? 1 : 0
        end
        cnt / m - cr
    end
    μc = find_zero(f, (minimum(y) - 3σy, maximum(y) + 3σy), Bisection())
    c  = μc .+ σy .* c₀
    # T̃ are right-censored times
    ỹ  = copy(y)
    isobserved = Vector{Bool}(undef, m)
    @inbounds @simd for i in 1:m
        isobserved[i] = y[i] ≤ c[i]
        if !isobserved[i]; ỹ[i] = c[i]; end
    end
    ỹ, isobserved
end

function simulate_right_censored_data_random(
    m           :: Integer;
    cr          :: Real = 0.2,
    Y_dist      :: Distributions.UnivariateDistribution = Normal(),
    Y_shift     = 0.0,
    Y_scale     = 1.0
)
    y = Vector{Float64}(undef, m)
    rand!(Y_dist, y)
    # standardize y
    y .= (y .- mean(Y_dist)) .* inv(sqrt(var(Y_dist)))
    # generate event times
    y .= y .* Y_scale .+ Y_shift
    
    if(cr == 0) return y, ones(Bool, m) end
    
    σy = std(y)

    censored_indices = sample(1:m, Int(round(cr * m)), replace = false)

    c = Vector{Float64}(undef, length(censored_indices))
    c = y[censored_indices] .- rand(Exponential(1/σy), length(censored_indices))  # random censoring times less than event times
    isobserved = trues(m)
    isobserved[censored_indices] .= false
    ỹ  = copy(y)
    ỹ[censored_indices] .= c
    ỹ, isobserved
end

function update_mean_component!(
    η::Vector{T},
    ηw::Vector{T},
    sa_array::Array{SnpArray, 1},
    sla_array::Array{SnpLinAlg{T}, 1},
    vc_array::Array{Vector{T}, 1};
    rowids::Union{Vector{Int}, Nothing} = nothing,
    D::Distributions.UnivariateDistribution = Normal(),
    sla::Bool = false
) where T
    fill!(η, zero(T))
    _update_mean_component!(
        η, 
        sa_array, 
        sla_array, 
        vc_array; 
        rowids = rowids,
        D = D, 
        sla = sla
    )
    LinearAlgebra.axpby!(1.0, ηw, 1.0, η)
end

function _update_mean_component!(
    out::Vector{T},
    sa_array::Array{SnpArray, 1},
    sla_array::Array{SnpLinAlg{T}, 1},
    vc_array::Array{Vector{T}, 1};
    rowids::Union{Vector{Int}, Nothing} = nothing,
    D::Distributions.UnivariateDistribution = Normal(),
    sla::Bool = false
) where T
    @assert length(sla_array) == length(vc_array) "length of sla_array and vc_array must be equal."
    if rowids !== nothing && sla
        error("rowids option is only available when sla=false")
    end
    K = length(sla_array)
    ηₖ = Vector{T}(undef, length(out))
    for k in 1:K
        fill!(ηₖ, zero(T))
        βₖ = Vector{T}(undef, length(vc_array[k]))
        rand!(D, βₖ)
        βₖ .= (βₖ .- mean(D)) .* inv(sqrt(var(D)))
        βₖ .= βₖ .* sqrt.(vc_array[k])
        if sla
            mul!(ηₖ, sla_array[k], βₖ)
        else
            if size(sla_array[k], 1) * size(sla_array[k], 2) ≤ 5e7
                # Always create full matrix first
                full_Xₖ = Matrix{T}(undef, size(sla_array[k], 1), size(sla_array[k], 2))
                mapsnp!(full_Xₖ, sa_array[k].data, sla_array[k].μ, sla_array[k].σinv)
                
                # Then use the appropriate subset for multiplication
                if isnothing(rowids)
                    mul!(ηₖ, full_Xₖ, βₖ)
                else
                    mul!(ηₖ, full_Xₖ[rowids, :], βₖ)
                end
            else
                f = 1000
                if mod(size(sla_array[k], 2), f) == 0
                    q, r = divrem(size(sla_array[k], 2), f)
                    full_Xₖₗ = Matrix{T}(undef, size(sla_array[k], 1), q)
                    
                    for i in 0:(f-1)
                        s = i * q + 1
                        e = (s - 1) + q   
                        mapsnp!(
                            full_Xₖₗ,
                            @view(sa_array[k].data[:,s:e]),
                            @view(sla_array[k].μ[s:e]),
                            @view(sla_array[k].σinv[s:e])
                        )
                        
                        if isnothing(rowids)
                            mul!(ηₖ, full_Xₖₗ, @view(βₖ[s:e]), one(T), one(T))
                        else
                            mul!(ηₖ, full_Xₖₗ[rowids, :], @view(βₖ[s:e]), one(T), one(T))
                        end
                    end
                else
                    q, r = divrem(size(sla_array[k], 2), f - 1)
                    full_Xₖₗ = Matrix{T}(undef, size(sla_array[k], 1), q)
                    full_Xₖₗ2 = Matrix{T}(undef, size(sla_array[k], 1), r)
                    
                    for i in 0:(f-2)
                        s = i * q + 1
                        e = (s - 1) + q
                        CVCData.mapsnp!(
                            full_Xₖₗ,
                            @view(sa_array[k].data[:,s:e]),
                            @view(sla_array[k].μ[s:e]),
                            @view(sla_array[k].σinv[s:e])
                        )
                        
                        if isnothing(rowids)
                            mul!(ηₖ, full_Xₖₗ, @view(βₖ[s:e]), one(T), one(T))
                        else
                            mul!(ηₖ, full_Xₖₗ[rowids, :], @view(βₖ[s:e]), one(T), one(T))
                        end
                    end
                    
                    s = (f - 1) * q + 1
                    e = (s - 1) + r
                    CVCData.mapsnp!(
                        full_Xₖₗ2,
                        @view(sa_array[k].data[:,s:e]),
                        @view(sla_array[k].μ[s:e]),
                        @view(sla_array[k].σinv[s:e])
                    )
                    
                    if isnothing(rowids)
                        mul!(ηₖ, full_Xₖₗ2, @view(βₖ[s:e]), one(T), one(T))
                    else
                        mul!(ηₖ, full_Xₖₗ2[rowids, :], @view(βₖ[s:e]), one(T), one(T))
                    end
                end
            end
        end
        LinearAlgebra.axpby!(1.0, ηₖ, 1.0, out)
    end
end

# mapsnp! functions are authored by Brendon Chau
"""
    mapsnp!(f, g, μ, σinv)

In-place additive imputation of a compressed column vector of a SnpArray `g` 
into vector of floating point numbers `f`. Centering and scaling are performed
automatically using scalars `μ` and `σinv`.
"""
@inline function mapsnp!(
    f::AbstractVector{T}, 
    g::AbstractVector{UInt8}, 
    μ::T, 
    σinv::T
    ) where {T<:AbstractFloat}
    K, rem = divrem(length(f), 4)
    @inbounds for j in 1:K
        # Base.Cartesian.@nexprs 4 i -> begin
        #     snp_i = (g[j] >> ((i - 1) * 2)) & 0x03
        #     f[(j - 1) << 2 + i] = (T((snp_i >= 0x02) * (snp_i - 0x01)) - !isone(snp_i) * μ) * σinv
        # end
        @simd for i in 1:4
            # rotate the correct value into the 0 and 1 place
            snp = (g[j] >> ((i - 1) << 1)) & 0x03
            # additive model imputation in-place
            f[(j - 1) << 2 + i] = (T((snp >= 0x02) * (snp - 0x01)) - !isone(snp) * μ) * σinv
        end
    end
    @inbounds if rem > 0
        @simd for i in 1:rem
            snp = (g[K + 1] >> ((i - 1) << 1)) & 0x03
            # additive model imputation in-place
            f[K << 2 + i] = (T((snp >= 0x02) * (snp - 0x01)) - !isone(snp) * μ) * σinv
        end
    end
end

"""
    mapsnp!(F, G, μ, σinv)

In-place additive imputation of a compressed SnpArray `G` into a matrix of 
floating point numbers `F`, parallelized over the columns of `F` and `G`. 
Centering and scaling are performed automatically using vectors `μ` and `σinv`.
"""
function mapsnp!(
    F::AbstractMatrix{T}, 
    G::AbstractMatrix{UInt8}, 
    μ::AbstractVector{T}, 
    σinv::AbstractVector{T}
    ) where {T<:AbstractFloat}
    @inbounds Threads.@threads for j in axes(F, 2)
        mapsnp!(view(F, :, j), view(G, :, j), μ[j], σinv[j])
    end
    return F
end

function mapsnp!(
    F::AbstractMatrix{T}, 
    G::AbstractSnpArray, 
    μ::AbstractVector{T}, 
    σinv::AbstractVector{T}
    ) where {T<:AbstractFloat}
    mapsnp!(F, G.data, μ, σinv)
    return F
end
