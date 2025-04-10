export GenoArch

struct GenoArch
    N::Int
    M::Int
    rho::Real
    qced_genodir::AbstractString
    ldmafmat::DataFrame
    ldtype::Symbol
    h2::Real
    maflb::Real
    mafub::Real
    cvr::Real
    a::Real
    b::Real
    h2hat::Real
    cvrhat::Real
    normalize::Bool
    phig::Real
    phie::Real
    rng::AbstractRNG
end

function GenoArch(
    g::Geno,
    h2::Real,
    maflb::Real,
    mafub::Real,
    cvr::Real,
    a::Real,
    b::Real;
    normalize::Bool=false,
    rng::AbstractRNG=Random.MerseeneTwister(2024)
    )
    ldmafmat = get_ldmafmat(g)
    ldtype = g.ldtype[]
    update_ldmafmat_vc!(
        rng, ldmafmat, h2, maflb, mafub, cvr, a, b; 
        ldtype = ldtype, normalize = normalize
    )
    phig = get_phig(ldmafmat)
    phie = get_phie(h2, ldmafmat, normalize)
    h2hat = phig / (phig + phie)
    cvrhat = get_numcv(ldmafmat) / g.M[1]
    GenoArch(g.N[1], g.M[1], g.rho, g.qced_genodir[1], ldmafmat,
            ldtype, h2, maflb, mafub, cvr, a, b, h2hat, 
            cvrhat, normalize, phig, phie, rng
    )
end

get_N(ga::GenoArch) = ga.N
get_M(ga::GenoArch) = ga.M
get_rho(ga::GenoArch) = ga.rho
get_qced_genodir(ga::GenoArch) = ga.qced_genodir
get_ldmafmat(ga::GenoArch) = deepcopy(ga.ldmafmat)
get_ldtype(ga::GenoArch) = ga.ldtype

function update_ldmafmat_vc!(
    rng::Random.AbstractRNG,
    ldmafmat::DataFrame,
    h2::Real,
    maflb::Real,
    mafub::Real,
    cvr::Real,
    a::Real,
    b::Real;
    ldtype::Symbol=:ldak,
    normalize::Bool=false    
    )
    rowmask = ldmafmat[:,:maf] .>= maflb .&& ldmafmat[:,:maf] .< mafub
    candcvsnprowids = ldmafmat[rowmask,:][:,:rowid]
    cvnum = Int(round(cvr * length(rowmask)))
    cvrowids = sample(rng, candcvsnprowids, cvnum, replace = false)
    cvrowmask = in.(ldmafmat[:,:rowid], Ref(cvrowids))
    vc = compute_vc.(
            cvrowmask, 
            ldmafmat[:,ldtype], 
            ldmafmat[:,:maf], 
            Ref(a), 
            Ref(b); 
            ldtype = ldtype
        )
    ldmafmat[!,:vc] = normalize ? normalize_vc(vc, h2) : vc
    ldmafmat[!,:isnormalized] .= normalize
end

function compute_vc(
    c::Union{Int,Bool}, 
    ld::Real, 
    f::Real, 
    a::Real, 
    b::Real;
    ldtype::Symbol=:ldak
    )
    if ldtype == :ldak
        w = ld
    elseif ldtype == :ldsc
        w = inv(ld)
    else
        error("Invalid ldtype")
    end
    c * w^b * (f * (1 - f))^a
end

function normalize_vc(vc::Vector, h2::Real)
    s = h2 / sum(vc)
    s .* vc
end

get_phig(ldmafmat::DataFrame)=sum(ldmafmat[:,:vc])
get_phie(h2, ldmafmat::DataFrame, normalize::Bool)=normalize ? 1 - h2 : (1 / h2 - 1) * get_phig(ldmafmat)
get_numcv(ldmafmat::DataFrame)=sum(ldmafmat[:,:vc] .> 0)
get_h2(ga::GenoArch)=ga.h2
get_h2hat(ga::GenoArch)=ga.h2hat
get_maflb(ga::GenoArch)=ga.maflb
get_mafub(ga::GenoArch)=ga.mafub
get_cvr(ga::GenoArch)=ga.cvr
get_cvrhat(ga::GenoArch)=ga.cvrhat
get_a(ga::GenoArch)=ga.a
get_b(ga::GenoArch)=ga.b
get_vcmat(ga::GenoArch)=ga.ldmafmat[:,[:chr, :snpid, :vc]]
get_normalize(ga::GenoArch)=ga.normalize
get_phig(ga::GenoArch)=ga.phig
get_phie(ga::GenoArch)=ga.phie

function Base.show(io::IO, ga::GenoArch)
    output = "N=$(ga.N), M=$(ga.M), h2=$(ga.h2), \
    ldtype=$(ga.ldtype), maflb=$(ga.maflb), mafub=$(ga.mafub), \
    cvr=$(ga.cvr), a=$(ga.a), b=$(ga.b), normalize=$(ga.normalize), \
    phig=$(ga.phig), phie=$(ga.phie), \
    h2hat=$(ga.h2hat), cvrhat=$(ga.cvrhat)"
    println("Showing genetic architecture setting:")
    println_wrapped(output, width=50)
    println(io)
end