export simulate_pheno!

function simulate_pheno!(
    s::CVCDataset,
    cr::AbstractFloat,
    rep::Int,
    phenodir_parent::AbstractString
)
    sa_array = Array{SnpArray, 1}(undef, s.K);
    T = eltype(s.ηw)
    sla_array = Array{SnpLinAlg{T}, 1}(undef, s.K)
    construct_sa_array!(sa_array, s.part_genodir);
    construct_sla_array!(sla_array, sa_array);
    η = Vector{T}(undef, s.N);
    s.cr[] = cr; s.rep[] = rep;
    s.phenodir_parent[] = phenodir_parent;
    s.phenodir[] = "$(phenodir_parent)/$(get_phenodirname(s))"
    mkpath(s.phenodir[])
    Random.seed!(rep)
    update_mean_component!(η, s.ηw, sa_array, sla_array, s.vc_array);
    u, delta =
        simulate_right_censored_data(
            length(η);
            Y_shift = η,
            Y_scale = sqrt(s.phie),
            cr = cr
        );
    save_pheno(u, delta, s.phenodir[], rep)
end

function get_phenodirname(s::CVCDataset)
    "pheno_N_$(s.N)_M_$(s.M)_ldtype_$(s.ldtype)_rho_$(s.rho)_\
    C_$(s.C[])_K_$(s.K)_nldbins_$(s.nldbins)_nmafbins_$(s.nmafbins)_\
    h2_$(s.h2)_cr_$(s.cr[])_maflb_$(s.maflb)_mafub_$(s.mafub)_cvr_$(s.cvr)_a_$(s.a)_b_$(s.b)"
end

function save_pheno(u::Vector, delta::Vector, phenodir::AbstractString, rep::Int)
    open("$phenodir/u$rep.txt", "w+") do io
        writedlm(io, u)
    end    
    open("$phenodir/delta$rep.txt", "w+") do io
        writedlm(io, Int.(delta))
    end
end