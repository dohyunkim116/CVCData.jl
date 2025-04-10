export simulate_cov

function simulate_cov(
    rng::Random.AbstractRNG, 
    N::Int,
    C::Int,
    covdir::AbstractString,
    covfile_basename::AbstractString="w"
    )
    Random.seed!(rng)
    mkpath(covdir)
    covfile_df = covfile_basename * "df.txt"
    w = [string.(1:N) ones(N) randn(N, C - 1)];
    wdf = DataFrame(w, ["f.eid"; "intercept"; string.(1:C-1)])
    CSV.write("$covdir/$(covfile_df)", wdf; delim = "\t", quotestrings = true, header=true)
    
    covfile = covfile_basename * ".txt"
    CSV.write("$(covdir)/$(covfile)", wdf[:,Not(:"f.eid")]; delim = "\t", header=false)
end