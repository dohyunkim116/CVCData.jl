# CVCData

[![Build Status](https://github.com/dohyunkim116/CVCData.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/dohyunkim116/CVCData.jl/actions/workflows/CI.yml)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- [![Coverage Status](https://coveralls.io/repos/github/dohyunkim116/CVCData.jl/badge.svg?branch=main)](https://coveralls.io/github/dohyunkim116/CVCData.jl?branch=main)
[![Documentation Status](https://readthedocs.org/projects/cvcdata/badge/?version=latest)](https://cvcdata.readthedocs.io/en/latest/?badge=latest) -->

A Julia package for preparing partitioned genetic data with different genetic architectures for heritability analysis using the Censored Variance Component (CVC) method.

## Installation

```julia
using Pkg
Pkg.add("CVCData")
```
The package also requires PLINK2, GCTA and LDAK software executables, as the user will specify the paths of these executables when constructing `Geno` object.

## Usage

```julia
using CVCData, Random

# 1. Set up directories
simdatadir = "path/to/simdata"
raw_genodir = "$simdatadir/geno/raw_geno_N_10000_M_20000"

# 2. Simulate genotype data
rng = Random.MersenneTwister(2025)
N = 10_000  # number of individuals
M = 20_000  # number of SNPs
nchrs = 22  # number of chromosomes
rho = 0.04  # LD correlation
mafs = 0.01 .+ (0.5 .- 0.01) .* rand(M)  # MAFs between 0.01 and 0.5
simulate_geno(rng, N, M, nchrs, rho, mafs, raw_genodir)

# 3. Create and process Geno object
qced_genodir_parent = "$simdatadir/qced_geno"
qced_genoobj_dir = "$simdatadir/geno_obj"
for chr in 1:nchrs
    g = Geno(
        "/path/to/plink", 
        "/path/to/gcta", 
        "/path/to/ldak", 
        raw_genodir, 
        "G", # raw genotype file basename
        qced_genodir_parent, # output directory for QCed genotype
        qced_genoobj_dir, # output directory for Geno object
        chr,  # chromosome
        0.01,  # minimum missing genptype rate (above this, SNPs are removed)
        0.01,  # minimum MAF threshold (below this, SNPs are removed)
        0.04   # LD correlation
    );
    update_ldmafmat!(g, ldtype=:ldak) # update LD scores using LDAK (LDAK scores)
end

# 4. Combine Geno objects per chromosome
gs = Geno[];
for chr in 1:nchrs
    gobjpath = glob("qced_geno_chrs_$(chr)_to_$(chr)_*.jls", qced_genoobj_dir)[]
    g = open(deserialize, gobjpath);
    push!(gs, g)
end
g = Geno(gs, qced_genodir_parent, qced_genoobj_dir);

# 5. Create partitioned genotype
gp = GenoPart(g, 
    "$simdatadir/part_geno",
    "$simdatadir/part_obj",
    nldbins=4, # number of LD bins; quartiles
    mafknots=[0.01, 0.05, 0.1, 0.2, 0.5]; # 6 MAF bins from 5 MAF knots
)

# 6. Partition by LD and MAF
chrs_unique = unique(gp.ldmafmat[:, :chr])
for chr in chrs_unique
    partition_by_ldmaf!(gp, chr)
end

# 7. Merge across chromosomes
ldmafbins_unique = unique(gp.ldmafmat[:, :ldmafbin]);
for ldmafbin in ldmafbins_unique
    _merge_across_chr!(gp, ldmafbin = ldmafbin)
end
# 7a). Check if the partitioning is correct
check_partition(gp)
# 7b). Convert bin names to more readable format.
# Bin mapping will also be created for each ldmafbin 
# in the directory where the partitioned genotypes are stored
convert_bin_names!(gp) 

# 8. Simulate and align covariates
# N is number of individuals, C is number of covariates
covdir = "$simdatadir/cov/cov_N_5000_C_24" # output directory for covariates
simulate_cov(rng, 5000, 24, covdir)
covdir_aligned, partgenodir_aligned, gp_aligned = align!(gp, covdir)

# 9. Create genetic architecture
# The causal variant ratio is the proportion of causal variants in the total number of variants.
# The MAF lower and upper bounds are used to filter the variants for causal variant selection.
# The alpha and beta parameters are used to modulate the impact of LD and MAF on the selection of causal variants.
ga = GenoArch(gp_aligned, 
    "$simdatadir/genoarch_obj", 
    0.5,   # heritability
    0.01,  # MAF lower bound
    0.5,   # MAF upper bound
    0.01,  # causal variant ratio
    0.75,  # alpha parameter modulating the impact of LD
    0.0    # beta parameter modulating the impact of MAF
) 

# 10. Create dataset and simulate phenotype
datasetobj_outputdir = "$simdatadir/dataset_obj"
dataset = CVCDataset(ga, gp_aligned, covdir_aligned, datasetobj_outputdir)

# The dataset object contains the genotype, covariate, and phenotype data.
cr = 0.2; rep = 1; phenodir_parent = "$simdatadir/pheno"
simulate_pheno!(dataset, cr, rep, phenodir_parent)
```