# phyloSMC
phyloSMC samples phylogenetic trees with annealed and online sequential Monte Carlo. Curupixa supplies the likelihood, tree priors, molecular clocks, tree operations and numerical weight calculations. This repository contains the particle population, temperature schedule, importance proposals and Metropolis–Hastings rejuvenation.
The initial implementation provides a C API and a small DNA command-line program. Annealed inference supports non-clock trees, strict clocks, independent lognormal or exponential relaxed clocks, and a branch-midpoint autocorrelated lognormal clock. Online inference adds one taxon at a time to a non-clock posterior sample, allowing existing relationships to change during rejuvenation. Dated online proposals are not implemented.
## Build
Clone the repository and its pinned Curupixa dependency over SSH:
```sh
git clone --recurse-submodules git@github.com:leomrtns/phyloSMC.git
cmake -S phyloSMC -B phyloSMC/build/release -DCMAKE_BUILD_TYPE=Release
cmake --build phyloSMC/build/release --parallel
ctest --test-dir phyloSMC/build/release --output-on-failure
```
A C17 compiler and CMake 3.16 or later are required. Nothing is downloaded during configuration. Add `-DBUILD_SHARED_LIBS=ON` to build shared libraries for CFFI or ctypes. During joint development, `-DPHYLOSMC_CURUPIXA_SOURCE=/path/to/curupixa` selects another checkout.
## Run
```sh
build/release/phyloSMC alignment.fasta --particles 512 --moves 16 --seed 42
build/release/phyloSMC alignment.fasta --particles 512 --clock lognormal
build/release/phyloSMC alignment.fasta --particles 512 --initial-taxa 4
```
Input is aligned DNA in FASTA or the NEXUS subset accepted by Curupixa. The command-line program uses JC69 with fixed substitution parameters. Non-clock trees have a rooted PDA topology prior and independent exponential branch lengths with mean 0.1 substitutions/site. Clock runs use contemporary tips, a constant-size coalescent with pairwise scale 1, and mean substitution rate 1 per time unit. Lognormal clocks use log standard deviation 0.3; the autocorrelated model uses the corresponding diffusion parameter. These defaults define an example analysis; use the C API to supply other model parameters and dated tree priors.
Standard output contains a weight and a Newick tree per particle, with branch lengths in substitutions/site. Standard error reports the evidence estimate, final ESS and move acceptance. Particle weights may differ; retain them when computing summaries. Clock ages and branch rates are available through the C API. The output is a posterior sample, and no consensus or best tree is selected automatically.
## Use the C API
Include `phylosmc.h`, initialise `psmc_options`, then call `psmc_anneal()` with Curupixa character data, a substitution model and a random callback. Pass the result and extended data to `psmc_online()` to add a taxon. Release results with `psmc_result_free()`. Results retain the observations; keep the substitution model alive until all results using it have been freed. See the [worked example](examples/api.c) and [algorithm definitions](docs/algorithms.md).
Embed this repository with `add_subdirectory()` and link `phyloSMC::phylosmc`. An installed copy supports `find_package(phyloSMC CONFIG REQUIRED)`. The executable and library start no worker threads. Independent analyses may run concurrently with separate random streams; your program controls parallelism.
## Scope and validation
This is an initial reference implementation: proposals use full tree copies and likelihood evaluations. Incremental partial updates, shared particle caches, guided online placement and serially dated online inference remain future work. Substitution and clock hyperparameters are supplied rather than fitted. Population-size and diversification-rate optimisation are also outside the initial engine. Resampling can reduce particle diversity; ESS alone does not establish adequate exploration, so compare independent runs and increase both particle count and rejuvenation when necessary.
Tests check analytic marginal likelihoods, prior branch means and topology probabilities, online importance weighting, clock proposals, reproducibility and invalid inputs. See [validation](docs/validation.md). Implementation follows [Wang et al.](https://doi.org/10.1093/sysbio/syz028) for annealing and [Fourment et al.](https://doi.org/10.1093/sysbio/syx090) for online SMC principles.
Copyright (C) 2026 Leonardo de Oliveira Martins. Licensed under [GPL-3.0-or-later](LICENSE).
