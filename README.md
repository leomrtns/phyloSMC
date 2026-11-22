# phyloSMC
phyloSMC samples phylogenetic trees and evolutionary parameters using sequential Monte Carlo (SMC). Annealed inference starts from a joint prior and gradually introduces the sequence likelihood. Online inference adds taxa to an existing posterior population, then allows the enlarged trees and their parameters to change during rejuvenation.
The command-line program reads aligned DNA. The C API also accepts other finite-state observations and supplied substitution models. Curupixa provides tree sampling and operations, substitution likelihoods, clock distributions and numerical weight calculations. phyloSMC manages particles, importance updates and Metropolis–Hastings moves.
For a mathematical account, including derivations of the likelihood, priors, importance weights and proposal corrections, read the [tutorial](docs/tutorial.md). The [validation guide](docs/validation.md) describes independent checks of the calculations.
## Build
Clone over SSH, including the pinned Curupixa dependency:
```sh
git clone --recurse-submodules git@github.com:leomrtns/phyloSMC.git
cmake -S phyloSMC -B phyloSMC/build/release -DCMAKE_BUILD_TYPE=Release
cmake --build phyloSMC/build/release --parallel
ctest --test-dir phyloSMC/build/release --output-on-failure
```
A C17 compiler and CMake 3.16 or later are required. Configuration downloads nothing. Add `-DBUILD_SHARED_LIBS=ON` for shared libraries. During joint development, configure with `-DPHYLOSMC_CURUPIXA_SOURCE=/path/to/curupixa` to use another checkout; that path stays in the build configuration.
## An ordinary unrooted analysis
```sh
build/release/phyloSMC alignment.fasta --particles 512 --moves 16 --seed 42
```
This runs unrooted, non-clock inference under HKY, jointly sampling nucleotide frequencies, the transition/transversion multiplier and the exponential branch-prior rate. Branch lengths are in substitutions per site. They need not be ultrametric, and no dates or biological root are required.
The topology prior is uniform over labelled unrooted binary trees, often called PDA. Every actual edge independently has an exponential length conditional on the sampled branch-prior rate. The stored binary root is placed at the first tip's endpoint, with one computational zero edge; that edge adds neither a parameter nor a prior density. Newick output is marked `[&U]`.
Use `--model jc69`, `f81`, `hky` or `gtr` to select the DNA process. F81 infers frequencies, HKY also infers kappa, and GTR infers frequencies and six exchangeabilities constrained to sum to one. All have stationary mean substitution rate one. Unrooted inference requires a stationary reversible process; the C API rejects incompatible supplied models.
To supply fixed values instead of sampling parameters:
```sh
build/release/phyloSMC alignment.fasta --model hky --parameters fixed \
  --frequencies 0.3,0.2,0.2,0.3 --kappa 3 --branch-rate 10
```
With branch-prior rate ten fixed, the mean edge length is 0.1. Under joint inference, the default rate prior is lognormal with median ten and log standard deviation one; its unconditional mean edge length is approximately 0.165. The distinction matters when interpreting or changing the prior.
## Rooted analyses
Choose `--root-mode movable` to sample rooted topologies and allow the root split to change. This uses uniform rooted PDA with independent exponential lengths on its actual branches. With a reversible model, the sequence likelihood alone does not locate the root; the rooted prior remains part of the analysis. Output is marked `[&R]`.
Choose `--root-mode fixed` to condition on two groups separated by the root. Relationships within both groups remain free. An explicit outgroup can contain one or several comma-separated names:
```sh
build/release/phyloSMC alignment.fasta --root-mode fixed --outgroup A,B
```
Alternatively, `--root-mode fixed --tree rooted.nwk` extracts the root's two groups from a supplied Newick tree. That tree must have a binary root and exactly the observation taxa. Only its root split is retained; its internal topology and branch lengths are not fixed.
If neither source is supplied, fixed-root mode constructs a neighbour-joining tree from pairwise observed mismatch proportions, clamps negative limbs through Curupixa's NJ interface, and applies midpoint rooting. Ambiguous observations are excluded separately for each pair. A pair with no jointly observed unambiguous sites is an error. The selected split is printed to standard error so that it can be recorded and checked. Selecting a fixed split from the data makes the analysis conditional on that selection.
## Clocks and dated-tree priors
```sh
build/release/phyloSMC alignment.fasta --clock lognormal --particles 512
```
Clock analyses infer a rooted tree in time and convert each time duration to a substitution length using its branch rate. Available clocks are `strict`, `lognormal`, `exponential` and `autocorrelated`. The last model is a mean-preserving lognormal process between branch midpoints, defined explicitly in the tutorial. All current tip ages are zero. Time trees are ultrametric; relaxed-clock trees usually have unequal root-to-tip substitution distances.
The default dated-tree prior is a constant-size coalescent. Its pairwise waiting-time scale is inferred jointly. `--tree-prior yule` or `--tree-prior birth-death` instead samples conditioned reconstructed trees with complete sampling, a fixed origin age supplied through `--origin-age` (default ten), and inferred birth/death rates. Yule fixes death at zero. Fixed root splits also work with these clock analyses, using the correctly normalised conditional prior.
Mean clock rate is fixed at one by default, establishing relative time units. To infer it under a supplied proper prior, use `--infer-clock-mean yes`. Contemporary sequences alone do not separately identify absolute time and mean rate; interpretation then depends on the supplied scale assumptions. The inferred lognormal or autocorrelated dispersion has default prior median 0.3. Exponential clocks have no dispersion parameter.
## Prior settings
| Inferred parameter | Default prior |
| --- | --- |
| Branch-prior rate | Lognormal, median 10 and log SD 1 |
| Mean clock rate, when explicitly inferred | Lognormal, median 1 and log SD 1 |
| Clock dispersion | Lognormal, median 0.3 and log SD 1 |
| Coalescent pairwise scale | Lognormal, median 1 and log SD 1 |
| Birth rate / death rate | Lognormal medians 1 / 0.5, log SD 1 for each |
| HKY kappa | Lognormal, median 2 and log SD 1 |
| DNA frequencies / GTR exchangeabilities | Symmetric Dirichlet, concentration 1 for each coordinate |
Only parameters applicable to the chosen model are inferred. `--prior-medians` and `--prior-log-sds` take six comma-separated values in the table's scalar order: branch rate, clock mean, dispersion, coalescent scale, birth rate, death rate. Frequency and exchangeability concentrations have separate options. Kappa has separate prior-median and prior-log-SD options. `--help` lists all settings. Values such as `--branch-rate` supply fixed parameters when `--parameters fixed` is selected; they do not change hyperprior medians in joint mode.
## Processes on one machine
```sh
build/release/phyloSMC alignment.fasta --particles 1024 --moves 32 --workers 4
```
`--workers` explicitly enables POSIX worker processes; the default is one, running serially. Workers distribute prior sampling, likelihood evaluations and rejuvenation across particles. Temperature selection and resampling remain coordinated. The random streams are assigned to particle indices before dispatch, so a fixed seed reproduces the same population when worker count changes on the same build and platform.
Random-stream allocation differs from the initial release, so an old seed does not reproduce that earlier version's run.
Workers are created at stage boundaries. Process creation and tree transport can outweigh calculation time on small datasets. Full likelihood evaluations and full tree copies are still used; no large-data speedup is promised. In an embedded application, start multiprocessing analyses from a single-threaded process because this backend uses POSIX `fork`. No worker threads or MPI communication are started. On non-POSIX systems use one worker.
## Online addition
```sh
build/release/phyloSMC alignment.fasta --initial-taxa 4 --particles 1024 --moves 32
```
The first four taxa, in input order, are analysed by annealing. Remaining taxa are added one at a time. Online inference supports non-clock unrooted, movable-root and fixed-root analyses, including joint substitution and branch-prior parameters. A fixed-root prefix must already contain at least one taxon from each side; the complete split is retained as taxa arrive.
An arrival uses one importance update followed by rejuvenation. A difficult placement can concentrate weights severely, so inspect repeated runs and consider larger populations. Clocked or serially dated online insertion is not implemented.
## Output and interpretation
Standard output is a tab-separated table containing each particle's weight, parameter coordinates and Newick tree with substitution lengths. Inapplicable substitution coordinates are `nan`; fixed model frequencies are reported when the model has four states. A fixed borrowed model's kappa or exchangeabilities are not recovered through this table. Scalar coordinates belonging to an inactive clock/tree-prior family are retained settings and must not be interpreted as inferred values. Keep particle weights when calculating posterior means or clade probabilities.
Standard error reports log evidence, final effective sample size, stages, move acceptance and the number of surviving initial ancestors. Equal weights after resampling can accompany many copied particles. Surviving ancestors measure genealogy, whereas distinct tree configurations depend on subsequent moves. Compare independent runs and check whether topology and parameter summaries stabilise as particle count and rejuvenation increase.
No consensus or best tree is selected automatically. The C API exposes clock ages and branch rates in addition to substitution-length output trees.
## C API
The public option and particle structures have grown; rebuild downstream clients when upgrading.
Include `phylosmc.h`, initialise `psmc_options` with `psmc_options_init()`, and call `psmc_anneal()` with observations, a substitution model and an explicit uniform-random callback. Unlike the DNA CLI, the API starts with a fixed supplied substitution model and fixed hyperparameters. Set `substitution` to F81, HKY or GTR and set the applicable `infer` flags to enable joint inference; `priors` contains the corresponding prior settings. Default root mode selects unrooted non-clock or movable-root clock inference.
For a fixed split, set `root_mode=PSMC_ROOT_FIXED` and supply one `root_side` boolean per taxon. Results copy those assignments. `psmc_root_split()` extracts assignments from a supplied rooted tree; `psmc_estimate_root_split()` performs the NJ/midpoint procedure. `psmc_online_split()` accepts an extended fixed split; ordinary `psmc_online()` handles unrooted and movable-root additions.
Read the returned particles, weights and parameters, then free them with `psmc_result_free()`. Results retain observations but borrow the supplied immutable substitution model, which must remain alive until all dependent results are freed. A returned tree's node IDs refer to that particle; process transport and subsequent inference calls may rebuild representations. See the [API example](examples/api.c).
Embed with `add_subdirectory()` and link `phyloSMC::phylosmc`, or use an installed package through `find_package(phyloSMC CONFIG REQUIRED)`. All mutable numerical workspaces belong to the analysis; independent serial analyses can run concurrently with separate random streams.
Copyright (C) 2026 Leonardo de Oliveira Martins. Licensed under [MIT](LICENSE). Curupixa and its incorporated permissive components retain their respective notices; see [the licence review](docs/licence-review.md).
