# Bayesian phylogenetic inference with sequential Monte Carlo
This chapter defines the probability distributions used by phyloSMC and derives the calculations needed to sample from them. The annealed sampler follows the SMC-sampler construction of Del Moral, Doucet and Jasra [1] and its phylogenetic application by Wang, Wang and Bouchard-Côté [2]. Online taxon addition follows the importance-sampling principles developed for phylogenetics by Fourment and colleagues [3]. The particular priors and proposals below specify this implementation; adopting the general SMC construction does not require adopting those choices.
## 1. State space and probability measures
Let $D$ be an alignment of $n$ named taxa. Write $T$ for a labelled, non-plane topology, $\mathbf b$ for substitution lengths, $\theta$ for substitution parameters and $\eta$ for parameters of the prior. For a non-clock analysis, the state is
$$
x=(T,\mathbf b,\theta,\eta).
$$
The topology has a discrete probability mass. Lengths have a density with respect to Lebesgue measure on positive edge lengths. Frequency vectors and GTR exchangeability vectors have densities on their respective simplices. Specifying these measures matters when counting topologies, conditioning on a split or changing coordinates in a proposal.
The posterior is
$$
\pi(x\mid D)=\frac{L(D\mid x)p(x)}{Z},\qquad
Z=\sum_T\int L(D\mid x)p(x)\,d\mathbf b\,d\theta\,d\eta.
$$
Here $Z$ is the marginal likelihood, or evidence. Every prior factor must be normalised if $Z$ is to have this interpretation. Constants that cancel in an MH acceptance ratio can remain essential for evidence calculations.
### 1.1 Unrooted and rooted topologies
For $n\geq3$, an unrooted binary tree has $2n-3$ edges and $(2n-5)!!$ possible labelled topologies. Its uniform, or PDA, topology probability is
$$
p_U(T)=\frac{1}{(2n-5)!!}.
$$
The two-tip tree has one edge and one topology. By comparison, a rooted binary tree has $2n-2$ non-stem edges and
$$
p_R(T)=\frac{1}{(2n-3)!!}.
$$
PDA abbreviates “proportional to distinguishable arrangements”. Uniformity refers to labelled topologies, rather than unlabelled shapes.
An orientation stored in a C tree object does not itself add a biological root to the statistical state. For unrooted inference, phyloSMC places its binary storage root at the endpoint of the first taxon's pendant edge. One incoming edge then has length zero and is solely computational. The remaining $2n-3$ edges are the actual unrooted edges. Sampling a uniform rooted topology on the other $n-1$ tips gives exactly $(2n-5)!!$ alternatives, each representing one unrooted topology.
This convention avoids an additional length parameter. Assigning independent exponential priors to both parts of an artificially split edge would instead give their sum a gamma distribution, changing the intended prior.
### 1.2 Conditioning on a root split
Let $C=A\mid B$ specify a root bipartition, with $|A|=a$ and $|B|=b=n-a$. Define $R(m)=(2m-3)!!$, with $R(1)=1$. Under rooted PDA,
$$
P_R(C)=\frac{R(a)R(b)}{R(n)},\qquad
p_R(T\mid C)=\frac{\mathbf1\{T\text{ has root split }C\}}{R(a)R(b)}.
$$
Independent uniform rooted resolutions within the two sides therefore sample the conditional topology directly. Conditioning on the singleton split containing the first tip gives the topology distribution used for unrooted inference; assigning the computational zero edge then gives its length parameterisation.
A fixed biological split conditions the analysis on $C$. A movable-root analysis samples rooted trees without this restriction. The initial estimate used to select a fixed split can come from midpoint-rooted neighbour joining, an explicit outgroup, or a supplied rooted tree. When the split is selected from the alignment, the resulting analysis is conditional on that data-dependent selection; it does not integrate uncertainty about the selection.
### 1.3 Selecting an automatic fixed division
For each taxon pair, let $S_{ij}$ contain the jointly observed, unambiguous site patterns and define the mismatch proportion
$$
D_{ij}=\frac{\sum_{s\in S_{ij}}m_s\mathbf1\{y_{is}\ne y_{js}\}}{\sum_{s\in S_{ij}}m_s}.
$$
A pair with an empty denominator is an error. Curupixa constructs a neighbour-joining tree from these distances [11], using its documented non-negative limb policy. Midpoint rooting then selects the halfway point of a diameter path, and the two groups on that edge define the constraint. The mismatch proportions provide a rough initial division; no correction for multiple substitutions is applied. This selection is only invoked when a fixed split is requested without an explicit source.
## 2. Sequence evolution and the pruning likelihood
For a finite state alphabet, a continuous-time Markov process has generator $Q$, with $q_{ij}\geq0$ for $i\ne j$ and
$$
q_{ii}=-\sum_{j\ne i}q_{ij},\qquad P(b)=\exp(Qb).
$$
The entry $P_{ij}(b)$ is the probability of child state $j$ given parent state $i$ across length $b$. For normalised stationary models,
$$
\boldsymbol\pi Q=0,\qquad -\sum_i\pi_iq_{ii}=1,
$$
so $b$ is the expected number of substitutions per site under stationarity.
### 2.1 JC69, F81, HKY and GTR
JC69 has equal DNA frequencies and equal exchangeabilities. F81 allows positive frequencies $\pi_A,\pi_C,\pi_G,\pi_T$ summing to one, with common exchangeability. HKY additionally multiplies the A–G and C–T exchangeabilities by $\kappa>0$ [4]. GTR uses six positive symmetric exchangeabilities $s_{ij}=s_{ji}$:
$$
\widetilde q_{ij}=s_{ij}\pi_j,\qquad
Q=\frac{\widetilde Q}{-\sum_i\pi_i\widetilde q_{ii}}.
$$
An overall multiplier on all $s_{ij}$ cancels. phyloSMC removes this redundant coordinate by constraining the six exchangeabilities to sum to one. Curupixa's GTR constructor performs the generator normalisation; its generic `rates4` constructor preserves the supplied rate scale instead.
For these stationary reversible processes,
$$
\pi_iP_{ij}(b)=\pi_jP_{ji}(b).
$$
Consequently, the sequence likelihood is independent of computational rooting. Unrooted inference rejects a supplied model without stationary reversibility. In a rooted analysis, a non-reversible model or a non-stationary root distribution can make root placement relevant to the likelihood.
### 2.2 Conditional likelihood recursion
For pattern $s$, let $\ell_{v,s}(i)$ be the likelihood of the observations descending from vertex $v$, conditional on state $i$ at $v$. At a tip, the observation vector supplies these values: an observed nucleotide has a single entry one; an ambiguity code has one for each allowed state; a fully missing observation has all entries one.
For an internal vertex,
$$
\ell_{v,s}(i)=\prod_{c\in\operatorname{children}(v)}
\left[\sum_j P_{ij}(b_c)\ell_{c,s}(j)\right].
$$
At a computational root $r$ with root distribution $\rho$,
$$
L_s=\sum_i\rho_i\ell_{r,s}(i),\qquad
\log L(D\mid x)=\sum_s m_s\log L_s,
$$
where $m_s$ is the multiplicity of pattern $s$. This is Felsenstein's pruning calculation [5]. Curupixa evaluates it in logarithmic coordinates to avoid underflow. Combining identical patterns changes the amount of work without changing the likelihood.
The implementation uses one substitution process across the alignment. Across-site gamma mixtures, invariant-site mixtures, codon processes and empirical protein matrices require additional model interfaces; they are not implied by having a generic SMC engine.
## 3. Non-clock branch and parameter priors
Conditional on the exponential rate $\lambda$, actual non-clock edges independently satisfy
$$
p(\mathbf b\mid T,\lambda)=\prod_{e\in E(T)}\lambda\exp(-\lambda b_e),\qquad b_e>0.
$$
The computational zero edge of an unrooted representation is excluded. Rooted analyses include both biological root edges and exclude the root stem.
When $\lambda$ is inferred, its prior is lognormal. For any inferred positive scalar $z$ with prior median $m$ and log standard deviation $s$,
$$
\log z\sim N(\log m,s^2),\qquad
\log p(z)=-\log z-\log s-\tfrac12\log(2\pi)
-\frac{(\log z-\log m)^2}{2s^2}.
$$
The $-\log z$ term belongs to the density in the original coordinate $z$ and must be included. The default prior for $\lambda$ has $m=10,s=1$. Thus the conditional mean branch length is $1/\lambda$, while the unconditional mean is
$$
E[b_e]=E[\lambda^{-1}]=\frac{\exp(s^2/2)}{m}.
$$
It is approximately $0.165$ under the joint default, rather than the $0.1$ obtained by fixing $\lambda=10$.
DNA frequencies have a symmetric Dirichlet prior:
$$
p(\boldsymbol\pi)=\frac{\Gamma(4\alpha)}{\Gamma(\alpha)^4}
\prod_{i=1}^4\pi_i^{\alpha-1}.
$$
GTR exchangeabilities have the analogous six-dimensional prior with concentration $\alpha_s$. Both concentrations default to one. HKY uses $\log\kappa\sim N(\log2,1)$. Curupixa's inverse gamma and normal CDFs supply the initial parameter draws; phyloSMC supplies the biological prior choices.
## 4. Clock states, tree times and rates
For a clock analysis, the state contains node ages $\mathbf t$ and branch rates $\mathbf r$:
$$
x=(T,\mathbf t,\mathbf r,\theta,\eta),\qquad
b_e=r_e\Delta t_e,\qquad \Delta t_e=t_{\mathrm{parent}(e)}-t_e>0.
$$
The target density is specified directly in ages and rates:
$$
p(x)=p(T,\mathbf t\mid\eta_T)\,
p(\mathbf r\mid T,\mathbf t,\eta_r)\,
p(\theta,\eta_T,\eta_r).
$$
There is no change-of-variables Jacobian for evaluating the likelihood at $b_e=r_e\Delta t_e$: the likelihood is a function evaluated at a deterministic transformation. A Jacobian would be needed if we instead transformed a density originally specified in substitution lengths into age/rate coordinates.
All tips currently have age zero. The time tree is therefore ultrametric. A strict clock also gives equal root-to-tip substitution distances. Relaxed rates generally give unequal substitution distances even though the time tree remains ultrametric.
### 4.1 Constant-size coalescent
With $k$ ancestral lineages and pairwise time scale $s_c$,
$$
\Delta t_k\sim\operatorname{Exp}\left(\frac{\binom{k}{2}}{s_c}\right).
$$
Choosing a merging pair uniformly gives the joint density of one specified pair and its waiting time:
$$
\frac1{s_c}\exp\left[-\frac{\binom{k}{2}}{s_c}\Delta t_k\right].
$$
For a labelled ranked history, multiply these factors from $k=n$ to $2$. The total root age has expectation
$$
E[t_{\mathrm{root}}]=s_c\sum_{k=2}^n\binom{k}{2}^{-1}
=2s_c\left(1-\frac1n\right).
$$
The scale corresponds, for example, to $2N_e$ generations under an appropriate diploid convention. Estimating this scale in relative clock units does not by itself estimate an absolute census population size. Kingman gives the underlying coalescent construction [6].
### 4.2 Conditioned Yule and birth–death trees
Curupixa's reconstructed birth–death prior conditions on $n$ contemporary sampled tips, complete sampling and a fixed origin age $O$. Extinct lineages and the origin stem are excluded. With birth rate $\lambda_b$, death rate $\mu_d$ and $d=\lambda_b-\mu_d\ne0$, define
$$
F(t)=\frac{\lambda_b(1-e^{-dt})}{\lambda_b-\mu_d e^{-dt}},\qquad
f(t)=\frac{\lambda_b d^2e^{-dt}}{(\lambda_b-\mu_d e^{-dt})^2}.
$$
At criticality, $\lambda_b=\mu_d$,
$$
F(t)=\frac{\lambda_b t}{1+\lambda_b t},\qquad
f(t)=\frac{\lambda_b}{(1+\lambda_b t)^2}.
$$
The ordered branching ages $0<t_1<\cdots<t_{n-1}<O$ have density $(n-1)!\prod_i f(t_i)/F(O)$. Multiplying by the uniform ranked-pair probabilities gives the density of the corresponding labelled tree. Yule is obtained by setting $\mu_d=0$. These conditioning choices follow the reconstructed-process treatment discussed by Gernhard [7]; other conditioning choices give different priors and evidence values.
Birth and death rates can be inferred jointly with the tree, using proper lognormal priors. The default medians are one and one-half, respectively, with log standard deviations one. The origin age remains supplied. Critical and subcritical parameter combinations are allowed by the implementation.
### 4.3 Dated trees conditional on a root split
For the ranked topology distribution shared by these dated priors, the probability of a particular root split $A\mid B$ is
$$
P(C)=\frac{2}{(n-1)\binom n a}.
$$
The conditional density is $p(T,\mathbf t)\mathbf1_C/P(C)$. Branching times are independent of the ranked labelled topology, so they retain their marginal distribution after conditioning.
A direct conditional sampler maintains $a'$ and $b'$ active lineages on the two sides. Until only one lineage remains on each side, choose the next within-side merger on side A with probability
$$
\frac{a'-1}{a'+b'-2},
$$
and otherwise on side B; then choose the pair uniformly within that side. Finally merge the two remaining lineages at the root. These are the ordinary pair-merger probabilities multiplied by the ratio of probabilities of eventual satisfaction of the split. This construction avoids rejection sampling for rare splits.
### 4.4 Branch-rate models
A strict clock sets $r_e=\bar r$ on every branch. Its rate distribution is a point mass, and only $\bar r$ can be inferred as a rate parameter.
For independent lognormal rates with arithmetic mean $\bar r$ and log standard deviation $\sigma$,
$$
\log r_e\sim N\left(\log\bar r-\frac{\sigma^2}{2},\sigma^2\right).
$$
Independent exponential rates satisfy $r_e\sim\operatorname{Exp}(1/\bar r)$. These families correspond to uncorrelated relaxed-clock choices discussed by Drummond and colleagues [8].
For the implemented autocorrelated model,
$$
\log r_e\mid r_p\sim N\left(\log r_p-\frac{v_e}{2},v_e\right),\qquad
v_e=\sigma^2\frac{\Delta t_e+\Delta t_p}{2}.
$$
At a root child, replace $r_p$ by $\bar r$ and use $v_e=\sigma^2\Delta t_e/2$. The root's rate entry is fixed at $\bar r$ and contributes no density factor. This branch-midpoint definition gives $E[r_e\mid r_p]=r_p$. It uses a constant rate along each edge; the substitution length is not an integral of a continuously varying geometric Brownian trajectory. Thorne, Kishino and Painter provide the historical autocorrelated-clock setting [9]; the exact mean-preserving convention above specifies our model.
The default inferred dispersion prior has median $0.3$ and log standard deviation one. Coalescent scale has median one and log standard deviation one. The mean clock rate remains fixed at one unless its inference is explicitly requested, in which case its lognormal prior also has median one and log standard deviation one.
### 4.5 Time-scale identifiability
For any $c>0$, the transformation $\Delta t_e\mapsto c\Delta t_e$, $r_e\mapsto r_e/c$ preserves every substitution length. Contemporary sequence observations alone therefore cannot generally identify separate absolute rate and time scales. Fixing the mean rate establishes relative time units. Inferring both scales with proper priors produces a proper statistical specification, but their absolute interpretation then depends on those priors or supplied temporal information.
## 5. Annealed SMC
Choose temperatures $0=\beta_0<\beta_1<\cdots<\beta_K=1$ and define
$$
\gamma_\beta(x)=p(x)L(D\mid x)^\beta,\qquad
Z_\beta=\int\gamma_\beta(x)\,dx,\qquad
\pi_\beta=\gamma_\beta/Z_\beta.
$$
With a normalised prior, $Z_0=1$ and $Z_1=Z$. Initialise $N$ particles independently from the joint prior and set $W_0^{(i)}=1/N$. The prior includes all inferred hyperparameters before any tree or branch is drawn conditionally on them.
### 5.1 Incremental importance weights
Suppose the current population approximates $\pi_\beta$. Before moving the particles, changing temperature by $\Delta$ gives incremental weight
$$
u_i=\frac{\gamma_{\beta+\Delta}(x_i)}{\gamma_\beta(x_i)}
=L(D\mid x_i)^\Delta.
$$
The normalised weights and evidence increment are
$$
c=\sum_i W_i u_i,\qquad W'_i=\frac{W_i u_i}{c},\qquad
\widehat Z\leftarrow\widehat Zc.
$$
In logarithmic coordinates,
$$
a_i=\log W_i+\Delta\log L_i,\qquad
\log c=m+\log\sum_i e^{a_i-m},\quad m=\max_i a_i.
$$
Adding $\log c$ to the running log evidence avoids direct multiplication of very small quantities.
### 5.2 ESS and conditional ESS
For normalised weights,
$$
\operatorname{ESS}(W)=\frac1{\sum_i W_i^2}.
$$
This measures weight concentration in the current population. For choosing the next temperature, conditional ESS measures the additional concentration caused by a candidate increment:
$$
\operatorname{CESS}(\Delta)=
N\frac{\left(\sum_i W_i L_i^\Delta\right)^2}
{\sum_i W_i L_i^{2\Delta}}.
$$
At $\Delta=0$, CESS equals $N$ even when the existing weights are unequal. ESS and CESS coincide after an increment only when the previous weights were equal. phyloSMC uses bisection to choose the largest allowed temperature satisfying $\operatorname{CESS}\geq\rho N$, with default $\rho=0.8$ [10].
### 5.3 Resampling
If ESS falls below $\tau N$, with default $\tau=0.5$, systematic resampling draws one offset $U\sim\operatorname{Uniform}(0,1/N)$ and places points
$$
U,\ U+1/N,\ldots,U+(N-1)/N
$$
in the cumulative weight intervals. Each selected particle is copied and receives weight $1/N$. Expected offspring counts are $NW_i$, preserving weighted expectations in the resampling step. This does not create new tree configurations; it duplicates existing particles.
The reported surviving-ancestor count tracks distinct members of the initial population represented by current descendants. Rejuvenation can create different trees among descendants of the same ancestor, so this count diagnoses genealogical loss rather than counting distinct current topologies.
### 5.4 Evidence and adaptation
For a fixed temperature schedule and suitable unbiased resampling, standard SMC normalising-constant estimators have established unbiasedness properties [1]. Here temperatures are selected from the evolving population. This implementation does not assert exact finite-population unbiasedness under that adaptation. Moreover, even when $\widehat Z$ is unbiased, Jensen's inequality gives $E[\log\widehat Z]\leq\log Z$. Independent runs help assess numerical variability; agreement does not prove adequate exploration of every posterior mode.
### 5.5 Order of operations
```text
Draw each particle's parameters from their hyperpriors.
Draw its topology and lengths, or times and rates, conditionally on those parameters.
Set beta = 0, all weights = 1/N, and log_evidence = 0.
While beta < 1:
  Select next_beta by the CESS criterion using the current likelihoods and weights.
  Form log_weights[i] = log(weights[i]) + (next_beta-beta)*log_likelihood[i].
  Normalise them and add their log normalising constant to log_evidence.
  Set beta = next_beta.
  If ESS is below the resampling threshold, resample and reset weights to 1/N.
  Apply the requested number of MH moves to each positive-weight particle at this beta.
Return particles, their retained weights, and log_evidence.
```
Recording the normalising constant before resampling preserves the evidence increment. Applying moves at the new temperature preserves the distribution to which the particles were just reweighted.
## 6. Metropolis–Hastings rejuvenation
After reweighting and optional resampling, use a transition kernel invariant for the current $\pi_\beta$. For proposal $q(x'\mid x)$,
$$
\alpha(x,x')=\min\left(1,
\exp\left[\beta(\log L'-\log L)+\log p'-\log p
+\log\frac{q(x\mid x')}{q(x'\mid x)}\right]\right).
$$
The likelihood term is tempered; prior and proposal terms are not. Invariance follows from detailed balance of the accepted transition probabilities. Finite rejuvenation does not imply that each particle has independently equilibrated at the new temperature.
### 6.1 Prior independence moves
Draw the entire proposed state from its joint prior. Since $q(x'\mid x)=p(x')$,
$$
\alpha=\min\left(1,(L'/L)^\beta\right).
$$
For a fixed split, use the conditional prior in both the target and proposal. All inferred parameters are redrawn along with the tree. Such moves establish access beyond a local neighbourhood, although their acceptance can be low for informative alignments.
### 6.2 Topology moves
NNI exchanges one child subtree across an eligible internal edge. Selecting an eligible edge and either alternative uniformly gives a symmetric proposal because the number of eligible edges is constant for each supported binary topology class. In unrooted mode, exclude moves across the computational root edge. In fixed-root mode, exclude exchanges across either root child so that the two groups remain intact. Other NNI moves can rearrange relationships within the groups. Movable-root proposals permit changes to the root split.
Clock proposals retain node-associated age and rate coordinates. Every topology change is followed by validation of temporal ordering and recalculation of the full tree and clock prior, including changed autocorrelations.
### 6.3 Positive scalar multipliers
For a branch length, relaxed branch rate or positive hyperparameter, propose
$$
z'=ze^\epsilon,\qquad \epsilon\sim\operatorname{Uniform}(-w,w).
$$
The proposal density with respect to $z'$ is $1/(2wz')$ on its support. Hence
$$
\frac{q(z\mid z')}{q(z'\mid z)}=\frac{z'}z=e^\epsilon.
$$
The log acceptance ratio must include $\epsilon$. For a hyperparameter, recompute every density that depends on it: changing $\lambda$ changes all branch priors, changing a clock dispersion changes rate densities, and changing coalescent scale changes the dated-tree density.
### 6.4 Simplex and age moves
For frequency or exchangeability coordinates, choose distinct positions $i,j$ uniformly. Their sum $s=z_i+z_j$ remains fixed. Propose $z'_i=z_i+\epsilon s$ and $z'_j=s-z'_i$, rejecting non-positive values. The reverse proposal has the same scale because $s$ is unchanged, so its Hastings factor is one.
An internal age uses an additive symmetric proposal $t'_v=t_v+\epsilon$. Reject it if it is not older than every child and younger than its parent, or if it violates the dated prior's support. The root has no parent bound but remains subject to its prior. Strict-clock mean-rate proposals update all deterministic branch-rate entries to the new mean; they do not treat these copies as independent random coordinates.
## 7. Online taxon addition
Let $x_n$ describe a posterior state for $n$ taxa. Retain its hyperparameters and substitute a proposed expanded tree $x_{n+1}=F(x_n,u)$, where $u$ contains the insertion choice and new continuous coordinates. The incremental importance factor is
$$
G=\frac{L_{n+1}(x_{n+1})p_{n+1}(x_{n+1})}
{L_n(x_n)p_n(x_n)q(u\mid x_n)}\,|J_F|.
$$
Equivalently, express $q$ directly as a density in the new edge coordinates and absorb the Jacobian there. phyloSMC uses this second convention.
### 7.1 Inserting on an existing edge
For an old edge of length $L$, choose $f\sim\operatorname{Uniform}(0,1)$ and split it into $a=fL$, $b=(1-f)L$. Draw pendant length $y\sim\operatorname{Exp}(\lambda)$. The transformation $(L,f,y)\mapsto(a,b,y)$ has absolute Jacobian $L$. Thus the proposal contribution in new coordinates is
$$
q_e(a,b,y\mid L)=\frac1L\lambda e^{-\lambda y}.
$$
The branch-prior ratio is
$$
\frac{\lambda^3 e^{-\lambda(a+b+y)}}{\lambda e^{-\lambda L}}
=\lambda^2e^{-\lambda y}.
$$
Dividing by the proposal leaves $\lambda L$. Omitting this factor changes the target distribution. When $\lambda$ is inferred, use that particle's retained $\lambda$.
### 7.2 Topology factors and root modes
For unrestricted rooted PDA, choose uniformly among $2n-1$ positions: $2n-2$ edges and the position above the root. The topology probability ratio is $1/(2n-1)$, cancelling that position-choice probability. Above the root, both new lengths are drawn independently from their exponential priors, leaving branch factor one.
For unrooted PDA, choose among its $2n-3$ actual edges. The topology ratio is $1/(2n-3)$ and again cancels the choice probability. There is no above-root insertion and no insertion on the computational zero edge.
For fixed-root inference, the extended side assignments must restrict to the previous assignments. The code retains the general rooted insertion proposal and assigns zero target weight to insertions violating the extended split. It evaluates the ratio of conditional topology probabilities explicitly; the conditioning constants change with taxon count and group sizes and cannot be omitted. Initial prefixes must already contain a member of each side.
### 7.3 Weight and evidence update
For each expanded particle,
$$
\log\widetilde W_i=\log W_i+\log L_{n+1,i}-\log L_{n,i}
+\log p_{n+1,i}-\log p_{n,i}-\log q_i.
$$
Normalising these weights estimates $Z_{n+1}/Z_n$. Resampling and rejuvenation then target the full expanded posterior. Existing relationships can change during rejuvenation.
The independent-edge tree priors for successive taxon counts are not generally projectively consistent: removing a tip can join two exponential edges into a gamma-distributed edge. Consequently, adding a completely unobserved tip need not leave the evidence unchanged for informative existing observations. The full prior ratio accounts for this difference.
Online addition currently supports non-clock trees, including joint substitution and branch-rate inference. Serially dated or clocked insertion needs age/rate proposals and their corresponding density calculations; it remains outside this implementation.
## 8. An analytic two-tip example
For a symmetric binary substitution process normalised to exit rate one,
$$
P_{01}(d)=\frac{1-e^{-2d}}2.
$$
At stationarity, observing tip states $(0,1)$ has likelihood $L(d)=\tfrac14(1-e^{-2d})$. Under an unrooted single-edge prior $d\sim\operatorname{Exp}(\lambda)$,
$$
Z_U=\frac14\left(1-\frac{\lambda}{\lambda+2}\right).
$$
Under a rooted two-edge prior, the total distance is $d=b_1+b_2$, so
$$
Z_R=\frac14\left[1-\left(\frac{\lambda}{\lambda+2}\right)^2\right].
$$
The likelihood is root invariant in both cases. The evidence differs because the priors on the total distance differ. These identities provide independent tests of the implemented probability measures.
As an exercise, derive the rooted three-tip evidence after adding a missing third tip. Under rooted PDA, the observed pair is a cherry with probability $1/3$, giving two intervening edges; otherwise its path has three edges. With $r=\lambda/(\lambda+2)$,
$$
Z_{R,3}=\frac14\left(1-\frac{r^2+2r^3}{3}\right).
$$
The corresponding unrooted three-tip tree has a two-edge observed path, giving $Z_{U,3}=\tfrac14(1-r^2)$. Both online calculations are checked against these expressions.
## 9. Multiprocessing and reproducibility
Within each stage, prior draws and particle rejuvenation can be assigned to different processes. Temperature selection, weight normalisation and resampling depend on the complete population and remain coordinated. Weight reductions use a fixed order in the parent process.
Before dispatch, the parent derives one seed from its random callback and allocates jumped xoshiro256++ streams to particle indices. Stream allocation is independent of the number of workers. Each particle receives the same stream and proposal order for a given seed and population, so changing the worker count does not change its result on the same platform and build. This does not promise identical floating-point arithmetic across different compilers or architectures.
POSIX workers receive private mutable workspaces through `fork()`. They return trees and parameter state through pipes. Trees are serialised as ordered Newick, with age/rate arrays attached by preorder; serial execution performs the same representation normalisation to preserve node ordering. Process creation and transport have a cost, and short analyses can be slower with several workers. There is no MPI or network transport. Because workers use POSIX `fork`, start a multiprocessing analysis from a single-threaded parent process; serial library calls remain available in other hosts.
Full tree copies and full likelihood calculations are retained in this version. Incremental likelihood caching would reduce work per proposal independently of process-level distribution.
## 10. From equations to source
| Calculation | Source |
| --- | --- |
| Tree priors, split conditioning and conditional sampling | Curupixa `tree_prior.c` |
| Substitution matrices and pruning | Curupixa `substitution_model.c`, `likelihood.c` |
| Clock probabilities and age/rate conversion | Curupixa `molecular_clock.c`, `clock.c` |
| CESS, normalisation and resampling indices | Curupixa `particle_weights.c` |
| Temperature stages, evidence and online weights | phyloSMC `src/phylosmc.c` |
| Particle priors and likelihood evaluation | phyloSMC `src/particle.c` |
| Hyperpriors and parameter proposals | phyloSMC `src/parameters.c` |
| Topology/age/branch moves and resampling | phyloSMC `src/moves.c` |
| Process dispatch and random streams | phyloSMC `src/workers.c` |
The tests establish analytic evidence identities, prior moments, topology probabilities and reproducibility contracts. Adequate exploration for a particular alignment must be assessed separately using independent runs, weighted clade/parameter summaries, and sensitivity to particle count and rejuvenation effort.
## References
[1] Del Moral P, Doucet A, Jasra A (2006). Sequential Monte Carlo samplers. *Journal of the Royal Statistical Society B* 68: 411–436. [doi:10.1111/j.1467-9868.2006.00553.x](https://doi.org/10.1111/j.1467-9868.2006.00553.x).

[2] Wang L, Wang S, Bouchard-Côté A (2020). An annealed sequential Monte Carlo method for Bayesian phylogenetics. *Systematic Biology* 69: 155–183. [doi:10.1093/sysbio/syz028](https://doi.org/10.1093/sysbio/syz028).

[3] Fourment M, Claywell BC, Dinh V, McCoy C, Matsen FA IV, Darling AE (2018). Effective online Bayesian phylogenetics via sequential Monte Carlo with guided proposals. *Systematic Biology* 67: 490–502. [doi:10.1093/sysbio/syx090](https://doi.org/10.1093/sysbio/syx090). Our uniform placement proposal implements the general importance construction, rather than that paper's guided-placement algorithms.

[4] Hasegawa M, Kishino H, Yano T (1985). Dating of the human–ape splitting by a molecular clock of mitochondrial DNA. *Journal of Molecular Evolution* 22: 160–174. [doi:10.1007/BF02101694](https://doi.org/10.1007/BF02101694).

[5] Felsenstein J (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. *Journal of Molecular Evolution* 17: 368–376. [doi:10.1007/BF01734359](https://doi.org/10.1007/BF01734359).

[6] Kingman JFC (1982). On the genealogy of large populations. *Journal of Applied Probability* 19: 27–43. [doi:10.2307/3213548](https://doi.org/10.2307/3213548).

[7] Gernhard T (2008). The conditioned reconstructed process. *Journal of Theoretical Biology* 253: 769–778. [doi:10.1016/j.jtbi.2008.04.005](https://doi.org/10.1016/j.jtbi.2008.04.005).

[8] Drummond AJ, Ho SYW, Phillips MJ, Rambaut A (2006). Relaxed phylogenetics and dating with confidence. *PLOS Biology* 4: e88. [doi:10.1371/journal.pbio.0040088](https://doi.org/10.1371/journal.pbio.0040088).

[9] Thorne JL, Kishino H, Painter IS (1998). Estimating the rate of evolution of the rate of molecular evolution. *Molecular Biology and Evolution* 15: 1647–1657. [doi:10.1093/oxfordjournals.molbev.a025892](https://doi.org/10.1093/oxfordjournals.molbev.a025892).

[10] Zhou Y, Johansen AM, Aston JAD (2016). Toward automatic model comparison: an adaptive sequential Monte Carlo approach. *Journal of Computational and Graphical Statistics* 25: 701–726. [doi:10.1080/10618600.2015.1060885](https://doi.org/10.1080/10618600.2015.1060885).

[11] Saitou N, Nei M (1987). The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Molecular Biology and Evolution* 4: 406–425. [doi:10.1093/oxfordjournals.molbev.a040454](https://doi.org/10.1093/oxfordjournals.molbev.a040454).
