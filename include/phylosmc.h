/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
/*! \file phylosmc.h
    \brief Annealed and online phylogenetic SMC using Curupixa likelihoods and priors.
    Wang et al. (2020), doi:10.1093/sysbio/syz028; Fourment et al. (2018), doi:10.1093/sysbio/syx090. */
#ifndef PHYLOSMC_H
#define PHYLOSMC_H
#include <inference.h>
#include <likelihood.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { PSMC_ROOT_DEFAULT, PSMC_UNROOTED, PSMC_ROOT_MOVABLE, PSMC_ROOT_FIXED } psmc_root_mode;
typedef enum { PSMC_MODEL_FIXED, PSMC_MODEL_F81, PSMC_MODEL_HKY, PSMC_MODEL_GTR } psmc_model_kind;
enum {
  PSMC_INFER_BRANCH_RATE = 1u, PSMC_INFER_CLOCK_MEAN = 2u, PSMC_INFER_CLOCK_DISPERSION = 4u,
  PSMC_INFER_COALESCENT_SCALE = 8u, PSMC_INFER_BIRTH_RATE = 16u, PSMC_INFER_DEATH_RATE = 32u
};
typedef struct {
  double branch_rate, clock_mean, clock_dispersion, coalescent_scale, birth_rate, death_rate;
  double frequencies[4], kappa, exchangeabilities[6];
} psmc_parameters;
typedef struct {
  /* Lognormal medians and log SDs for the six positive parameters, in flag order. */
  double median[6], log_sd[6];
  double frequency_concentration, exchangeability_concentration, kappa_median, kappa_log_sd;
} psmc_hyperpriors;
typedef struct {
  size_t particles, mcmc_steps, max_stages;
  double conditional_ess, resample_ess, proposal_width, independence_probability, branch_prior_rate;
  bool use_clock;
  crpx_clock_model clock;
  crpx_tree_prior tree_prior;
  psmc_root_mode root_mode;
  /* One assignment per observation taxon, both sides nonempty. Required for FIXED; copied into results.
     For clock runs a fixed split conditions the dated prior, without fixing any branching ages. */
  const bool *root_side;
  psmc_model_kind substitution;
  unsigned infer;
  psmc_hyperpriors priors;
  size_t workers;
} psmc_options;
typedef struct {
  crpx_tree *tree;
  double *ages, *rates;
  double log_likelihood, log_prior;
  psmc_parameters parameters;
} psmc_particle;
typedef struct {
  size_t count, stages, proposals, accepted;
  double beta, ess, log_evidence;
  psmc_particle *particles;
  double *weights;
  psmc_options options;
  const crpx_character_data *data;
  const crpx_substitution_model *model;
  /* Distinct initial particle IDs still represented after resampling, not a count of distinct current trees. */
  size_t surviving_ancestors;
  size_t *ancestor_ids;
} psmc_result;
/* Default non-clock: unrooted PDA topology, independent Exp(10) substitution lengths.
   Clock mode uses the supplied dated tree prior and absolute branch-rate model.
   Set substitution/infer to sample parameters jointly. Default clock prior is
   constant-size coalescent, with a strict rate of one. workers=1 starts no processes or threads. */
void psmc_options_init (psmc_options *options);
/* Draw from the prior, adapt temperatures by conditional ESS, optionally resample,
   then rejuvenate with MH moves. Models may use any Curupixa state alphabet.
   DEFAULT selects unrooted non-clock or movable-root clock inference. Explicit UNROOTED requires a stationary
   reversible model and use_clock=false. Its binary storage root is at the first taxon's endpoint: one zero edge
   is computational, and the other 2n-3 edges have independent priors. FIXED conditions on root_side.
   Workers>1 explicitly enables POSIX fork workers; call from a single-threaded process. Random substreams are
   allocated by particle before work is dispatched, so worker count does not change results.
   Result retains data; model is borrowed and must remain alive and immutable until
   result destruction. Trees store substitution lengths in non-clock mode and time
   lengths in clock mode. Result arrays own their memory and may be read directly; treat retained results as read-only.
   Failed calls set *out=NULL. Random draws consumed on failure are not restored. */
crpx_status psmc_anneal (const crpx_character_data *data, const crpx_substitution_model *model,
                         const psmc_options *options, crpx_uniform_function uniform, void *state, psmc_result **out);
/* Add one taxon (last in data). Earlier names, sites, observation weights and model
   must be unchanged; unit site weights are required. new data may compact patterns
   differently. Non-clock mode only. Unrooted insertion chooses an actual edge; rooted insertion also includes
   above-root. Split an old length uniformly and draw new lengths from the particle's exponential prior.
   Correct weights include the split-length density/Jacobian. Rejuvenate at beta=1.
   Fixed-root online calls require psmc_online_split() with the extended split; its restriction must equal the old split.
   Input remains unchanged, including on failure. log_evidence adds the predictive
   evidence increment to the previous estimate. No new annealing bridge is claimed. */
crpx_status psmc_online (const psmc_result *previous, const crpx_character_data *data, crpx_uniform_function uniform,
                         void *state, psmc_result **out);
crpx_status psmc_online_split (const psmc_result *previous, const crpx_character_data *data, const bool *root_side,
                               crpx_uniform_function uniform, void *state, psmc_result **out);
/* Extract the two groups at the biological root of a supplied tree. Names must match observations exactly. */
crpx_status psmc_root_split (const crpx_tree *tree, const crpx_character_data *data, bool *side);
/* Construct NJ from observed mismatch proportions (ambiguous observations excluded per pair), root at midpoint,
   then extract the split. No jointly observed sites for a pair is an error. */
crpx_status psmc_estimate_root_split (const crpx_character_data *data, bool *side);
void psmc_result_free (psmc_result *result);
/* Return a cloned tree with substitution lengths, including in clock analyses. */
crpx_status psmc_substitution_tree (const psmc_result *result, size_t particle, crpx_tree **out);
/* Select the first taxa observations, preserving sites and ambiguity. Unit site
   weights required. Useful for sequential addition from one complete alignment. */
crpx_status psmc_data_prefix (const crpx_character_data *data, size_t taxa, crpx_character_data **out);
#ifdef __cplusplus
}
#endif
#endif
