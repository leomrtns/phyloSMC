/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
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
typedef struct {
  size_t particles, mcmc_steps, max_stages;
  double conditional_ess, resample_ess, proposal_width, independence_probability, branch_prior_rate;
  bool use_clock;
  crpx_clock_model clock;
  crpx_tree_prior tree_prior;
} psmc_options;
typedef struct {
  crpx_tree *tree;
  double *ages, *rates;
  double log_likelihood, log_prior;
} psmc_particle;
typedef struct {
  size_t count, stages, proposals, accepted;
  double beta, ess, log_evidence;
  psmc_particle *particles;
  double *weights;
  psmc_options options;
  const crpx_character_data *data;
  const crpx_substitution_model *model;
} psmc_result;
/* Default non-clock: rooted PDA topology, independent Exp(10) substitution lengths.
   Clock mode uses the supplied dated tree prior and absolute branch-rate model.
   Hyperparameters and substitution model remain fixed. Default clock prior is
   constant-size coalescent, with a strict rate of one. No worker threads start. */
void psmc_options_init (psmc_options *options);
/* Draw from the prior, adapt temperatures by conditional ESS, optionally resample,
   then rejuvenate with MH moves. Models may use any Curupixa state alphabet.
   Result retains data; model is borrowed and must remain alive and immutable until
   result destruction. Trees store substitution lengths in non-clock mode and time
   lengths in clock mode. Result arrays own their memory and may be read directly; treat retained results as read-only.
   Failed calls set *out=NULL. Random draws consumed on failure are not restored. */
crpx_status psmc_anneal (const crpx_character_data *data, const crpx_substitution_model *model,
                         const psmc_options *options, crpx_uniform_function uniform, void *state, psmc_result **out);
/* Add one taxon (last in data). Earlier names, sites, observation weights and model
   must be unchanged; unit site weights are required. new data may compact patterns
   differently. Non-clock mode only. Propose an edge uniformly including above-root,
   split the old length uniformly and draw new lengths from their exponential prior.
   Correct weights include the split-length density/Jacobian. Rejuvenate at beta=1.
   Input remains unchanged, including on failure. log_evidence adds the predictive
   evidence increment to the previous estimate. No new annealing bridge is claimed. */
crpx_status psmc_online (const psmc_result *previous, const crpx_character_data *data, crpx_uniform_function uniform,
                         void *state, psmc_result **out);
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
