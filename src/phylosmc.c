/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
#include "smc_internal.h"
#include <particle_weights.h>
#include <topology_randomise.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void
psmc_options_init (psmc_options *options)
{
  if (options)
    *options = (psmc_options){
      256, 8, 1000, .8, .5, .5, .05, 10, false, {CRPX_CLOCK_STRICT, 1, .3}, {CRPX_TREE_PRIOR_COALESCENT, 1, 0, 10, 1}};
}

void
psmc_result_free (psmc_result *result)
{
  if (!result) return;
  if (result->particles)
    for (size_t i = 0; i < result->count; i++)
      psmc_private_free_particle (result->particles + i);
  free (result->particles);
  free (result->weights);
  crpx_character_data_free ((crpx_character_data *) result->data);
  free (result);
}

crpx_status
psmc_anneal (const crpx_character_data *data, const crpx_substitution_model *model, const psmc_options *options,
             crpx_uniform_function uniform, void *state, psmc_result **out)
{ /* Each stage changes likelihood power, records its normaliser, then moves states at that power. */
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  psmc_evaluation context = {0};
  psmc_result *result = NULL;
  double *logs = NULL, *scores = NULL;
  crpx_status status = psmc_private_prepare (&context, data, model, options, uniform, state);
  if (status != CRPX_OK) goto done;
  status = psmc_private_result_new (&context, &result);
  if (status != CRPX_OK) goto done;
  size_t count = result->count;
  logs = malloc (count * sizeof (*logs));
  scores = malloc (count * sizeof (*scores));
  if (!logs || !scores) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  for (size_t i = 0; i < count; i++) {
    status = psmc_private_draw_prior (&context, result->particles + i);
    if (status != CRPX_OK) goto done;
  }
  while (result->beta < 1 && result->stages < options->max_stages) {
    for (size_t i = 0; i < count; i++)
      scores[i] = result->particles[i].log_likelihood;
    double next, normaliser;
    status = crpx_particle_temperature (result->weights, scores, count, result->beta, options->conditional_ess, &next);
    if (status != CRPX_OK) goto done;
    for (size_t i = 0; i < count; i++)
      logs[i] = log (result->weights[i]) + (next - result->beta) * scores[i];
    status = crpx_particle_normalise (logs, count, result->weights, &normaliser, &result->ess);
    if (status != CRPX_OK) goto done;
    result->log_evidence += normaliser;
    result->beta = next;
    result->stages++;
    if (result->ess < options->resample_ess * count) status = psmc_private_resample (&context, result);
    if (status == CRPX_OK) status = psmc_private_rejuvenate (&context, result);
    if (status != CRPX_OK) goto done;
  }
  if (result->beta != 1) {
    status = CRPX_ERROR_CONVERGENCE;
    goto done;
  }
  *out = result;
  result = NULL;
done:
  free (logs);
  free (scores);
  psmc_result_free (result);
  psmc_private_clear (&context);
  return status;
}

crpx_status
psmc_online (const psmc_result *previous, const crpx_character_data *data, crpx_uniform_function uniform, void *state,
             psmc_result **out)
{ /* A dimension-changing importance proposal, followed by full-posterior MH rejuvenation. */
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  crpx_status status = psmc_private_validate_online (previous, data);
  if (status != CRPX_OK) return status;
  psmc_evaluation context = {0};
  psmc_result *result = NULL;
  double *logs = NULL;
  status = psmc_private_prepare (&context, data, previous->model, &previous->options, uniform, state);
  if (status != CRPX_OK) goto done;
  status = psmc_private_result_new (&context, &result);
  if (status != CRPX_OK) goto done;
  result->beta = 1;
  result->stages = previous->stages + 1;
  logs = malloc (result->count * sizeof (*logs));
  if (!logs) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  double rate = result->options.branch_prior_rate;
  const char *label = crpx_character_data_name (data, crpx_character_data_taxa (data) - 1);
  for (size_t i = 0; i < result->count; i++) {
    const psmc_particle *old = previous->particles + i;
    size_t nodes = crpx_tree_node_count (old->tree), edge;
    double u, fraction, pendant, stem = 0, log_proposal = -log ((double) nodes);
    status = psmc_private_index_draw (&context, nodes, &edge);
    if (status == CRPX_OK) status = psmc_private_draw (&context, &fraction);
    if (status == CRPX_OK) status = psmc_private_draw (&context, &u);
    if (status != CRPX_OK) goto done;
    pendant = -log1p (-u) / rate;
    log_proposal += log (rate) - rate * pendant;
    if (edge == crpx_tree_root (old->tree)) {
      status = psmc_private_draw (&context, &u);
      if (status != CRPX_OK) goto done;
      stem = -log1p (-u) / rate;
      log_proposal += log (rate) - rate * stem;
    } else {
      double length;
      status = crpx_tree_branch_length (old->tree, (crpx_node_id) edge, &length);
      if (status != CRPX_OK) goto done;
      if (length <= 0) {
        status = CRPX_ERROR_INVALID_ARGUMENT;
        goto done;
      }
      log_proposal -= log (length);
    }
    psmc_particle *particle = result->particles + i;
    status = crpx_tree_insert_tip (old->tree, (crpx_node_id) edge, fraction, pendant, stem, label, &particle->tree);
    if (status == CRPX_OK) status = psmc_private_evaluate (&context, particle);
    if (status != CRPX_OK) goto done;
    logs[i] = log (previous->weights[i]) + particle->log_prior - old->log_prior - log_proposal +
              particle->log_likelihood - old->log_likelihood;
  }
  double increment;
  status = crpx_particle_normalise (logs, result->count, result->weights, &increment, &result->ess);
  if (status != CRPX_OK) goto done;
  result->log_evidence = previous->log_evidence + increment;
  if (result->ess < result->options.resample_ess * result->count) status = psmc_private_resample (&context, result);
  if (status == CRPX_OK) status = psmc_private_rejuvenate (&context, result);
  if (status == CRPX_OK) {
    *out = result;
    result = NULL;
  }
done:
  free (logs);
  psmc_result_free (result);
  psmc_private_clear (&context);
  return status;
}

crpx_status
psmc_substitution_tree (const psmc_result *result, size_t particle, crpx_tree **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!result || particle >= result->count) return CRPX_ERROR_INVALID_ARGUMENT;
  const psmc_particle *input = result->particles + particle;
  crpx_status status = crpx_tree_clone (input->tree, out);
  if (status != CRPX_OK || !result->options.use_clock) return status;
  size_t nodes = crpx_tree_node_count (input->tree);
  double *lengths = malloc (nodes * sizeof (*lengths));
  if (!lengths)
    status = CRPX_ERROR_ALLOCATION;
  else
    status = crpx_clock_lengths (input->tree, input->ages, input->rates, 1, lengths, nodes);
  for (size_t i = 0; i < nodes && status == CRPX_OK; i++)
    if (i != crpx_tree_root (*out)) status = crpx_tree_set_branch_length (*out, (crpx_node_id) i, lengths[i]);
  free (lengths);
  if (status != CRPX_OK) {
    crpx_tree_free (*out);
    *out = NULL;
  }
  return status;
}
