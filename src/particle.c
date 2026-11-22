/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
#include "smc_internal.h"
#include <particle_weights.h>
#include <topology_randomise.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

crpx_status
psmc_private_prepare (psmc_evaluation *context, const crpx_character_data *data, const crpx_substitution_model *model,
                      const psmc_options *options, crpx_uniform_function uniform, void *state)
{
  if (!data || !model || !options || !uniform || !options->particles || !options->max_stages ||
      crpx_character_data_taxa (data) < 2 ||
      crpx_character_data_states (data) != crpx_substitution_model_states (model) ||
      !isfinite (options->conditional_ess) || options->conditional_ess <= 0 || options->conditional_ess >= 1 ||
      !isfinite (options->resample_ess) || options->resample_ess <= 0 || options->resample_ess > 1 ||
      !isfinite (options->proposal_width) || options->proposal_width <= 0 ||
      !isfinite (options->independence_probability) || options->independence_probability < 0 ||
      options->independence_probability > 1 || !isfinite (options->branch_prior_rate) ||
      options->branch_prior_rate <= 0 || (options->use_clock && options->tree_prior.kind == CRPX_TREE_PRIOR_UNIFORM))
    return CRPX_ERROR_INVALID_ARGUMENT;
  size_t taxa = crpx_character_data_taxa (data);
  if (taxa > UINT32_MAX / 2 || options->particles > SIZE_MAX / sizeof (psmc_particle)) return CRPX_ERROR_OVERFLOW;
  context->data = data;
  context->model = model;
  context->options = options;
  context->uniform = uniform;
  context->state = state;
  context->nodes = 2 * taxa - 1;
  if (context->nodes > SIZE_MAX / sizeof (double)) return CRPX_ERROR_OVERFLOW;
  context->lengths = malloc (context->nodes * sizeof (double));
  if (!context->lengths) return CRPX_ERROR_ALLOCATION;
  crpx_status status = crpx_likelihood_workspace_new (&context->likelihood);
  if (status == CRPX_OK) status = crpx_tree_workspace_new (CRPX_TREE_WORKSPACE_TRAVERSAL, &context->traversal);
  return status;
}

void
psmc_private_clear (psmc_evaluation *context)
{
  free (context->lengths);
  crpx_likelihood_workspace_free (context->likelihood);
  crpx_tree_workspace_free (context->traversal);
}

crpx_status
psmc_private_result_new (psmc_evaluation *context, psmc_result **out)
{
  psmc_result *result = calloc (1, sizeof (*result));
  if (!result) return CRPX_ERROR_ALLOCATION;
  result->count = context->options->particles;
  result->options = *context->options;
  result->model = context->model;
  result->particles = calloc (result->count, sizeof (*result->particles));
  result->weights = malloc (result->count * sizeof (*result->weights));
  crpx_status status = crpx_character_data_retain (context->data);
  if (status == CRPX_OK) result->data = context->data;
  if (!result->particles || !result->weights) status = CRPX_ERROR_ALLOCATION;
  if (status != CRPX_OK) {
    psmc_result_free (result);
    return status;
  }
  for (size_t i = 0; i < result->count; i++)
    result->weights[i] = 1.0 / result->count;
  result->ess = (double) result->count;
  *out = result;
  return CRPX_OK;
}

crpx_status
psmc_private_draw_prior (psmc_evaluation *context, psmc_particle *particle)
{
  size_t taxa = crpx_character_data_taxa (context->data), nodes = context->nodes;
  const char **names = malloc (taxa * sizeof (*names));
  if (!names) return CRPX_ERROR_ALLOCATION;
  for (size_t i = 0; i < taxa; i++)
    names[i] = crpx_character_data_name (context->data, i);
  crpx_tree_prior uniform_prior = {CRPX_TREE_PRIOR_UNIFORM, 0, 0, 0, 0};
  crpx_status status =
    crpx_tree_prior_sample (context->options->use_clock ? &context->options->tree_prior : &uniform_prior, taxa, names,
                            context->uniform, context->state, &particle->tree);
  free (names);
  if (status != CRPX_OK) return status;
  if (context->options->use_clock) {
    particle->ages = malloc (nodes * sizeof (double));
    particle->rates = malloc (nodes * sizeof (double));
    if (!particle->ages || !particle->rates) return CRPX_ERROR_ALLOCATION;
    status = crpx_tree_ages (particle->tree, NULL, 1e-9, context->traversal, particle->ages, nodes);
    if (status == CRPX_OK)
      status = crpx_clock_sample (&context->options->clock, particle->tree, particle->ages, context->uniform,
                                  context->state, context->traversal, particle->rates, nodes);
  } else
    for (size_t i = 0; i < nodes && status == CRPX_OK; i++)
      if (i != crpx_tree_root (particle->tree)) {
        double u;
        status = psmc_private_draw (context, &u);
        if (status == CRPX_OK)
          status = crpx_tree_set_branch_length (particle->tree, (crpx_node_id) i,
                                                -log1p (-u) / context->options->branch_prior_rate);
      }
  if (status == CRPX_OK) status = psmc_private_evaluate (context, particle);
  return status;
}

crpx_status
psmc_private_evaluate (psmc_evaluation *context, psmc_particle *particle)
{
  crpx_likelihood_options options;
  crpx_likelihood_options_init (&options);
  options.retain_partials = false;
  crpx_status status;
  if (context->options->use_clock) {
    status = crpx_tree_prior_log_density (&context->options->tree_prior, particle->tree, &particle->log_prior);
    double clock_density;
    if (status == CRPX_OK)
      status = crpx_clock_log_density (&context->options->clock, particle->tree, particle->ages, particle->rates,
                                       context->traversal, context->nodes, &clock_density);
    if (status != CRPX_OK) return status;
    particle->log_prior += clock_density;
    status = crpx_clock_lengths (particle->tree, particle->ages, particle->rates, 1, context->lengths, context->nodes);
    if (status != CRPX_OK) return status;
    options.edge_values = context->lengths;
    options.edge_count = context->nodes;
  } else {
    crpx_tree_prior uniform_prior = {CRPX_TREE_PRIOR_UNIFORM, 0, 0, 0, 0};
    status = crpx_tree_prior_log_density (&uniform_prior, particle->tree, &particle->log_prior);
    if (status != CRPX_OK) return status;
    for (size_t i = 0; i < context->nodes; i++)
      if (i != crpx_tree_root (particle->tree)) {
        double length;
        status = crpx_tree_branch_length (particle->tree, (crpx_node_id) i, &length);
        if (status != CRPX_OK) return status;
        if (!isfinite (length) || length <= 0) return CRPX_ERROR_INVALID_ARGUMENT;
        particle->log_prior += log (context->options->branch_prior_rate) - context->options->branch_prior_rate * length;
      }
  }
  return crpx_likelihood_evaluate (particle->tree, context->data, context->model, &options, context->likelihood,
                                   &particle->log_likelihood);
}

crpx_status
psmc_private_copy_particle (const psmc_particle *source, psmc_particle *copy)
{
  memset (copy, 0, sizeof (*copy));
  crpx_status status = crpx_tree_clone (source->tree, &copy->tree);
  if (status != CRPX_OK) return status;
  copy->log_likelihood = source->log_likelihood;
  copy->log_prior = source->log_prior;
  if (source->ages) {
    size_t bytes = crpx_tree_node_count (source->tree) * sizeof (double);
    copy->ages = malloc (bytes);
    copy->rates = malloc (bytes);
    if (!copy->ages || !copy->rates) {
      psmc_private_free_particle (copy);
      return CRPX_ERROR_ALLOCATION;
    }
    memcpy (copy->ages, source->ages, bytes);
    memcpy (copy->rates, source->rates, bytes);
  }
  return CRPX_OK;
}

void
psmc_private_free_particle (psmc_particle *particle)
{
  crpx_tree_free (particle->tree);
  free (particle->ages);
  free (particle->rates);
  memset (particle, 0, sizeof (*particle));
}

crpx_status
psmc_private_draw (psmc_evaluation *context, double *value)
{
  crpx_status status = crpx_sample_uniform (context->uniform, context->state, value);
  if (status == CRPX_OK && !*value) *value = 0x1p-54;
  return status;
}

crpx_status
psmc_private_index_draw (psmc_evaluation *context, size_t count, size_t *index)
{
  double u;
  crpx_status status = psmc_private_draw (context, &u);
  if (status == CRPX_OK) *index = (size_t) ((long double) u * count);
  return status;
}
