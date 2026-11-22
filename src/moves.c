/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
#include "smc_internal.h"
#include <particle_weights.h>
#include <topology_randomise.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
static crpx_status local_proposal (psmc_evaluation *context, psmc_particle *proposal, double *reverse_forward,
                                   bool *valid);
static bool eligible_node (const crpx_tree *tree, crpx_node_id node, unsigned kind, const psmc_evaluation *context);
static crpx_status set_times (psmc_particle *particle);

crpx_status
psmc_private_resample (psmc_evaluation *context, psmc_result *result)
{
  size_t *indices = malloc (result->count * sizeof (*indices));
  psmc_particle *copies = calloc (result->count, sizeof (*copies));
  size_t *ancestors = malloc (result->count * sizeof (*ancestors));
  crpx_status status = CRPX_OK;
  if (!ancestors || !indices || !copies) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  status = crpx_particle_resample (result->weights, result->count, result->count, false, context->uniform,
                                   context->state, indices);
  for (size_t i = 0; i < result->count && status == CRPX_OK; i++)
    status = psmc_private_copy_particle (result->particles + indices[i], copies + i);
  if (status == CRPX_OK) {
    for (size_t i = 0; i < result->count; i++) ancestors[i] = result->ancestor_ids[indices[i]];
    memcpy (result->ancestor_ids, ancestors, result->count * sizeof (*ancestors));
    result->surviving_ancestors = 0;
    memset (indices, 0, result->count * sizeof (*indices));
    for (size_t i = 0; i < result->count; i++)
      if (!indices[ancestors[i]]++) result->surviving_ancestors++;
    for (size_t i = 0; i < result->count; i++) {
      psmc_private_free_particle (result->particles + i);
      result->weights[i] = 1.0 / result->count;
    }
    free (result->particles);
    result->particles = copies;
    copies = NULL;
    result->ess = (double) result->count;
  }
done:
  if (copies) {
    for (size_t i = 0; i < result->count; i++)
      psmc_private_free_particle (copies + i);
    free (copies);
  }
  free (indices);
  free (ancestors);
  return status;
}

crpx_status
psmc_private_rejuvenate (psmc_evaluation *context, psmc_result *result)
{
  return psmc_private_execute (context, result, false);
}

crpx_status
psmc_private_rejuvenate_one (psmc_evaluation *context, psmc_result *result, size_t i)
{ /* Reweighting selects the target; these MH moves restore diversity at that target. */
  if (result->weights[i] == 0) return CRPX_OK;
    for (size_t step = 0; step < result->options.mcmc_steps; step++) {
      psmc_particle proposal = {0}, *current = result->particles + i;
      double u, reverse_forward = 0;
      crpx_status status = psmc_private_draw (context, &u);
      if (status != CRPX_OK) return status;
      bool independent = u < result->options.independence_probability, valid = true;
      if (independent)
        status = psmc_private_draw_prior (context, &proposal);
      else {
        status = psmc_private_copy_particle (current, &proposal);
        if (status == CRPX_OK) status = local_proposal (context, &proposal, &reverse_forward, &valid);
        if (status == CRPX_OK && valid) status = psmc_private_evaluate (context, &proposal);
      }
      if (status != CRPX_OK) {
        psmc_private_free_particle (&proposal);
        return status;
      }
      result->proposals++;
      if (!valid) {
        psmc_private_free_particle (&proposal);
        continue;
      }
      double difference = result->beta * (proposal.log_likelihood - current->log_likelihood);
      if (!independent) difference += proposal.log_prior - current->log_prior;
      status = psmc_private_draw (context, &u);
      if (status != CRPX_OK) {
        psmc_private_free_particle (&proposal);
        return status;
      }
      bool accepted;
      status = crpx_metropolis_accept (difference, independent ? 0 : reverse_forward, log (u), &accepted);
      if (status != CRPX_OK) {
        psmc_private_free_particle (&proposal);
        return status;
      }
      if (accepted) {
        psmc_private_free_particle (current);
        *current = proposal;
        result->accepted++;
      } else
        psmc_private_free_particle (&proposal);
    }
  return CRPX_OK;
}
static crpx_status
local_proposal (psmc_evaluation *context, psmc_particle *proposal, double *reverse_forward, bool *valid)
{ /* Selection probabilities are topology-independent on labelled binary trees. */
  const psmc_options *options = context->options;
  double u;
  crpx_status status = psmc_private_draw (context, &u);
  if (status != CRPX_OK) return status;
  if ((options->infer || options->substitution != PSMC_MODEL_FIXED) && u < .25) {
    status = psmc_private_parameter_proposal (context, proposal, reverse_forward, valid);
    if (status == CRPX_OK && *valid && options->use_clock) {
      if (options->clock.kind == CRPX_CLOCK_STRICT)
        for (size_t i = 0; i < context->nodes; i++) proposal->rates[i] = proposal->parameters.clock_mean;
      proposal->rates[crpx_tree_root (proposal->tree)] = proposal->parameters.clock_mean;
    }
    return status;
  }
  status = psmc_private_draw (context, &u);
  if (status != CRPX_OK) return status;
  unsigned kind = (unsigned) (u * (options->use_clock ? 3 : 2));
  if (options->use_clock && options->clock.kind == CRPX_CLOCK_STRICT && kind == 2) kind = 1;
  size_t eligible = 0, selected;
  for (size_t node = 0; node < context->nodes; node++)
    eligible += eligible_node (proposal->tree, (crpx_node_id) node, kind, context);
  if (!eligible) {
    *valid = false;
    return CRPX_OK;
  }
  status = psmc_private_index_draw (context, eligible, &selected);
  if (status != CRPX_OK) return status;
  crpx_node_id node = 0;
  for (; node < context->nodes; node++)
    if (eligible_node (proposal->tree, node, kind, context)) {
      if (selected == 0) break;
      selected--;
    }
  status = psmc_private_draw (context, &u);
  if (status != CRPX_OK) return status;
  if (kind == 0)
    status = crpx_tree_nni (proposal->tree, node, u < .5 ? 0 : 1);
  else if (options->use_clock && kind == 1)
    proposal->ages[node] += (2 * u - 1) * options->proposal_width;
  else {
    double change = (2 * u - 1) * options->proposal_width, value;
    *reverse_forward = change;
    if (options->use_clock)
      value = proposal->rates[node];
    else {
      status = crpx_tree_branch_length (proposal->tree, node, &value);
      if (status != CRPX_OK) return status;
    }
    value *= exp (change);
    if (!isfinite (value) || value <= 0) {
      *valid = false;
      return CRPX_OK;
    }
    if (options->use_clock)
      proposal->rates[node] = value;
    else
      status = crpx_tree_set_branch_length (proposal->tree, node, value);
  }
  if (status == CRPX_OK && options->use_clock) {
    status = set_times (proposal);
    if (status == CRPX_ERROR_INVALID_ARGUMENT) {
      *valid = false;
      status = CRPX_OK;
    }
  }
  return status;
}
static bool
eligible_node (const crpx_tree *tree, crpx_node_id node, unsigned kind, const psmc_evaluation *context)
{
  bool clock = context->options->use_clock;
  psmc_root_mode mode = psmc_private_root_mode (context->options);
  if (mode == PSMC_UNROOTED && crpx_tree_first_child (tree, node) == CRPX_NODE_NONE &&
      !strcmp (crpx_tree_label (tree, node), crpx_character_data_name (context->data, 0))) return false;
  if (kind == 0 && (mode == PSMC_UNROOTED || mode == PSMC_ROOT_FIXED) &&
      crpx_tree_parent (tree, node) == crpx_tree_root (tree)) return false;
  bool internal = crpx_tree_first_child (tree, node) != CRPX_NODE_NONE;
  if (kind == 0) return internal && node != crpx_tree_root (tree);
  if (kind == 1 && clock) return internal;
  return node != crpx_tree_root (tree);
}

static crpx_status
set_times (psmc_particle *particle)
{
  size_t nodes = crpx_tree_node_count (particle->tree);
  for (size_t i = 0; i < nodes; i++) {
    crpx_node_id parent = crpx_tree_parent (particle->tree, (crpx_node_id) i);
    if (parent != CRPX_NODE_NONE && particle->ages[parent] <= particle->ages[i]) return CRPX_ERROR_INVALID_ARGUMENT;
  }
  for (size_t i = 0; i < nodes; i++) {
    crpx_node_id parent = crpx_tree_parent (particle->tree, (crpx_node_id) i);
    if (parent != CRPX_NODE_NONE) {
      crpx_status status =
        crpx_tree_set_branch_length (particle->tree, (crpx_node_id) i, particle->ages[parent] - particle->ages[i]);
      if (status != CRPX_OK) return status;
    }
  }
  return CRPX_OK;
}
