/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "phylosmc.h"
#include <global/random_number.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#define OK(call) assert ((call) == CRPX_OK)
static double
uniform (void *state)
{
  return crpx_rng_stream_double (state);
}
static double
invalid_uniform (void *state)
{
  (void) state;
  return 1;
}
int
main (void)
{
  const char *names[] = {"A", "B", "C", "D"};
  uint8_t unknown[] = {255, 255, 255, 255};
  crpx_character_data *data, *prefix;
  crpx_substitution_model *model;
  crpx_rng_stream stream;
  crpx_rng_stream_init (&stream, 1234, 0);
  OK (crpx_character_data_from_binary (names, 4, 1, unknown, NULL, true, &data));
  OK (crpx_substitution_model_poisson (2, &model));
  psmc_options options;
  psmc_options_init (&options);
  options.particles = 3000;
  options.mcmc_steps = 3;
  options.independence_probability = .15;
  OK (psmc_data_prefix (data, 3, &prefix));
  psmc_result *result, *updated;
  OK (psmc_anneal (prefix, model, &options, uniform, &stream, &result));
  assert (fabs (result->log_evidence) < 1e-12 && result->beta == 1 && result->stages == 1);
  /* With unknown observations, prior branch means must survive MH and online weighting. */
  double mean = 0;
  for (size_t i = 0; i < result->count; i++) {
    size_t nodes = crpx_tree_node_count (result->particles[i].tree);
    for (size_t node = 0; node < nodes; node++)
      if (node != crpx_tree_root (result->particles[i].tree)) {
        double length;
        OK (crpx_tree_branch_length (result->particles[i].tree, (crpx_node_id) node, &length));
        mean += result->weights[i] * length / (nodes - 1);
      }
  }
  assert (fabs (mean - .1) < .007);
  OK (psmc_online (result, data, uniform, &stream, &updated));
  assert (fabs (updated->log_evidence) < .08 && updated->count == result->count);
  size_t balanced = 0;
  double weighted_balanced = 0, updated_mean = 0;
  for (size_t i = 0; i < updated->count; i++) {
    crpx_tree *tree = updated->particles[i].tree;
    assert (crpx_tree_leaf_count (tree) == 4);
    crpx_node_id first = crpx_tree_first_child (tree, crpx_tree_root (tree));
    crpx_node_id second = crpx_tree_next_sibling (tree, first);
    bool shape =
      crpx_tree_first_child (tree, first) != CRPX_NODE_NONE && crpx_tree_first_child (tree, second) != CRPX_NODE_NONE;
    balanced += shape;
    weighted_balanced += updated->weights[i] * shape;
    for (size_t node = 0; node < crpx_tree_node_count (tree); node++)
      if (node != crpx_tree_root (tree)) {
        double length;
        OK (crpx_tree_branch_length (tree, (crpx_node_id) node, &length));
        updated_mean += updated->weights[i] * length / 6;
      }
  }
  (void) balanced;
  assert (fabs (weighted_balanced - .2) < .04 && fabs (updated_mean - .1) < .008);
  psmc_result_free (updated);
  psmc_result_free (result);
  crpx_character_data_free (prefix);
  crpx_character_data_free (data);
  /* Two binary observations: integrate Exp(lambda) edges analytically. */
  uint8_t observed[] = {0, 1};
  OK (crpx_character_data_from_binary (names, 2, 1, observed, NULL, true, &data));
  options.particles = 16000;
  options.mcmc_steps = 0;
  OK (psmc_anneal (data, model, &options, uniform, &stream, &result));
  double exact = .25 * (1 - pow (10.0 / 12, 2));
  assert (fabs (exp (result->log_evidence) - exact) < .0015);
  uint8_t third[] = {0, 1, 255};
  crpx_character_data *extended;
  OK (crpx_character_data_from_binary (names, 3, 1, third, NULL, true, &extended));
  OK (psmc_online (result, extended, uniform, &stream, &updated));
  double r = 10.0 / 12;
  double exact_three = .25 * (1 - (r * r + 2 * r * r * r) / 3);
  assert (fabs (exp (updated->log_evidence) - exact_three) < .003);
  psmc_result_free (updated);
  crpx_character_data_free (extended);
  psmc_result_free (result);
  crpx_character_data_free (data);
  /* Identical streams give identical complete runs; failed calls have no result. */
  options.particles = 32;
  options.mcmc_steps = 3;
  crpx_rng_stream first_stream, second_stream;
  crpx_rng_stream_init (&first_stream, 567, 0);
  crpx_rng_stream_init (&second_stream, 567, 0);
  uint8_t repeat[] = {0, 1};
  OK (crpx_character_data_from_binary (names, 2, 1, repeat, NULL, true, &data));
  OK (psmc_anneal (data, model, &options, uniform, &first_stream, &result));
  OK (psmc_anneal (data, model, &options, uniform, &second_stream, &updated));
  assert (result->log_evidence == updated->log_evidence && result->accepted == updated->accepted);
  for (size_t i = 0; i < result->count; i++) {
    char *a, *b;
    OK (crpx_tree_to_newick (result->particles[i].tree, &a));
    OK (crpx_tree_to_newick (updated->particles[i].tree, &b));
    assert (!strcmp (a, b) && result->weights[i] == updated->weights[i]);
    crpx_tree_string_free (a);
    crpx_tree_string_free (b);
  }
  psmc_result_free (result);
  psmc_result_free (updated);
  result = (void *) 1;
  assert (psmc_anneal (data, model, &options, invalid_uniform, NULL, &result) == CRPX_ERROR_INVALID_ARGUMENT &&
          !result);
  crpx_character_data_free (data);
  /* Every clock family exercises posterior age/rate moves and substitution conversion. */
  uint8_t traits[] = {0, 0, 1, 1, 0, 1, 0, 1};
  OK (crpx_character_data_from_binary (names, 4, 2, traits, NULL, true, &data));
  for (unsigned kind = 0; kind < 4; kind++) {
    options.particles = 64;
    options.mcmc_steps = 8;
    options.use_clock = true;
    options.clock.kind = (crpx_clock_kind) kind;
    OK (psmc_anneal (data, model, &options, uniform, &stream, &result));
    assert (isfinite (result->log_evidence) && result->beta == 1 && result->accepted > 0);
    crpx_tree *substitutions;
    OK (psmc_substitution_tree (result, 0, &substitutions));
    OK (crpx_tree_validate_policy (substitutions, CRPX_TREE_REQUIRE_BINARY | CRPX_TREE_REQUIRE_COMPLETE_LENGTHS |
                                                    CRPX_TREE_REQUIRE_POSITIVE_LENGTHS));
    crpx_tree_free (substitutions);
    updated = (void *) 1;
    assert (psmc_online (result, data, uniform, &stream, &updated) == CRPX_ERROR_UNSUPPORTED && !updated);
    psmc_result_free (result);
  }
  crpx_character_data_free (data);
  crpx_substitution_model_free (model);
  return 0;
}
