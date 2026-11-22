/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
#include "phylosmc.h"
#include <global/random_number.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#define OK(x) assert ((x) == CRPX_OK)
static double uniform (void *state) { return crpx_rng_stream_double (state); }
static void
same (const psmc_result *a, const psmc_result *b)
{
  assert (a->log_evidence == b->log_evidence && a->accepted == b->accepted && a->stages == b->stages);
  for (size_t i = 0; i < a->count; i++) {
    char *x, *y;
    OK (crpx_tree_to_newick (a->particles[i].tree, &x));
    OK (crpx_tree_to_newick (b->particles[i].tree, &y));
    assert (!strcmp (x, y) && a->weights[i] == b->weights[i]);
    assert (!memcmp (&a->particles[i].parameters, &b->particles[i].parameters, sizeof (psmc_parameters)));
    if (a->particles[i].ages) {
      size_t bytes = crpx_tree_node_count (a->particles[i].tree) * sizeof (double);
      assert (!memcmp (a->particles[i].ages, b->particles[i].ages, bytes));
      assert (!memcmp (a->particles[i].rates, b->particles[i].rates, bytes));
    }
    crpx_tree_string_free (x); crpx_tree_string_free (y);
  }
}
int
main (void)
{
  const char *names[] = {"A","B","C","D"};
  uint8_t unknown[] = {255,255,255,255};
  crpx_character_data *data, *two;
  crpx_substitution_model *model;
  crpx_rng_stream stream;
  crpx_rng_stream_init (&stream, 555, 0);
  OK (crpx_character_data_from_binary (names, 4, 1, unknown, NULL, true, &data));
  OK (crpx_substitution_model_poisson (2, &model));
  psmc_options options;
  psmc_options_init (&options);
  options.particles = 6000; options.mcmc_steps = 4;
  psmc_result *result, *extended;
  OK (psmc_anneal (data, model, &options, uniform, &stream, &result));
  double mean = 0, cherry = 0;
  for (size_t i = 0; i < result->count; i++) {
    crpx_tree *tree = result->particles[i].tree;
    crpx_node_id a = CRPX_NODE_NONE, b = CRPX_NODE_NONE, c = CRPX_NODE_NONE;
    for (crpx_node_id v = 0; v < crpx_tree_node_count (tree); v++) {
      const char *label = crpx_tree_label (tree, v);
      if (label && !strcmp (label, "A")) a = v;
      if (label && !strcmp (label, "B")) b = v;
      if (label && !strcmp (label, "C")) c = v;
      if (v != crpx_tree_root (tree)) {
        double length;
        OK (crpx_tree_branch_length (tree, v, &length));
        mean += result->weights[i] * length / 5;
      }
    }
    double zero;
    OK (crpx_tree_branch_length (tree, a, &zero));
    assert (zero == 0 && crpx_tree_parent (tree, a) == crpx_tree_root (tree));
    cherry += result->weights[i] * (crpx_tree_parent (tree, b) == crpx_tree_parent (tree, c));
  }
  assert (fabs (mean - .1) < .006 && fabs (cherry - 1./3) < .035);
  psmc_result_free (result);
  uint8_t observed[] = {0,1};
  OK (crpx_character_data_from_binary (names, 2, 1, observed, NULL, true, &two));
  options.particles = 20000; options.mcmc_steps = 0;
  OK (psmc_anneal (two, model, &options, uniform, &stream, &result));
  assert (fabs (exp (result->log_evidence) - .25 * (1 - 10./12)) < .001);
  uint8_t obs3[] = {0,1,255};
  crpx_character_data *three;
  OK (crpx_character_data_from_binary (names, 3, 1, obs3, NULL, true, &three));
  OK (psmc_online (result, three, uniform, &stream, &extended));
  assert (fabs (exp (extended->log_evidence) - .25 * (1 - pow (10./12, 2))) < .002);
  psmc_result_free (result); psmc_result_free (extended);
  crpx_character_data_free (two); crpx_character_data_free (three);
  bool side[] = {true,true,false,false};
  options.root_mode = PSMC_ROOT_FIXED; options.root_side = side;
  options.particles = 1200; options.mcmc_steps = 3;
  for (unsigned clock = 0; clock < 2; clock++) {
    options.use_clock = clock;
    OK (psmc_anneal (data, model, &options, uniform, &stream, &result));
    assert (fabs (result->log_evidence) < 1e-12);
    for (size_t i = 0; i < result->count; i++) {
      bool actual[4];
      OK (psmc_root_split (result->particles[i].tree, data, actual));
      for (size_t j = 0; j < 4; j++) assert ((actual[j] != actual[0]) == (side[j] != side[0]));
    }
    psmc_result_free (result);
  }
  options.use_clock = false;
  options.particles = 3000; options.resample_ess = .01;
  OK (psmc_data_prefix (data, 3, &three));
  OK (psmc_anneal (three, model, &options, uniform, &stream, &result));
  OK (psmc_online_split (result, data, side, uniform, &stream, &extended));
  assert (fabs (extended->log_evidence) < .08);
  psmc_result_free (result); psmc_result_free (extended); crpx_character_data_free (three);
  const char *all_names[] = {"A","B","C","D","E","F"};
  uint8_t all_unknown[] = {255,255,255,255,255,255};
  bool all_side[] = {true,true,false,false,false,false};
  crpx_character_data *all;
  OK (crpx_character_data_from_binary (all_names,6,1,all_unknown,NULL,true,&all));
  OK (psmc_data_prefix (all,3,&three));
  OK (psmc_anneal (three,model,&options,uniform,&stream,&result));
  crpx_character_data_free (three);
  for (size_t n=4;n<=6;n++) {
    OK (psmc_data_prefix (all,n,&three));
    OK (psmc_online_split (result,three,all_side,uniform,&stream,&extended));
    assert (isfinite (extended->log_evidence));
    psmc_result_free (result); result=extended;
    crpx_character_data_free (three);
  }
  psmc_result_free (result); crpx_character_data_free (all);
  crpx_substitution_model_free (model); crpx_character_data_free (data);
  /* Joint parameter priors survive zero-information likelihoods and rejuvenation. */
  double tips[16]; for (size_t i = 0; i < 16; i++) tips[i] = 1;
  OK (crpx_character_data_from_tip_likelihoods (names, 4, 1, 4, tips, NULL, true, &data));
  OK (crpx_substitution_model_poisson (4, &model));
  psmc_options_init (&options);
  options.particles = 2500; options.mcmc_steps = 8;
  options.substitution = PSMC_MODEL_GTR; options.infer = PSMC_INFER_BRANCH_RATE;
  OK (psmc_anneal (data, model, &options, uniform, &stream, &result));
  double logmean = 0, frequency = 0, exchange = 0;
  for (size_t i = 0; i < result->count; i++) {
    const psmc_parameters *p = &result->particles[i].parameters;
    logmean += log (p->branch_rate) / result->count;
    frequency += p->frequencies[0] / result->count;
    exchange += p->exchangeabilities[0] / result->count;
  }
  assert (fabs (logmean - log (10)) < .08 && fabs (frequency - .25) < .02 && fabs (exchange - 1./6) < .02);
  assert (fabs (result->log_evidence) < 1e-10);
  psmc_result_free (result);
#ifndef _WIN32
  /* Complete populations, including clock coordinates, reproduce across process counts. */
  options.particles = 24; options.mcmc_steps = 5;
  options.substitution = PSMC_MODEL_HKY;
  for (unsigned clock = 0; clock < 5; clock++) {
    options.use_clock = clock != 0;
    options.clock.kind = clock ? (crpx_clock_kind)(clock - 1) : CRPX_CLOCK_STRICT;
    options.infer = clock ? PSMC_INFER_COALESCENT_SCALE : PSMC_INFER_BRANCH_RATE;
    if (clock == 2 || clock == 4) options.infer |= PSMC_INFER_CLOCK_DISPERSION;
    crpx_rng_stream a, b;
    crpx_rng_stream_init (&a, 123, 0); b = a;
    options.workers = 1;
    OK (psmc_anneal (data, model, &options, uniform, &a, &result));
    options.workers = 3;
    OK (psmc_anneal (data, model, &options, uniform, &b, &extended));
    same (result, extended);
    psmc_result_free (result); psmc_result_free (extended);
  }
#endif
  options.use_clock = false; options.root_mode = PSMC_UNROOTED;
  options.infer = PSMC_INFER_BRANCH_RATE; options.workers = 0;
  result = (void *)1;
  assert (psmc_anneal (data, model, &options, uniform, &stream, &result) == CRPX_ERROR_INVALID_ARGUMENT && !result);
  options.workers = 1; options.use_clock = true;
  assert (psmc_anneal (data, model, &options, uniform, &stream, &result) == CRPX_ERROR_INVALID_ARGUMENT);
  options.use_clock = false; options.priors.log_sd[0] = 0;
  assert (psmc_anneal (data, model, &options, uniform, &stream, &result) == CRPX_ERROR_INVALID_ARGUMENT);
  crpx_substitution_model_free (model); crpx_character_data_free (data);
  return 0;
}
