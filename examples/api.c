/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
#include <phylosmc.h>
#include <global/random_number.h>
#include <stdio.h>
static double
uniform (void *stream)
{
  return crpx_rng_stream_double (stream);
}
int
main (void)
{
  const char *names[] = {"A", "B", "C"};
  uint8_t traits[] = {0, 0, 1, 0, 1, 1};
  crpx_character_data *data = NULL;
  crpx_substitution_model *model = NULL;
  psmc_result *result = NULL;
  psmc_options options;
  psmc_options_init (&options);
  options.use_clock = true;
  options.clock = (crpx_clock_model){CRPX_CLOCK_LOGNORMAL, .1, .4};
  options.tree_prior.coalescent_scale = 10;
  crpx_rng_stream stream;
  crpx_rng_stream_init (&stream, 42, 0);
  crpx_status status = crpx_character_data_from_binary (names, 3, 2, traits, NULL, true, &data);
  if (status == CRPX_OK) status = crpx_substitution_model_poisson (2, &model);
  if (status == CRPX_OK) status = psmc_anneal (data, model, &options, uniform, &stream, &result);
  if (status == CRPX_OK)
    printf ("Log evidence: %.9g\n", result->log_evidence);
  else
    fprintf (stderr, "%s\n", crpx_status_string (status));
  psmc_result_free (result);
  crpx_substitution_model_free (model);
  crpx_character_data_free (data);
  return status == CRPX_OK ? 0 : 1;
}
