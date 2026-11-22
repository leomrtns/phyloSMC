/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
#include "phylosmc.h"
#include <global/random_number.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
static double uniform (void *state);
static int usage (void);
static bool number (const char *text, size_t *value);
int
main (int argc, char **argv)
{ /* A small DNA front end; the C API also accepts traits and user-supplied generators. */
  if (argc < 2) return usage ();
  if (!strcmp (argv[1], "--help")) {
    (void) usage ();
    return 0;
  }
  psmc_options options;
  psmc_options_init (&options);
  size_t seed = 1, initial = 0;
  for (int i = 2; i < argc; i++) {
    if (i + 1 == argc) return usage ();
    const char *name = argv[i++], *value = argv[i];
    if (!strcmp (name, "--clock")) {
      options.use_clock = strcmp (value, "none") != 0;
      if (!strcmp (value, "strict"))
        options.clock.kind = CRPX_CLOCK_STRICT;
      else if (!strcmp (value, "lognormal"))
        options.clock.kind = CRPX_CLOCK_LOGNORMAL;
      else if (!strcmp (value, "exponential"))
        options.clock.kind = CRPX_CLOCK_EXPONENTIAL;
      else if (!strcmp (value, "autocorrelated"))
        options.clock.kind = CRPX_CLOCK_AUTOCORRELATED;
      else if (strcmp (value, "none"))
        return usage ();
    } else {
      size_t parsed;
      if (!number (value, &parsed)) return usage ();
      if (!strcmp (name, "--particles"))
        options.particles = parsed;
      else if (!strcmp (name, "--moves"))
        options.mcmc_steps = parsed;
      else if (!strcmp (name, "--seed"))
        seed = parsed;
      else if (!strcmp (name, "--initial-taxa"))
        initial = parsed;
      else
        return usage ();
    }
  }
  crpx_alignment *alignment = NULL;
  crpx_character_data *data = NULL, *subset = NULL;
  crpx_substitution_model *model = NULL;
  psmc_result *result = NULL;
  crpx_rng_stream stream;
  crpx_rng_stream_init (&stream, seed, 0);
  crpx_status status = crpx_alignment_read (argv[1], true, &alignment);
  if (status == CRPX_OK) status = crpx_character_data_from_alignment_dna (alignment, true, &data);
  if (status == CRPX_OK) status = crpx_substitution_model_poisson (4, &model);
  size_t taxa = crpx_character_data_taxa (data);
  if (status == CRPX_OK && initial && (initial < 2 || initial > taxa || options.use_clock))
    status = CRPX_ERROR_INVALID_ARGUMENT;
  if (status == CRPX_OK && initial) status = psmc_data_prefix (data, initial, &subset);
  if (status == CRPX_OK) status = psmc_anneal (subset ? subset : data, model, &options, uniform, &stream, &result);
  crpx_character_data_free (subset);
  subset = NULL;
  for (size_t count = initial + 1; initial && count <= taxa && status == CRPX_OK; count++) {
    psmc_result *updated = NULL;
    status = psmc_data_prefix (data, count, &subset);
    if (status == CRPX_OK) status = psmc_online (result, subset, uniform, &stream, &updated);
    crpx_character_data_free (subset);
    subset = NULL;
    if (status == CRPX_OK) {
      psmc_result_free (result);
      result = updated;
    }
  }
  if (status == CRPX_OK) {
    fprintf (stderr, "log evidence %.12g; ESS %.6g/%zu; stages %zu; accepted %zu/%zu\n", result->log_evidence,
             result->ess, result->count, result->stages, result->accepted, result->proposals);
    puts ("weight\tnewick_substitutions");
    for (size_t i = 0; i < result->count && status == CRPX_OK; i++) {
      crpx_tree *tree = NULL;
      char *newick = NULL;
      status = psmc_substitution_tree (result, i, &tree);
      if (status == CRPX_OK) status = crpx_tree_to_newick (tree, &newick);
      if (status == CRPX_OK) printf ("%.17g\t%s\n", result->weights[i], newick);
      crpx_tree_string_free (newick);
      crpx_tree_free (tree);
    }
  }
  if (status != CRPX_OK) fprintf (stderr, "phyloSMC: %s\n", crpx_status_string (status));
  psmc_result_free (result);
  crpx_substitution_model_free (model);
  crpx_character_data_free (data);
  crpx_alignment_free (alignment);
  return status == CRPX_OK ? 0 : 1;
}
static double
uniform (void *state)
{
  return crpx_rng_stream_double (state);
}
static int
usage (void)
{
  fputs ("Usage: phyloSMC ALIGNMENT [--particles N] [--moves N] [--seed N]\n"
         "       [--clock none|strict|lognormal|exponential|autocorrelated] [--initial-taxa N]\n"
         "DNA FASTA/NEXUS input; JC69 model. --initial-taxa enables non-clock online updates.\n"
         "Clock runs use contemporary tips, mean rate 1, coalescent scale 1; relaxed log SD 0.3.\n"
         "Use the C API for other model parameters. Output is a weighted tree sample, not a consensus.\n",
         stderr);
  return 2;
}
static bool
number (const char *text, size_t *value)
{
  if (!text || !*text || *text == '-') return false;
  char *end;
  errno = 0;
  unsigned long long parsed = strtoull (text, &end, 10);
  if (*end || errno || parsed > SIZE_MAX) return false;
  *value = (size_t) parsed;
  return true;
}
