/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
#include "phylosmc.h"
#include <global/random_number.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
static double uniform (void *state);
static int usage (void);
static bool number (const char *text, size_t *value);
static bool real (const char *text, double *value);
static bool vector (const char *text, double *values, size_t n);
static crpx_status read_tree (const char *path, crpx_tree **out);

int
main (int argc, char **argv)
{ /* A DNA front end with explicit statistical settings; general character models remain available through C. */
  if (argc < 2) return usage ();
  if (!strcmp (argv[1], "--help")) { (void)usage (); return 0; }
  psmc_options options;
  psmc_options_init (&options);
  size_t seed = 1, initial = 0;
  bool joint = true, infer_mean = false;
  psmc_model_kind kind = PSMC_MODEL_HKY;
  const char *tree_file = NULL, *outgroup = NULL;
  double pi[4] = {.25,.25,.25,.25}, exchange[6] = {1,1,1,1,1,1}, kappa = 2;
  for (int i = 2; i < argc; i++) {
    if (i + 1 == argc) return usage ();
    const char *name = argv[i++], *value = argv[i];
    if (!strcmp (name, "--clock")) {
      options.use_clock = strcmp (value, "none") != 0;
      if (!strcmp (value, "strict")) options.clock.kind = CRPX_CLOCK_STRICT;
      else if (!strcmp (value, "lognormal")) options.clock.kind = CRPX_CLOCK_LOGNORMAL;
      else if (!strcmp (value, "exponential")) options.clock.kind = CRPX_CLOCK_EXPONENTIAL;
      else if (!strcmp (value, "autocorrelated")) options.clock.kind = CRPX_CLOCK_AUTOCORRELATED;
      else if (strcmp (value, "none")) return usage ();
    } else if (!strcmp (name, "--root-mode")) {
      if (!strcmp (value, "unrooted")) options.root_mode = PSMC_UNROOTED;
      else if (!strcmp (value, "fixed")) options.root_mode = PSMC_ROOT_FIXED;
      else if (!strcmp (value, "movable")) options.root_mode = PSMC_ROOT_MOVABLE;
      else return usage ();
    } else if (!strcmp (name, "--model")) {
      if (!strcmp (value, "jc69")) kind = PSMC_MODEL_FIXED;
      else if (!strcmp (value, "f81")) kind = PSMC_MODEL_F81;
      else if (!strcmp (value, "hky")) kind = PSMC_MODEL_HKY;
      else if (!strcmp (value, "gtr")) kind = PSMC_MODEL_GTR;
      else return usage ();
    } else if (!strcmp (name, "--parameters")) {
      if (strcmp (value, "joint") && strcmp (value, "fixed")) return usage ();
      joint = !strcmp (value, "joint");
    } else if (!strcmp (name, "--infer-clock-mean")) {
      if (strcmp (value, "yes") && strcmp (value, "no")) return usage ();
      infer_mean = !strcmp (value, "yes");
    } else if (!strcmp (name, "--tree-prior")) {
      if (!strcmp (value, "coalescent")) options.tree_prior.kind = CRPX_TREE_PRIOR_COALESCENT;
      else if (!strcmp (value, "yule")) options.tree_prior.kind = CRPX_TREE_PRIOR_YULE;
      else if (!strcmp (value, "birth-death")) options.tree_prior.kind = CRPX_TREE_PRIOR_BIRTH_DEATH;
      else return usage ();
    } else if (!strcmp (name, "--tree")) tree_file = value;
    else if (!strcmp (name, "--outgroup")) outgroup = value;
    else if (!strcmp (name, "--frequencies")) { if (!vector (value, pi, 4)) return usage (); }
    else if (!strcmp (name, "--exchangeabilities")) { if (!vector (value, exchange, 6)) return usage (); }
    else if (!strcmp (name, "--prior-medians")) { if (!vector (value, options.priors.median, 6)) return usage (); }
    else if (!strcmp (name, "--prior-log-sds")) { if (!vector (value, options.priors.log_sd, 6)) return usage (); }
    else {
      size_t parsed;
      double x;
      if (!strcmp (name, "--particles") || !strcmp (name, "--moves") || !strcmp (name, "--seed") ||
          !strcmp (name, "--initial-taxa") || !strcmp (name, "--workers")) {
        if (!number (value, &parsed)) return usage ();
        if (!strcmp (name, "--particles")) options.particles = parsed;
        else if (!strcmp (name, "--moves")) options.mcmc_steps = parsed;
        else if (!strcmp (name, "--seed")) seed = parsed;
        else if (!strcmp (name, "--workers")) options.workers = parsed;
        else initial = parsed;
      } else {
        if (!real (value, &x)) return usage ();
        if (!strcmp (name, "--branch-rate")) options.branch_prior_rate = x;
        else if (!strcmp (name, "--clock-mean")) options.clock.mean_rate = x;
        else if (!strcmp (name, "--clock-dispersion")) options.clock.dispersion = x;
        else if (!strcmp (name, "--coalescent-scale")) options.tree_prior.coalescent_scale = x;
        else if (!strcmp (name, "--birth-rate")) options.tree_prior.birth_rate = x;
        else if (!strcmp (name, "--death-rate")) options.tree_prior.death_rate = x;
        else if (!strcmp (name, "--origin-age")) options.tree_prior.origin_age = x;
        else if (!strcmp (name, "--kappa")) kappa = x;
        else if (!strcmp (name, "--frequency-concentration")) options.priors.frequency_concentration = x;
        else if (!strcmp (name, "--exchangeability-concentration")) options.priors.exchangeability_concentration = x;
        else if (!strcmp (name, "--kappa-prior-median")) options.priors.kappa_median = x;
        else if (!strcmp (name, "--kappa-prior-log-sd")) options.priors.kappa_log_sd = x;
        else return usage ();
      }
    }
  }
  if ((tree_file && outgroup) || ((tree_file || outgroup) && options.root_mode != PSMC_ROOT_FIXED) ||
      (infer_mean && (!joint || !options.use_clock))) return usage ();
  options.substitution = joint ? kind : PSMC_MODEL_FIXED;
  if (joint) {
    if (!options.use_clock) options.infer = PSMC_INFER_BRANCH_RATE;
    else {
      options.infer = options.tree_prior.kind == CRPX_TREE_PRIOR_COALESCENT ?
                      PSMC_INFER_COALESCENT_SCALE : PSMC_INFER_BIRTH_RATE;
      if (options.tree_prior.kind == CRPX_TREE_PRIOR_BIRTH_DEATH) options.infer |= PSMC_INFER_DEATH_RATE;
      if (options.clock.kind == CRPX_CLOCK_LOGNORMAL || options.clock.kind == CRPX_CLOCK_AUTOCORRELATED)
        options.infer |= PSMC_INFER_CLOCK_DISPERSION;
      if (infer_mean) options.infer |= PSMC_INFER_CLOCK_MEAN;
    }
  }
  crpx_alignment *alignment = NULL;
  crpx_character_data *data = NULL, *subset = NULL;
  crpx_substitution_model *model = NULL;
  psmc_result *result = NULL;
  bool *side = NULL;
  crpx_rng_stream stream;
  crpx_rng_stream_init (&stream, seed, 0);
  crpx_status status = crpx_alignment_read (argv[1], true, &alignment);
  if (status == CRPX_OK) status = crpx_character_data_from_alignment_dna (alignment, true, &data);
  if (status == CRPX_OK) {
    if (kind == PSMC_MODEL_FIXED) status = crpx_substitution_model_poisson (4, &model);
    else if (kind == PSMC_MODEL_F81) status = crpx_substitution_model_f81 (4, pi, &model);
    else if (kind == PSMC_MODEL_HKY) status = crpx_substitution_model_hky (pi, kappa, &model);
    else status = crpx_substitution_model_gtr (pi, exchange, &model);
  }
  size_t taxa = crpx_character_data_taxa (data);
  if (status == CRPX_OK && options.root_mode == PSMC_ROOT_FIXED) {
    side = calloc (taxa, sizeof (*side));
    if (!side) status = CRPX_ERROR_ALLOCATION;
    else if (tree_file) {
      crpx_tree *tree = NULL;
      status = read_tree (tree_file, &tree);
      if (status == CRPX_OK) status = psmc_root_split (tree, data, side);
      crpx_tree_free (tree);
    } else if (outgroup) {
      char *names = malloc (strlen (outgroup) + 1);
      if (!names) status = CRPX_ERROR_ALLOCATION;
      else {
        strcpy (names, outgroup);
        for (char *name = strtok (names, ","); name && status == CRPX_OK; name = strtok (NULL, ",")) {
          size_t i = 0;
          while (i < taxa && strcmp (name, crpx_character_data_name (data, i))) i++;
          if (i == taxa || side[i]) status = CRPX_ERROR_INVALID_ARGUMENT;
          else side[i] = true;
        }
        free (names);
      }
    } else status = psmc_estimate_root_split (data, side);
    options.root_side = side;
    if (status == CRPX_OK) {
      fputs ("Conditioning on root split:", stderr);
      for (size_t i = 0; i < taxa; i++) if (side[i]) fprintf (stderr, " %s", crpx_character_data_name (data, i));
      fputs (" |", stderr);
      for (size_t i = 0; i < taxa; i++) if (!side[i]) fprintf (stderr, " %s", crpx_character_data_name (data, i));
      fputc ('\n', stderr);
    }
  }
  if (status == CRPX_OK && initial && (initial < 2 || initial > taxa || options.use_clock)) status = CRPX_ERROR_INVALID_ARGUMENT;
  if (status == CRPX_OK && initial) status = psmc_data_prefix (data, initial, &subset);
  if (status == CRPX_OK) status = psmc_anneal (subset ? subset : data, model, &options, uniform, &stream, &result);
  crpx_character_data_free (subset); subset = NULL;
  for (size_t count = initial + 1; initial && count <= taxa && status == CRPX_OK; count++) {
    psmc_result *updated = NULL;
    status = psmc_data_prefix (data, count, &subset);
    if (status == CRPX_OK) status = psmc_online_split (result, subset, side, uniform, &stream, &updated);
    crpx_character_data_free (subset); subset = NULL;
    if (status == CRPX_OK) { psmc_result_free (result); result = updated; }
  }
  if (status == CRPX_OK) {
    fprintf (stderr, "log evidence %.12g; ESS %.6g/%zu; stages %zu; accepted %zu/%zu; initial ancestors %zu\n",
             result->log_evidence, result->ess, result->count, result->stages, result->accepted,
             result->proposals, result->surviving_ancestors);
    puts ("weight\tbranch_rate\tclock_mean\tclock_dispersion\tcoalescent_scale\tbirth_rate\tdeath_rate\tpi_A\tpi_C\tpi_G\tpi_T\tkappa\tr_AC\tr_AG\tr_AT\tr_CG\tr_CT\tr_GT\tnewick_substitutions");
    for (size_t i = 0; i < result->count && status == CRPX_OK; i++) {
      crpx_tree *tree = NULL;
      char *newick = NULL;
      status = psmc_substitution_tree (result, i, &tree);
      if (status == CRPX_OK) status = crpx_tree_to_newick (tree, &newick);
      if (status == CRPX_OK) {
        const psmc_parameters *p = &result->particles[i].parameters;
        printf ("%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g", result->weights[i], p->branch_rate,
                p->clock_mean, p->clock_dispersion, p->coalescent_scale, p->birth_rate, p->death_rate);
        for (size_t j = 0; j < 4; j++) printf ("\t%.17g", p->frequencies[j]);
        printf ("\t%.17g", p->kappa);
        for (size_t j = 0; j < 6; j++) printf ("\t%.17g", p->exchangeabilities[j]);
        bool unrooted = options.root_mode == PSMC_UNROOTED ||
                        (options.root_mode == PSMC_ROOT_DEFAULT && !options.use_clock);
        printf ("\t%s%s\n", unrooted ? "[&U]" : "[&R]", newick);
      }
      crpx_tree_string_free (newick); crpx_tree_free (tree);
    }
  }
  if (status != CRPX_OK) fprintf (stderr, "phyloSMC: %s\n", crpx_status_string (status));
  psmc_result_free (result); crpx_substitution_model_free (model);
  crpx_character_data_free (data); crpx_alignment_free (alignment); free (side);
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
  fputs ("Usage: phyloSMC ALIGNMENT [--particles N] [--moves N] [--workers N] [--seed N]\n"
    "  [--model jc69|f81|hky|gtr] [--parameters joint|fixed] [--initial-taxa N]\n"
    "  [--root-mode unrooted|fixed|movable] [--outgroup NAME,NAME | --tree FILE]\n"
    "  [--clock none|strict|lognormal|exponential|autocorrelated]\n"
    "  [--tree-prior coalescent|yule|birth-death] [--origin-age X]\n"
    "  [--branch-rate X] [--clock-mean X] [--clock-dispersion X] [--coalescent-scale X]\n"
    "  [--birth-rate X] [--death-rate X] [--infer-clock-mean yes|no]\n"
    "  [--frequencies A,C,G,T] [--kappa X] [--exchangeabilities AC,AG,AT,CG,CT,GT]\n"
    "  [--prior-medians BRANCH,MEAN,DISPERSION,COALESCENT,BIRTH,DEATH]\n"
    "  [--prior-log-sds BRANCH,MEAN,DISPERSION,COALESCENT,BIRTH,DEATH]\n"
    "  [--frequency-concentration X] [--exchangeability-concentration X]\n"
    "  [--kappa-prior-median X] [--kappa-prior-log-sd X]\n"
    "Defaults: unrooted HKY with joint inference; clocks use contemporary-tip coalescent trees.\n"
    "Fixed rooting without --tree/--outgroup uses midpoint-rooted NJ mismatch distances.\n"
    "Workers are processes on one POSIX machine; 1 is serial. See README.md and docs/tutorial.md.\n", stderr);
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
  *value = (size_t)parsed;
  return true;
}
static bool
real (const char *text, double *value)
{
  char *end;
  errno = 0;
  double x = strtod (text, &end);
  if (end == text || *end || errno || !isfinite (x)) return false;
  *value = x;
  return true;
}
static bool
vector (const char *text, double *values, size_t n)
{
  for (size_t i = 0; i < n; i++) {
    char *end;
    errno = 0;
    double x = strtod (text, &end);
    if (end == text || errno || !isfinite (x) || (i + 1 == n ? *end != 0 : *end != ',')) return false;
    values[i] = x; text = end + 1;
  }
  return true;
}
static crpx_status
read_tree (const char *path, crpx_tree **out)
{
  FILE *file = fopen (path, "rb");
  if (!file) return CRPX_ERROR_IO;
  if (fseek (file, 0, SEEK_END)) { fclose (file); return CRPX_ERROR_IO; }
  long bytes = ftell (file);
  if (bytes < 0 || (unsigned long)bytes >= SIZE_MAX) { fclose (file); return CRPX_ERROR_IO; }
  rewind (file);
  char *text = malloc ((size_t)bytes + 1);
  if (!text) { fclose (file); return CRPX_ERROR_ALLOCATION; }
  crpx_status status = CRPX_ERROR_IO;
  if (fread (text, 1, bytes, file) == (size_t)bytes) {
    text[bytes] = 0;
    status = crpx_tree_from_newick (text, out, NULL);
  }
  fclose (file); free (text);
  return status;
}
