/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
#include "smc_internal.h"
#include <particle_weights.h>
#include <topology_randomise.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
static bool unit_weights (const crpx_character_data *data);

crpx_status
psmc_data_prefix (const crpx_character_data *data, size_t taxa, crpx_character_data **out)
{ /* Re-expand compact patterns before selecting taxa, then compact the new observations. */
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!data || !taxa || taxa > crpx_character_data_taxa (data) || !unit_weights (data))
    return CRPX_ERROR_INVALID_ARGUMENT;
  size_t sites = crpx_character_data_sites (data), states = crpx_character_data_states (data);
  if (taxa > SIZE_MAX / states / sizeof (double) || (sites && taxa * states > SIZE_MAX / sites / sizeof (double)))
    return CRPX_ERROR_OVERFLOW;
  double *values = sites ? malloc (sites * taxa * states * sizeof (*values)) : NULL;
  const char **names = malloc (taxa * sizeof (*names));
  if ((sites && !values) || !names) {
    free (values);
    free (names);
    return CRPX_ERROR_ALLOCATION;
  }
  const size_t *patterns = crpx_character_data_site_patterns (data);
  const double *observed = crpx_character_data_pattern_data (data);
  for (size_t i = 0; i < taxa; i++)
    names[i] = crpx_character_data_name (data, i);
  for (size_t site = 0; site < sites; site++)
    memcpy (values + site * taxa * states, observed + patterns[site] * crpx_character_data_taxa (data) * states,
            taxa * states * sizeof (double));
  crpx_status status = crpx_character_data_from_tip_likelihoods (names, taxa, sites, states, values, NULL, true, out);
  free (values);
  free (names);
  return status;
}

static bool
unit_weights (const crpx_character_data *data)
{
  size_t patterns = crpx_character_data_patterns (data), sites = crpx_character_data_sites (data);
  size_t *counts = patterns ? calloc (patterns, sizeof (*counts)) : NULL;
  if (patterns && !counts) return false;
  const size_t *mapping = crpx_character_data_site_patterns (data),
               *weights = crpx_character_data_pattern_weights (data);
  for (size_t site = 0; site < sites; site++)
    counts[mapping[site]]++;
  bool valid = true;
  for (size_t i = 0; i < patterns; i++)
    if (counts[i] != weights[i]) valid = false;
  free (counts);
  return valid;
}

crpx_status
psmc_private_validate_online (const psmc_result *previous, const crpx_character_data *data)
{
  if (!previous || !data || previous->beta != 1 || !previous->count) return CRPX_ERROR_INVALID_ARGUMENT;
  if (previous->options.use_clock) return CRPX_ERROR_UNSUPPORTED;
  const crpx_character_data *old = previous->data;
  size_t taxa = crpx_character_data_taxa (old), sites = crpx_character_data_sites (old),
         states = crpx_character_data_states (old);
  if (crpx_character_data_taxa (data) != taxa + 1 || crpx_character_data_sites (data) != sites ||
      crpx_character_data_states (data) != states || !unit_weights (old) || !unit_weights (data))
    return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < taxa; i++)
    if (strcmp (crpx_character_data_name (old, i), crpx_character_data_name (data, i)))
      return CRPX_ERROR_INVALID_ARGUMENT;
  const size_t *a = crpx_character_data_site_patterns (old), *b = crpx_character_data_site_patterns (data);
  const double *x = crpx_character_data_pattern_data (old), *y = crpx_character_data_pattern_data (data);
  for (size_t site = 0; site < sites; site++)
    for (size_t taxon = 0; taxon < taxa; taxon++)
      for (size_t state = 0; state < states; state++)
        if (x[(a[site] * taxa + taxon) * states + state] != y[(b[site] * (taxa + 1) + taxon) * states + state])
          return CRPX_ERROR_INVALID_ARGUMENT;
  return CRPX_OK;
}
