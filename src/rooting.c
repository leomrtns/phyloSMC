/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
#include "smc_internal.h"
#include <rooting.h>
#include <upgma.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

psmc_root_mode
psmc_private_root_mode (const psmc_options *options)
{
  return options->root_mode == PSMC_ROOT_DEFAULT ?
    (options->use_clock ? PSMC_ROOT_MOVABLE : PSMC_UNROOTED) : options->root_mode;
}

crpx_node_id
psmc_private_anchor (const psmc_evaluation *context, const crpx_tree *tree)
{
  const char *name = crpx_character_data_name (context->data, 0);
  for (crpx_node_id v = 0; v < crpx_tree_node_count (tree); v++)
    if (crpx_tree_first_child (tree, v) == CRPX_NODE_NONE && crpx_tree_label (tree, v) &&
        !strcmp (name, crpx_tree_label (tree, v))) return v;
  return CRPX_NODE_NONE;
}

crpx_status
psmc_root_split (const crpx_tree *tree, const crpx_character_data *data, bool *side)
{ /* Resolve tip correspondence explicitly; an arbitrary internal storage orientation is not a rooted split. */
  if (!tree || !data || !side) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_tree_validate_policy (tree, CRPX_TREE_REQUIRE_UNIQUE_TIPS);
  if (status != CRPX_OK) return status;
  size_t n = crpx_character_data_taxa (data);
  if (n != crpx_tree_leaf_count (tree)) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_node_id root = crpx_tree_root (tree), first = crpx_tree_first_child (tree, root);
  crpx_node_id second = crpx_tree_next_sibling (tree, first);
  if (first == CRPX_NODE_NONE || second == CRPX_NODE_NONE || crpx_tree_next_sibling (tree, second) != CRPX_NODE_NONE)
    return CRPX_ERROR_INVALID_ARGUMENT;
  bool *temporary = malloc (n * sizeof (*temporary));
  if (!temporary) return CRPX_ERROR_ALLOCATION;
  for (size_t i = 0; i < n; i++) {
    crpx_node_id tip = CRPX_NODE_NONE;
    for (crpx_node_id v = 0; v < crpx_tree_node_count (tree); v++)
      if (crpx_tree_first_child (tree, v) == CRPX_NODE_NONE &&
          !strcmp (crpx_tree_label (tree, v), crpx_character_data_name (data, i))) { tip = v; break; }
    if (tip == CRPX_NODE_NONE) { free (temporary); return CRPX_ERROR_INVALID_ARGUMENT; }
    while (crpx_tree_parent (tree, tip) != root) tip = crpx_tree_parent (tree, tip);
    temporary[i] = tip == first;
  }
  memcpy (side, temporary, n * sizeof (*side));
  free (temporary);
  return CRPX_OK;
}

crpx_status
psmc_estimate_root_split (const crpx_character_data *data, bool *side)
{ /* Build an explicitly data-dependent constraint: mismatch-distance NJ followed by midpoint rooting. */
  if (!data || !side) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t n = crpx_character_data_taxa (data), states = crpx_character_data_states (data);
  if (n < 2) return CRPX_ERROR_INVALID_ARGUMENT;
  const double *observations = crpx_character_data_pattern_data (data);
  const size_t *weights = crpx_character_data_pattern_weights (data);
  const char **names = malloc (n * sizeof (*names));
  crpx_distance_matrix_t distances = NULL;
  crpx_tree *tree = NULL;
  if (!names) return CRPX_ERROR_ALLOCATION;
  crpx_status status = crpx_distance_matrix_new (n, &distances);
  for (size_t i = 0; i < n; i++) names[i] = crpx_character_data_name (data, i);
  for (size_t i = 0; i < n && status == CRPX_OK; i++) for (size_t j = 0; j < i && status == CRPX_OK; j++) {
    double mismatch = 0, total = 0;
    for (size_t p = 0; p < crpx_character_data_patterns (data); p++) {
      size_t a = SIZE_MAX, b = SIZE_MAX, ca = 0, cb = 0;
      for (size_t k = 0; k < states; k++) {
        if (observations[(p * n + i) * states + k] > 0) { a = k; ca++; }
        if (observations[(p * n + j) * states + k] > 0) { b = k; cb++; }
      }
      if (ca == 1 && cb == 1) { total += weights[p]; mismatch += weights[p] * (a != b); }
    }
    if (!total) status = CRPX_ERROR_INVALID_ARGUMENT;
    else status = crpx_distance_matrix_set (distances, i, j, mismatch / total);
  }
  if (status == CRPX_OK) status = crpx_tree_neighbour_joining (distances, names, false, &tree);
  crpx_root_location location;
  crpx_node_id root;
  if (status == CRPX_OK) status = crpx_tree_root_midpoint (tree, NULL, &location);
  if (status == CRPX_OK) status = crpx_tree_reroot_edge (tree, location.child, location.fraction,
                                                        CRPX_TREE_KEEP_OLD_ROOT, &root);
  if (status == CRPX_OK) status = psmc_root_split (tree, data, side);
  crpx_tree_free (tree); crpx_distance_matrix_free (distances); free (names);
  return status;
}
