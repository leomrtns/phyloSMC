/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file smc_internal.h
    \brief State and functions shared by the SMC stages, particle calculations and proposals. */
#ifndef PHYLOSMC_INTERNAL_H
#define PHYLOSMC_INTERNAL_H
#include "phylosmc.h"
typedef struct {
  const crpx_character_data *data;
  const crpx_substitution_model *model;
  const psmc_options *options;
  crpx_uniform_function uniform;
  void *state;
  crpx_likelihood_workspace *likelihood;
  crpx_tree_workspace *traversal;
  double *lengths;
  size_t nodes;
} psmc_evaluation;
crpx_status psmc_private_prepare (psmc_evaluation *context, const crpx_character_data *data,
                                  const crpx_substitution_model *model, const psmc_options *options,
                                  crpx_uniform_function uniform, void *state);
void psmc_private_clear (psmc_evaluation *context);
crpx_status psmc_private_result_new (psmc_evaluation *context, psmc_result **out);
crpx_status psmc_private_draw_prior (psmc_evaluation *context, psmc_particle *particle);
crpx_status psmc_private_evaluate (psmc_evaluation *context, psmc_particle *particle);
crpx_status psmc_private_copy_particle (const psmc_particle *source, psmc_particle *copy);
void psmc_private_free_particle (psmc_particle *particle);
crpx_status psmc_private_draw (psmc_evaluation *context, double *value);
crpx_status psmc_private_index_draw (psmc_evaluation *context, size_t count, size_t *index);
crpx_status psmc_private_resample (psmc_evaluation *context, psmc_result *result);
crpx_status psmc_private_rejuvenate (psmc_evaluation *context, psmc_result *result);
crpx_status psmc_private_validate_online (const psmc_result *previous, const crpx_character_data *data);
#endif
