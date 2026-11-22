/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
#include "smc_internal.h"
#include <prob_distribution.h>
#include <math.h>
#include <string.h>

static crpx_status simplex (psmc_evaluation *context, double *values, size_t n, double alpha, bool sample,
                             double *density);

void
psmc_private_models (const psmc_options *options, const psmc_parameters *p,
                      crpx_tree_prior *prior, crpx_clock_model *clock)
{
  *prior = options->tree_prior;
  *clock = options->clock;
  prior->coalescent_scale = p->coalescent_scale;
  prior->birth_rate = p->birth_rate;
  prior->death_rate = p->death_rate;
  clock->mean_rate = p->clock_mean;
  clock->dispersion = p->clock_dispersion;
}

crpx_status
psmc_private_parameters (psmc_evaluation *context, psmc_particle *particle, bool sample, double *log_prior)
{ /* Proper hyperpriors are part of every annealed target, including temperature zero. */
  const psmc_options *o = context->options;
  psmc_parameters *p = &particle->parameters;
  double density = 0;
  if (sample) {
    *p = (psmc_parameters){.branch_rate=o->branch_prior_rate, .clock_mean=o->clock.mean_rate,
      .clock_dispersion=o->clock.dispersion, .coalescent_scale=o->tree_prior.coalescent_scale,
      .birth_rate=o->tree_prior.birth_rate, .death_rate=o->tree_prior.death_rate,
      .frequencies={.25,.25,.25,.25}, .kappa=2,
      .exchangeabilities={1./6,1./6,1./6,1./6,1./6,1./6}};
  }
  if (sample && o->substitution == PSMC_MODEL_FIXED) {
    const double *pi = crpx_substitution_model_frequencies (context->model);
    for (size_t i = 0; i < 4; i++) p->frequencies[i] = crpx_substitution_model_states (context->model) == 4 ? pi[i] : NAN;
  }
  if (sample && o->substitution != PSMC_MODEL_HKY) p->kappa = NAN;
  if (sample && o->substitution != PSMC_MODEL_GTR)
    for (size_t i = 0; i < 6; i++) p->exchangeabilities[i] = NAN;
  double *positive[] = {&p->branch_rate,&p->clock_mean,&p->clock_dispersion,&p->coalescent_scale,&p->birth_rate,&p->death_rate};
  for (unsigned i = 0; i < 6; i++) if (o->infer & (1u << i)) {
    if (sample) {
      double u;
      crpx_status status = psmc_private_draw (context, &u);
      if (status != CRPX_OK) return status;
      *positive[i] = exp (crpx_qnorm (u, log (o->priors.median[i]), o->priors.log_sd[i], false));
    }
    double value = crpx_dlnorm (*positive[i], log (o->priors.median[i]), o->priors.log_sd[i], true);
    if (!isfinite (value)) return CRPX_ERROR_OVERFLOW;
    density += value;
  }
  if (o->substitution != PSMC_MODEL_FIXED) {
    crpx_status status = simplex (context, p->frequencies, 4, o->priors.frequency_concentration, sample, &density);
    if (status != CRPX_OK) return status;
    if (o->substitution == PSMC_MODEL_HKY) {
      if (sample) {
        double u;
        status = psmc_private_draw (context, &u);
        if (status != CRPX_OK) return status;
        p->kappa = exp (crpx_qnorm (u, log (o->priors.kappa_median), o->priors.kappa_log_sd, false));
      }
      density += crpx_dlnorm (p->kappa, log (o->priors.kappa_median), o->priors.kappa_log_sd, true);
    }
    if (o->substitution == PSMC_MODEL_GTR) {
      status = simplex (context, p->exchangeabilities, 6, o->priors.exchangeability_concentration, sample, &density);
      if (status != CRPX_OK) return status;
    }
  }
  if (!isfinite (density)) return CRPX_ERROR_OVERFLOW;
  *log_prior = density;
  return CRPX_OK;
}

crpx_status
psmc_private_parameter_proposal (psmc_evaluation *context, psmc_particle *particle,
                                 double *reverse_forward, bool *valid)
{ /* Multipliers for positive scalars; symmetric pair transfers within each simplex. */
  const psmc_options *o = context->options;
  psmc_parameters *p = &particle->parameters;
  double *positive[] = {&p->branch_rate,&p->clock_mean,&p->clock_dispersion,&p->coalescent_scale,&p->birth_rate,&p->death_rate};
  size_t count = o->substitution == PSMC_MODEL_FIXED ? 0 : o->substitution == PSMC_MODEL_F81 ? 1 : 2;
  for (unsigned i = 0; i < 6; i++) if (o->infer & (1u << i)) count++;
  if (!count) { *valid = false; return CRPX_OK; }
  size_t selected;
  crpx_status status = psmc_private_index_draw (context, count, &selected);
  if (status != CRPX_OK) return status;
  double *scalar = NULL, *vector = NULL;
  size_t n = 0;
  for (unsigned i = 0; i < 6; i++) if (o->infer & (1u << i)) {
    if (!selected) { scalar = positive[i]; break; }
    selected--;
  }
  if (!scalar) {
    if (!selected) { vector = p->frequencies; n = 4; }
    else if (o->substitution == PSMC_MODEL_HKY) scalar = &p->kappa;
    else { vector = p->exchangeabilities; n = 6; }
  }
  double u;
  status = psmc_private_draw (context, &u);
  if (status != CRPX_OK) return status;
  double change = (2 * u - 1) * o->proposal_width;
  if (scalar) {
    *scalar *= exp (change);
    *reverse_forward = change;
    *valid = isfinite (*scalar) && *scalar > 0;
  } else {
    size_t a, b;
    status = psmc_private_index_draw (context, n, &a);
    if (status == CRPX_OK) status = psmc_private_index_draw (context, n - 1, &b);
    if (status != CRPX_OK) return status;
    if (b >= a) b++;
    double total = vector[a] + vector[b];
    vector[a] += change * total;
    vector[b] = total - vector[a];
    *valid = vector[a] > 0 && vector[b] > 0;
  }
  return status;
}

static crpx_status
simplex (psmc_evaluation *context, double *values, size_t n, double alpha, bool sample, double *density)
{
  double total = 0;
  if (sample) {
    for (size_t i = 0; i < n; i++) {
      double u;
      crpx_status status = psmc_private_draw (context, &u);
      if (status != CRPX_OK) return status;
      values[i] = crpx_qgamma (u, alpha, 1, false);
      if (!isfinite (values[i]) || values[i] <= 0) return CRPX_ERROR_OVERFLOW;
      total += values[i];
    }
    for (size_t i = 0; i < n; i++) values[i] /= total;
  }
  *density += crpx_lgammafn (n * alpha, NULL) - n * crpx_lgammafn (alpha, NULL);
  for (size_t i = 0; i < n; i++) {
    if (!isfinite (values[i]) || values[i] <= 0) return CRPX_ERROR_INVALID_ARGUMENT;
    *density += (alpha - 1) * log (values[i]);
  }
  return CRPX_OK;
}
