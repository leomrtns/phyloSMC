/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
#define _POSIX_C_SOURCE 200809L
#include "smc_internal.h"
#include <global/random_number.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifndef _WIN32
#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>
#include <signal.h>
#include <errno.h>
#endif

static double
worker_uniform (void *state)
{
  return crpx_rng_stream_double (state);
}


static crpx_status
canonicalise (psmc_particle *particle)
{ /* Match the process transport representation also in serial runs, keeping proposal node ordering reproducible. */
  char *newick = NULL;
  crpx_tree *tree = NULL;
  crpx_tree_workspace *order = NULL;
  double *ages = NULL, *rates = NULL;
  size_t n = crpx_tree_node_count (particle->tree);
  crpx_status status = crpx_tree_to_newick (particle->tree, &newick);
  if (status == CRPX_OK) status = crpx_tree_from_newick (newick, &tree, NULL);
  if (status != CRPX_OK) goto done;
  if (particle->ages) {
    ages = malloc (n * sizeof (*ages)); rates = malloc (n * sizeof (*rates));
    if (!ages || !rates) { status = CRPX_ERROR_ALLOCATION; goto done; }
    status = crpx_tree_workspace_new (CRPX_TREE_WORKSPACE_TRAVERSAL, &order);
    if (status == CRPX_OK) status = crpx_tree_workspace_refresh (order, particle->tree);
    if (status != CRPX_OK) goto done;
    const crpx_node_id *nodes = crpx_tree_workspace_preorder (order);
    for (size_t i = 0; i < n; i++) { ages[i] = particle->ages[nodes[i]]; rates[i] = particle->rates[nodes[i]]; }
    status = crpx_tree_workspace_refresh (order, tree);
    if (status != CRPX_OK) goto done;
    nodes = crpx_tree_workspace_preorder (order);
    for (size_t i = 0; i < n; i++) { particle->ages[nodes[i]] = ages[i]; particle->rates[nodes[i]] = rates[i]; }
  }
  crpx_tree_free (particle->tree);
  particle->tree = tree; tree = NULL;
done:
  crpx_tree_workspace_free (order); crpx_tree_free (tree); crpx_tree_string_free (newick);
  free (ages); free (rates);
  return status;
}

#ifndef _WIN32
static bool
transfer (FILE *file, void *data, size_t size, bool writing)
{
  return writing ? fwrite (data, 1, size, file) == size : fread (data, 1, size, file) == size;
}

static crpx_status
particle_transfer (FILE *file, psmc_particle *particle, bool writing)
{ /* Newick preserves ordered topology; attach clock arrays by preorder rather than assuming parser node IDs. */
  char *newick = NULL;
  size_t bytes = 0;
  crpx_status status = CRPX_OK;
  crpx_tree_workspace *order = NULL;
  if (writing) {
    status = crpx_tree_to_newick (particle->tree, &newick);
    if (status != CRPX_OK) return status;
    bytes = strlen (newick) + 1;
  }
  if (!transfer (file, &bytes, sizeof (bytes), writing) || !bytes) { status = CRPX_ERROR_IO; goto done; }
  if (!writing) {
    newick = malloc (bytes);
    if (!newick) { status = CRPX_ERROR_ALLOCATION; goto done; }
  }
  if (!transfer (file, newick, bytes, writing)) { status = CRPX_ERROR_IO; goto done; }
  if (!writing) {
    if (newick[bytes - 1]) { status = CRPX_ERROR_IO; goto done; }
    status = crpx_tree_from_newick (newick, &particle->tree, NULL);
    if (status != CRPX_OK) goto done;
  }
  if (!transfer (file, &particle->log_likelihood, sizeof (double), writing) ||
      !transfer (file, &particle->log_prior, sizeof (double), writing) ||
      !transfer (file, &particle->parameters, sizeof (particle->parameters), writing)) {
    status = CRPX_ERROR_IO; goto done;
  }
  bool clock = particle->ages != NULL;
  if (!transfer (file, &clock, sizeof (clock), writing)) { status = CRPX_ERROR_IO; goto done; }
  if (clock) {
    size_t n = crpx_tree_node_count (particle->tree);
    if (!writing) {
      particle->ages = malloc (n * sizeof (double));
      particle->rates = malloc (n * sizeof (double));
      if (!particle->ages || !particle->rates) { status = CRPX_ERROR_ALLOCATION; goto done; }
    }
    status = crpx_tree_workspace_new (CRPX_TREE_WORKSPACE_TRAVERSAL, &order);
    if (status == CRPX_OK) status = crpx_tree_workspace_refresh (order, particle->tree);
    if (status != CRPX_OK) goto done;
    const crpx_node_id *nodes = crpx_tree_workspace_preorder (order);
    for (size_t i = 0; i < n; i++)
      if (!transfer (file, particle->ages + nodes[i], sizeof (double), writing) ||
          !transfer (file, particle->rates + nodes[i], sizeof (double), writing)) {
        status = CRPX_ERROR_IO; break;
      }
  }
done:
  crpx_tree_workspace_free (order);
  crpx_tree_string_free (newick);
  return status;
}
#endif

crpx_status
psmc_private_execute (psmc_evaluation *context, psmc_result *result, bool prior)
{ /* Each particle receives a jumped stream before scheduling. Coordinating draws never occur in a child. */
  size_t workers = result->options.workers;
  if (workers > result->count) workers = result->count;
#ifdef _WIN32
  if (workers > 1) return CRPX_ERROR_UNSUPPORTED;
#endif
  double a, b;
  crpx_status status = psmc_private_draw (context, &a);
  if (status == CRPX_OK) status = psmc_private_draw (context, &b);
  if (status != CRPX_OK) return status;
  uint64_t seed = ((uint64_t)(a * 0x1p32) << 32) | (uint64_t)(b * 0x1p32);
  crpx_rng_stream *streams = malloc (result->count * sizeof (*streams));
  if (!streams) return CRPX_ERROR_ALLOCATION;
  crpx_rng_stream_init (streams, seed, 0);
  for (size_t i = 1; i < result->count; i++) { streams[i] = streams[i - 1]; crpx_rng_stream_jump (streams + i); }
  psmc_evaluation local = *context;
  local.uniform = worker_uniform;
  if (workers == 1) {
    for (size_t i = 0; i < result->count && status == CRPX_OK; i++) {
      local.state = streams + i;
      status = prior ? psmc_private_draw_prior (&local, result->particles + i) :
                       psmc_private_rejuvenate_one (&local, result, i);
      if (status == CRPX_OK) status = canonicalise (result->particles + i);
    }
    free (streams);
    return status;
  }
#ifndef _WIN32
  pid_t *processes = calloc (workers, sizeof (*processes));
  FILE **channels = calloc (workers, sizeof (*channels));
  size_t started = 0;
  if (!processes || !channels) { status = CRPX_ERROR_ALLOCATION; goto done; }
  for (size_t w = 0; w < workers; w++) {
    int pipefd[2];
    if (pipe (pipefd)) { status = CRPX_ERROR_IO; break; }
    pid_t pid = fork ();
    if (pid < 0) { close (pipefd[0]); close (pipefd[1]); status = CRPX_ERROR_IO; break; }
    if (!pid) {
      close (pipefd[0]);
      for (size_t j = 0; j < started; j++) fclose (channels[j]);
      FILE *channel = fdopen (pipefd[1], "wb");
      if (!channel) _exit (1);
      for (size_t i = w; i < result->count; i += workers) {
        local.state = streams + i;
        result->proposals = result->accepted = 0;
        status = prior ? psmc_private_draw_prior (&local, result->particles + i) :
                         psmc_private_rejuvenate_one (&local, result, i);
        if (!transfer (channel, &status, sizeof (status), true) || status != CRPX_OK) break;
        status = particle_transfer (channel, result->particles + i, true);
        if (status != CRPX_OK ||
            !transfer (channel, &result->proposals, sizeof (size_t), true) ||
            !transfer (channel, &result->accepted, sizeof (size_t), true)) { status = CRPX_ERROR_IO; break; }
        psmc_private_free_particle (result->particles + i);
      }
      if (fclose (channel)) status = CRPX_ERROR_IO;
      _exit (status == CRPX_OK ? 0 : 1);
    }
    close (pipefd[1]);
    processes[w] = pid;
    channels[w] = fdopen (pipefd[0], "rb");
    started++;
    if (!channels[w]) { close (pipefd[0]); status = CRPX_ERROR_IO; break; }
  }
  for (size_t w = 0; w < started && status == CRPX_OK; w++)
    for (size_t i = w; i < result->count; i += workers) {
      crpx_status child_status;
      psmc_particle received = {0};
      size_t proposed = 0, accepted = 0;
      if (!transfer (channels[w], &child_status, sizeof (child_status), false)) { status = CRPX_ERROR_IO; break; }
      if (child_status != CRPX_OK) { status = child_status; break; }
      status = particle_transfer (channels[w], &received, false);
      if (status == CRPX_OK &&
          (!transfer (channels[w], &proposed, sizeof (size_t), false) ||
           !transfer (channels[w], &accepted, sizeof (size_t), false))) status = CRPX_ERROR_IO;
      if (status != CRPX_OK) { psmc_private_free_particle (&received); break; }
      psmc_private_free_particle (result->particles + i);
      result->particles[i] = received;
      result->proposals += proposed;
      result->accepted += accepted;
    }
done:
  for (size_t w = 0; w < started; w++) {
    if (status != CRPX_OK) kill (processes[w], SIGTERM);
    if (channels[w]) fclose (channels[w]);
  }
  for (size_t w = 0; w < started; w++) {
    int state;
    pid_t waited;
    do { waited = waitpid (processes[w], &state, 0); } while (waited < 0 && errno == EINTR);
    if (waited < 0 || !WIFEXITED (state) || WEXITSTATUS (state)) if (status == CRPX_OK) status = CRPX_ERROR_IO;
  }
  free (processes); free (channels);
#endif
  free (streams);
  return status;
}
