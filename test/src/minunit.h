/*===========================================================================*/
/*
 * BEGIN ROOUNFOLD COPYRIGHT
 * RooUnfold — Unfolding library for particle-physics inverse problems
 *
 * Copyright © 2021–2025 CERN and the authors’ respective research institutions
 * Please refer to the CONTRIBUTORS file for details.
 *
 * License: BSD-3-Clause
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * END ROOUNFOLD COPYRIGHT
 */
/*===========================================================================*/

#undef NDEBUG
#ifndef _minunit_h
#define _minunit_h

#include <stdio.h>
#include "dbg.h"
#include <stdlib.h>

#define mu_suite_start() const char *message = NULL

#define mu_assert(test, message) \
   if (!(test)) {                \
      log_err(message);          \
      return message;            \
   }
#define mu_run_test(test)         \
   debug("\n-----%s", " " #test); \
   message = test();              \
   if (message)                   \
      return message;

#define RUN_TESTS(name)                       \
   int main(int, char *argv[])                \
   {                                          \
      debug("----- RUNNING: %s", argv[0]);    \
      printf("----\nRUNNING: %s\n", argv[0]); \
      const char *result = name();            \
      if (result != 0) {                      \
         printf("FAILED: %s\n", result);      \
      } else {                                \
         printf("ALL TESTS PASSED\n");        \
      }                                       \
      exit(result != 0);                      \
   }

#endif
