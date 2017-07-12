#pragma once

#include "table2.h"
#include "util.h"
#include <stdbool.h>

struct pparams {
  uint16_t K_padded;
  uint16_t S;
  uint16_t H;
  uint16_t W;
  uint16_t L;
  uint16_t P;
  uint16_t P1;
    uint16_t B;
  uint16_t J;
};

struct pparams params_init(uint16_t symbols);
uint16_vec params_get_idxs(struct pparams *prm, uint32_t X);
