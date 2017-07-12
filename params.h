#ifndef NANORQ_PARAMS_H
#define NANORQ_PARAMS_H

#include "table2.h"
#include "util.h"
#include <stdbool.h>

struct pparams {
  /*u*/short K_padded;
  /*u*/short S;
  /*u*/short H;
  /*u*/short W;
  /*u*/short L;
  /*u*/short P;
  /*u*/short P1;
  /*u*/short U;
  /*u*/short B;
  /*u*/short J;
};

struct pparams params_init(/*u*/short symbols);
uint16_vec params_get_idxs(struct pparams *prm, /*u*/int X);

#endif
