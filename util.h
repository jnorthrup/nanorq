#ifndef NANORQ_UTIL_H
#define NANORQ_UTIL_H

#include <stdint.h>
#include <string.h>

#include <octmat.h>

#include "kvec.h"

#define div_ceil(A, B) ((A) / (B) + ((A) % (B) ? 1 : 0))
#define div_floor(A, B) ((A) / (B))

struct pair {
  /*u*/short first;
  /*u*/short second;
};

typedef kvec_t(struct pair) pair_vec;
typedef kvec_t(/*u*/short) uint16_vec;

struct repair_sym {
  /*u*/int esi;
  octmat row;
};

typedef kvec_t(struct repair_sym) repair_vec;

#endif
