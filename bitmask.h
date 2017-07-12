#ifndef BITMASK_H
#define BITMASK_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "kvec.h"

struct bitmask {
  kvec_t(/*u*/int) mask;
};
/*

struct bitmask *bitmask_new(int initial);
void bitmask_set(struct bitmask *bm, int id);
void bitmask_clear(struct bitmask *bm, int id);
bool bitmask_check(struct bitmask *bm, int id);
int bitmask_popcount(struct bitmask *bm);
int bitmask_gaps(struct bitmask *bm, int until);
void bitmask_free(struct bitmask *bm);
void bitmask_print(struct bitmask *bm);
*/

#endif
