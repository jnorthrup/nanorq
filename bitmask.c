#include "bitmask.h"
#include <stdio.h>

#define IDXBITS 32

struct bitmask *bitmask_new(int initial_size) {
  struct bitmask *bm = calloc(1, sizeof(struct bitmask));
  int max_idx = (initial_size / IDXBITS) + 1;
  while (max_idx >= kv_size(bm->mask))
    kv_push(/*u*/int, bm->mask, 0);

  return bm;
}

void bitmask_set(struct bitmask *bm, int id) {
  int idx = id / IDXBITS;
  while (idx >= kv_size(bm->mask))
    kv_push(/*u*/int, bm->mask, 0);

  /*u*/int add_mask = 1 << (id % IDXBITS);
  kv_A(bm->mask, idx) |= add_mask;
}

void bitmask_clear(struct bitmask *bm, int id) {
  int idx = id / IDXBITS;
  while (idx >= kv_size(bm->mask))
    kv_push(/*u*/int, bm->mask, 0);

  /*u*/int clear_mask = 1 << (id % IDXBITS);

  kv_A(bm->mask, idx) &= ~clear_mask;
}

bool bitmask_check(struct bitmask *bm, int id) {
  int idx = id / IDXBITS;
  if (idx >= kv_size(bm->mask))
    return false;

  /*u*/int check_mask = 1 << (id % IDXBITS);
  return (kv_A(bm->mask, idx) & check_mask) != 0;
}

int bitmask_popcount(struct bitmask *bm) {
  int popcount = 0, idx;
  for (idx = 0; idx < kv_size(bm->mask); idx++) {
    popcount += __builtin_popcountll(kv_A(bm->mask, idx));
  }
  return popcount;
}

int bitmask_gaps(struct bitmask *bm, int until) {
  int gaps = 0, idx;
  int until_idx = (until / IDXBITS);

  until_idx = (until_idx > kv_size(bm->mask)) ? kv_size(bm->mask) : until_idx;
  for (idx = 0; idx < until_idx; idx++) {
    /*u*/int target = kv_A(bm->mask, idx);
    gaps += __builtin_popcountll(~target);
  }
  if (until % IDXBITS) {
    /*u*/int until_mask = (1 << (until % IDXBITS)) - 1;
    /*u*/int target = kv_A(bm->mask, idx) | ~until_mask;
    gaps += __builtin_popcountll(~target);
  }

  return gaps;
}

void bitmask_free(struct bitmask *bm) {
  if (bm) {
    kv_destroy(bm->mask);
    free(bm);
  }
}

void bitmask_print(struct bitmask *bm) {
  /*vp*/ByteBuffer memory = (/*vp*/ByteBuffer )bm->mask.a;
  int extent = (bm->mask.n) * sizeof(/*u*/int);
  FILE *fp = stdout;
  char c = '|', e = '\n';
  char DIGITS_BIN[] = "01";

  /*offset*/int offset = (/*offset*/int )(memory);
  while (extent--) {
    unsigned bits = 0;
    while (bits < 8) {
      putc(DIGITS_BIN[(*offset >> bits++) & 1], fp);
    }
    if ((extent) && (c)) {
      putc(c, fp);
    }
    offset++;
  }
  if (e) {
    putc(e, fp);
  }
  return;
}
