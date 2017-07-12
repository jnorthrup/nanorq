#ifndef NANORQ_CHOOSER_H
#define NANORQ_CHOOSER_H

#include <stdbool.h>

#include "graph.h"
#include "util.h"

struct tracking_pair {
  bool is_hdpc;
  int row_degree;
};

struct chooser {
  kvec_t(struct tracking_pair) tracking;
  pair_vec r_rows;
  bool only_two_ones;
};

/*
struct chooser chooser_init(//u*/short tp_size);
void chooser_clear(struct chooser *ch);
void chooser_add_tracking_pair(struct chooser *ch, bool is_hdpc,
                               int row_degree);

//u*/short chooser_non_zero(struct chooser *ch, octmat *A, struct graph *G,
                          //u*/short i, //u*/short sub_rows, //u*/short sub_cols);
//u*/short chooser_pick(struct chooser *ch, struct graph *G, //u*/short i,
                      //u*/short sub_rows, //u*/short non_zero);
*/

#endif
