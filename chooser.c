#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "chooser.h"

struct chooser chooser_init(/*u*/short tp_size) {
  struct chooser ch = {0};

  kv_init(ch.tracking);
  kv_resize(struct tracking_pair, ch.tracking, tp_size);

  kv_init(ch.r_rows);
  kv_size(ch.r_rows) = 0;

  ch.only_two_ones = false;
  return ch;
}

void chooser_clear(struct chooser *ch) {
  kv_destroy(ch->r_rows);
  kv_destroy(ch->tracking);
}

void chooser_add_tracking_pair(struct chooser *ch, bool is_hdpc,
                               int row_degree) {
  struct tracking_pair tp = {is_hdpc, row_degree};
  kv_push(struct tracking_pair, ch->tracking, tp);
}

/*u*/short chooser_non_zero(struct chooser *ch, octmat *A, struct graph *G,
                          /*u*/short i, /*u*/short sub_rows, /*u*/short sub_cols) {
  /*u*/short non_zero = sub_cols + 1;

  for (/*u*/short row = 0; row < sub_rows; row++) {
    /*u*/short non_zero_tmp = 0;
    /*u*/short ones = 0;
    /*u*/short ones_idx[] = {0, 0};
    bool next_row = false;

    for (/*u*/short col = 0; col < sub_cols; col++) {
      if ((/*u*/byte)(om_A(*A, row + i, col + i)) != 0) {
        if (++non_zero_tmp > non_zero) {
          next_row = true;
          break;
        }
      }
      if ((/*u*/byte)(om_A(*A, row + i, col + i)) == 1) {
        if (++ones <= 2)
          ones_idx[ones - 1] = col;
      }
    }

    if (next_row || non_zero_tmp == 0)
      continue;

    if (non_zero == non_zero_tmp) {
      if (!ch->only_two_ones || ones == 2) {
        struct pair rp = {row, ones_idx[0]};
        kv_push(struct pair, ch->r_rows, rp);
      }
    } else {
      non_zero = non_zero_tmp;
      kv_size(ch->r_rows) = 0;
      struct pair rp = {row, ones_idx[0]};
      kv_push(struct pair, ch->r_rows, rp);
    }

    if (ones == 2) {
      if (non_zero == 2) {
        if (kv_A(ch->tracking, row).is_hdpc == 0) {
          graph_link(G, ones_idx[0], ones_idx[1]);
        }
        if (!ch->only_two_ones) {
          ch->only_two_ones = true;
          kv_size(ch->r_rows) = 0;
          struct pair rp = {row, ones_idx[0]};
          kv_push(struct pair, ch->r_rows, rp);
        }
      }
    }
  }
  return non_zero;
}

/*u*/short chooser_pick(struct chooser *ch, struct graph *G, /*u*/short i,
                      /*u*/short sub_rows, /*u*/short non_zero) {
  /*u*/short chosen = sub_rows;

  if (non_zero != 2) {
    /*u*/short min_row = sub_rows;
    /*u*/short min_row_hdpc = min_row;
    int min_degree = SIZE_MAX;
    int min_degree_hdpc = min_degree;
    for (int rp_idx = 0; rp_idx < kv_size(ch->r_rows); rp_idx++) {
      /*u*/short row = kv_A(ch->r_rows, rp_idx).first;
      if (kv_A(ch->tracking, row + i).is_hdpc) {
        if (kv_A(ch->tracking, row + i).row_degree < min_degree_hdpc) {
          min_degree_hdpc = kv_A(ch->tracking, row + i).row_degree;
          min_row_hdpc = row;
        }
      } else {
        if (kv_A(ch->tracking, row + i).row_degree < min_degree) {
          min_degree = kv_A(ch->tracking, row + i).row_degree;
          min_row = row;
        }
      }
    }
    if (min_row != sub_rows) {
      chosen = min_row;
    } else {
      chosen = min_row_hdpc;
    }
  } else {
    if (ch->only_two_ones) {
      for (int rp_idx = 0; rp_idx < kv_size(ch->r_rows); rp_idx++) {
        if (graph_is_max(G, kv_A(ch->r_rows, rp_idx).second)) {
          chosen = kv_A(ch->r_rows, rp_idx).first;
          break;
        }
      }
    }
    if (chosen == sub_rows) {
      chosen = kv_A(ch->r_rows, 0).first;
    }
  }
  return chosen;
}
