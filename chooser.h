#pragma once

struct tracking_pair {
  bool is_hdpc;
  size_t row_degree;
};

struct chooser {
  kvec_t(struct tracking_pair) tracking;
   pair_vec r_rows;
  bool only_two_ones;
};

struct chooser chooser_init(uint16_t  );
void chooser_clear(struct chooser * );
void chooser_add_tracking_pair(struct chooser * , bool  ,
                               size_t  );

uint16_t chooser_non_zero(struct chooser * , octmat *  , struct graph * ,
                          uint16_t  , uint16_t  , uint16_t  );
uint16_t chooser_pick(struct chooser * , struct graph *  , uint16_t  ,
                      uint16_t  , uint16_t  );

