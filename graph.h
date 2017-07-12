#ifndef NANORQ_GRAPH_H
#define NANORQ_GRAPH_H

#include <stdbool.h>
#include <stdint.h>

#include "util.h"

struct graph {
  pair_vec edges;
  /*u*/short max_edges;
};
/*

struct graph *graph_new(//u*/short size);
void graph_link(struct graph *g, //u*/short node_a, //u*/short node_b);
bool graph_is_max(struct graph *g, //u*/short id);
//u*/short graph_find(struct graph *g, //u*/short id);
void graph_free(struct graph *g);
*/

#endif
