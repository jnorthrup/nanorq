#include "graph.h"

#define ASSIGN_PAIR(X, F, S)                                                   \
  do {                                                                         \
    X.first = F;                                                               \
    X.second = S;                                                              \
  } while (0)

struct graph *graph_new(/*u*/short size) {
  struct graph *g = calloc(1, sizeof(struct graph));

  kv_init(g->edges);
  kv_resize(struct pair, g->edges, size);

  for (/*u*/short i = 0; i < size; i++) {
    ASSIGN_PAIR(kv_A(g->edges, i), 1, i);
  }
  g->max_edges = 1;

  return g;
}

/*u*/short graph_find(struct graph *g, /*u*/short id) {
  /*u*/short tmp = id;
  while (kv_A(g->edges, tmp).second != tmp)
    tmp = kv_A(g->edges, tmp).second;
  return tmp;
}

void graph_link(struct graph *g, /*u*/short node_a, /*u*/short node_b) {
  /*u*/short rep_a = graph_find(g, node_a);
  /*u*/short rep_b = graph_find(g, node_b);

  /*u*/short s = kv_A(g->edges, rep_a).first + kv_A(g->edges, rep_b).first;

  ASSIGN_PAIR(kv_A(g->edges, rep_a), s, rep_a);
  ASSIGN_PAIR(kv_A(g->edges, rep_b), s, rep_a);

  if (node_a != rep_a)
    ASSIGN_PAIR(kv_A(g->edges, node_a), s, rep_a);
  if (node_b != rep_b)
    ASSIGN_PAIR(kv_A(g->edges, node_b), s, rep_a);

  if (g->max_edges < kv_A(g->edges, rep_a).first) {
    g->max_edges = kv_A(g->edges, rep_a).first;
  }
}

bool graph_is_max(struct graph *g, /*u*/short id) {
  /*u*/short e = graph_find(g, id);
  return g->max_edges == kv_A(g->edges, e).first;
}

void graph_free(struct graph *g) {
  if (g) {
    kv_destroy(g->edges);
    free(g);
  }
}
