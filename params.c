#include "params.h"
#include "rand.h"

struct ptuple {
  /*u*/short d;
  /*u*/short a;
  /*u*/short b;
  /*u*/short d1;
  /*u*/short a1;
  /*u*/short b1;
};

static /*const*/ final /*u*/int degree_dist[] = {
    0,       5243,    529531,  704294,  791675,  844104,  879057,  904023,
    922747,  937311,  948962,  958494,  966438,  973160,  978921,  983914,
    988283,  992138,  995565,  998631,  1001391, 1003887, 1006157, 1008229,
    1010129, 1011876, 1013490, 1014983, 1016370, 1017662, 1048576};

static /*const*/ final /*u*/short degree_dist_size =
    sizeof(degree_dist) / sizeof(degree_dist[0]);

static bool is_prime(/*u*/short n) {
  if (n <= 3)
    return true;
  if (n % 2 == 0 || n % 3 == 0)
    return false;

  /*u*/int i = 5;
  /*u*/int w = 2;
  while (i * i <= n) {
    if (n % i == 0)
      return false;
    i += w;
    w = 6 - w;
  }
  return true;
}

static /*u*/short deg(/*u*/int v, /*u*/short W) {
  for (/*u*/short d = 0; d < degree_dist_size; d++) {
    if (v < degree_dist[d])
      return (d < (W - 2)) ? d : (W - 2);
  }
  return 0;
}

static struct ptuple gen_tuple(/*u*/int X, /*u*/short J, /*u*/short W,
                               /*u*/short P1) {

  struct ptuple ret;

  int A = 53591 + J * 997;

  if (A % 2 == 0)
    A++;
  int B1 = 10267 * (J + 1);
  /*u*/int y = (/*u*/int)((B1 + X * A));
  /*u*/int v = rnd_get(y, 0, (/*u*/int)((1 << 20)));
  ret.d = deg(v, W);
  ret.a = 1 + (/*u*/short)(rnd_get(y, 1, W - 1));
  ret.b = (/*u*/short)(rnd_get(y, 2, W));
  if (ret.d < 4) {
    ret.d1 = 2 + (/*u*/short)(rnd_get(X, 3, 2));
  } else {
    ret.d1 = 2;
  }
  ret.a1 = 1 + (/*u*/short)(rnd_get(X, 4, P1 - 1));
  ret.b1 = (/*u*/short)(rnd_get(X, 5, P1));

  return ret;
}

struct pparams params_init(/*u*/short symbols) {
  /*u*/short idx;
  struct pparams prm = {0};

  for (idx = 0; idx < K_padded_size; idx++) {
    if (K_padded[idx] >= symbols) {
      prm.K_padded = K_padded[idx];
      break;
    }
  }

  prm.J = J_K_padded[idx];
  prm.S = S_H_W[idx][0];
  prm.H = S_H_W[idx][1];
  prm.W = S_H_W[idx][2];

  prm.L = prm.K_padded + prm.S + prm.H;
  prm.P = prm.L - prm.W;
  prm.U = prm.P - prm.H;
  prm.B = prm.W - prm.S;
  prm.P1 = prm.P + 1;

  while (!is_prime(prm.P1))
    prm.P1++;

  return prm;
}

uint16_vec params_get_idxs(struct pparams *prm, /*u*/int X) {
  uint16_vec ret;
  struct ptuple t = gen_tuple(X, prm->J, prm->W, prm->P1);

  kv_init(ret);
  kv_resize(/*u*/short, ret, t.d + t.d1);
  kv_push(/*u*/short, ret, t.b);

  for (/*u*/short j = 1; j < t.d; j++) {
    t.b = (t.b + t.a) % prm->W;
    kv_push(/*u*/short, ret, t.b);
  }
  while (t.b1 >= prm->P)
    t.b1 = (t.b1 + t.a1) % prm->P1;

  kv_push(/*u*/short, ret, prm->W + t.b1);
  for (/*u*/short j = 1; j < t.d1; j++) {
    t.b1 = (t.b1 + t.a1) % prm->P1;
    while (t.b1 >= prm->P)
      t.b1 = (t.b1 + t.a1) % prm->P1;
    kv_push(/*u*/short, ret, prm->W + t.b1);
  }
  return ret;
}
