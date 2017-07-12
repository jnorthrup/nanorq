#ifndef NANORQ_PRECODE_H
#define NANORQ_PRECODE_H

#include "bitmask.h"
#include "params.h"

void precode_matrix_gen(struct pparams *prm, octmat *A, /*u*/short overhead);

octmat precode_matrix_intermediate1(struct pparams *prm, octmat *A, octmat *D);
bool precode_matrix_intermediate2(octmat *M, octmat *A, octmat *D,
                                  struct pparams *prm, repair_vec *repair_bin,
                                  struct bitmask *mask, /*u*/short num_symbols,
                                  /*u*/short overhead);

octmat precode_matrix_encode(struct pparams *prm, octmat *C, /*u*/int isi);

bool precode_matrix_decode(struct pparams *prm, octmat *X,
                           repair_vec *repair_bin, struct bitmask *mask);

#endif
