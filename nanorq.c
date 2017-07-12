#include <stdio.h>

#include "nanorq.h"
#include "precode.h"

struct oti_common {
  int F;    /* input size in bytes */
  /*u*/short T;  /* the symbol size in octets, which MUST be a multiple of Al */
  /*u*/short Al; /* byte alignment, 0 < Al <= 8, 4 is recommended */
  /*u*/short SS; /* sub symbol size (multiple of alignment) */
  int WS;   /* max sub block size decode */
};

struct oti_scheme {
  /*u*/short Z; /* number of source blocks */
  /*u*/int N; /* number of sub-blocks in each source block */
  int Kt;  /* the total number of symbols required to represent input */
};

struct partition {
  /*u*/short IL; /* size of long blocks */
  /*u*/short IS; /* size of short blocks*/
  /*u*/short JL; /* number of long blocks */
  /*u*/short JS; /* number of short blocks */
};

struct source_block {
  int sbloc;
  int part_tot;
  struct partition part;
  /*u*/short al;
};

struct encoder_core {
  /*u*/byte sbn;
  /*u*/short num_symbols;
  /*u*/short symbol_size;
  struct pparams prm;
  octmat symbolmat;
};

struct decoder_core {
  /*u*/byte sbn;
  /*u*/short num_symbols;
  /*u*/short symbol_size;
  struct pparams prm;
  octmat symbolmat;
  repair_vec repair_bin;
  struct bitmask *mask;
};

struct nanorq {
  struct oti_common common;
  struct oti_scheme scheme;

  struct partition src_part; /* (KL, KS, ZL, ZS) = Partition[Kt, Z] */
  struct partition sub_part; /* (TL, TS, NL, NS) = Partition[T/Al, N] */

  struct encoder_core *encoders[UINT8_MAX];
  struct decoder_core *decoders[UINT8_MAX];
};

static struct oti_scheme gen_scheme_specific(struct oti_common *common) {
  struct oti_scheme ret = {0};
  int N_max = (int)(div_floor(common->T, common->SS)), n;
  ret.Kt = div_ceil(common->F, common->T);

  uint16_vec KL;
  kv_init(KL);
  kv_resize(/*u*/short, KL, N_max);
  memset(KL.a, 0, kv_size(KL));

  for (n = 1; n <= N_max; n++) {
    int KL_max =
        (int)common->WS /
        ((int)common->Al * div_ceil(common->T, (int)(common->Al * n)));
    if (KL_max > UINT16_MAX)
      KL_max = K_max;
    /*u*/short idx;
    for (idx = 0; idx < K_padded_size; idx++) {
      if (K_padded[idx] > KL_max)
        break;
    }
    kv_push(/*u*/short, KL, K_padded[idx == 0 ? 0 : (idx - 1)]);
  }

  int Z_tmp = (int)(div_ceil(ret.Kt, kv_A(KL, N_max - 1)));
  if (Z_tmp > (UINT8_MAX + 1)) {
    return ret;
  }

  ret.Z = Z_tmp;
  /*u*/short tmp = (/*u*/short)(div_ceil(ret.Kt, ret.Z));
  for (n = 0; n < kv_size(KL); n++) {
    if (tmp <= kv_A(KL, n)) {
      ret.N = (/*u*/short)(n + 1);
      break;
    }
  }
  kv_destroy(KL);

  return ret;
}

static struct partition fill_partition(int I, /*u*/short J) {
  struct partition p = {0, 0, 0, 0};
  if (J == 0)
    return p;
  p.IL = (/*u*/short)(div_ceil(I, J));
  p.IS = (/*u*/short)(div_floor(I, J));
  p.JL = (/*u*/short)(I - p.IS * J);
  p.JS = J - p.JL;

  if (p.JL == 0)
    p.IL = 0;
  return p;
}

static struct source_block get_source_block(nanorq *rq, /*u*/byte sbn,
                                            /*u*/short symbol_size) {
  struct source_block ret;
  ret.part = rq->sub_part;
  ret.al = rq->common.Al;
  ret.sbloc = 0;
  ret.part_tot = rq->sub_part.IL * rq->sub_part.JL;

  if (sbn < rq->src_part.JL) {
    ret.sbloc = sbn * rq->src_part.IL * symbol_size;
    return ret;
  } else if (sbn - rq->src_part.JL < rq->src_part.JS) {
    ret.sbloc = (rq->src_part.IL * rq->src_part.JL) * symbol_size +
                (sbn - rq->src_part.JL) * rq->src_part.IS * symbol_size;
    return ret;
  }

  return ret;
}

static int get_symbol_offset(struct source_block *blk, int pos,
                                /*u*/short K, /*u*/int symbol_id) {
  int i;

  if (pos < blk->part_tot) {
    int sub_blk_id = pos / blk->part.IL;
    i = blk->sbloc + sub_blk_id * K * blk->part.IL + symbol_id * blk->part.IL +
        pos % blk->part.IL;
  } else {
    int pos_part2 = pos - blk->part_tot;
    int sub_blk_id = pos_part2 / blk->part.IS;
    i = blk->sbloc + (blk->part_tot * K) + sub_blk_id * K * blk->part.IS +
        symbol_id * blk->part.IS + pos_part2 % blk->part.IS;
  }

  return i * blk->al;
}

static struct encoder_core *nanorq_block_encoder(nanorq *rq, /*u*/byte sbn) {
  /*u*/short num_symbols = nanorq_block_symbols(rq, sbn);
  /*u*/short symbol_size = rq->common.T / rq->common.Al;

  if (rq->encoders[sbn])
    return rq->encoders[sbn];

  if (num_symbols == 0 || symbol_size == 0)
    return NULL;

  struct encoder_core *enc = calloc(1, sizeof(struct encoder_core));
  enc->sbn = sbn;
  enc->num_symbols = num_symbols;
  enc->symbol_size = symbol_size;
  enc->prm = params_init(num_symbols);
  rq->encoders[sbn] = enc;
  return enc;
}

bool nanorq_generate_symbols(nanorq *rq, /*u*/byte sbn, struct ioctx *io) {
  octmat A = OM_INITIAL, D = OM_INITIAL;

  struct encoder_core *enc = nanorq_block_encoder(rq, sbn);
  struct pparams *prm = NULL;

  if (enc == NULL)
    return false;

  if (enc->symbolmat.rows > 0)
    return true;

  prm = &enc->prm;
  precode_matrix_gen(prm, &A, 0);

  om_resize(&D, prm->K_padded + prm->S + prm->H,
            enc->symbol_size * rq->common.Al);

  /*u*/short row = 0, col = 0;
  for (row = 0; row < prm->S + prm->H; row++) {
    for (col = 0; col < D.cols; col++) {
      om_A(D, row, col) = 0;
    }
  }

  struct source_block blk = get_source_block(rq, sbn, enc->symbol_size);
  for (; row < prm->S + prm->H + enc->num_symbols; row++) {
    /*u*/int symbol_id = row - (prm->S + prm->H);
    col = 0;
    for (/*u*/short i = 0; i < enc->symbol_size;) {
      int offset = get_symbol_offset(&blk, i, enc->num_symbols, symbol_id);
      /*u*/short sublen = (i < blk.part_tot) ? blk.part.IL : blk.part.IS;
      /*u*/short stride = sublen * rq->common.Al;
      /*u*/byte buf[stride];
      i += sublen;

      int got = 0;
      if (io->seek(io, offset)) {
        got = io->read(io, buf, stride);
      }
      for (int byte = 0; byte < got; byte++) {
        om_A(D, row, col++) = buf[byte];
      }
      for (int byte = got; byte < stride; byte++) {
        om_A(D, row, col++) = 0;
      }
    }
  }

  for (; row < D.rows; row++) {
    for (/*u*/short col = 0; col < D.cols; col++)
      om_A(D, row, col) = 0;
  }

  enc->symbolmat = precode_matrix_intermediate1(prm, &A, &D);
  if (enc->symbolmat.rows == 0)
    return false;

  om_destroy(&A);
  om_destroy(&D);

  return true;
}

nanorq *nanorq_encoder_new(/*u*/long len, /*u*/short T, /*u*/short SS, /*u*/byte Al,
                           int WS) {
  nanorq *rq = NULL;

  if (T == 0 || Al == 0 || T < Al || T % Al != 0 || SS < Al || (SS % Al) != 0 ||
      SS > T) {
    return NULL;
  }

  rq = calloc(1, sizeof(nanorq));
  rq->common.F = len;
  rq->common.T = T;
  rq->common.Al = Al;
  rq->common.SS = SS;
  rq->common.WS = WS;

  rq->scheme = gen_scheme_specific(&rq->common);

  if (rq->scheme.Z == 0 || rq->scheme.N == 0 ||
      div_ceil(rq->scheme.Kt, rq->scheme.Z) > K_max) {
    free(rq);
    return NULL;
  }

  rq->src_part = fill_partition(rq->scheme.Kt, rq->scheme.Z);
  rq->sub_part = fill_partition(rq->common.T / rq->common.Al, rq->scheme.N);

#ifdef NANORQ_DEBUG
  fprintf(stderr, "T: %06d SS: %06d AL: %d WS: %06lu\n", T, SS, Al, WS);
  fprintf(stderr, "Z: %06d N : %06d       Kt: %06lu\n", rq->scheme.Z,
          rq->scheme.N, rq->scheme.Kt);

  fprintf(stderr, "P1 %dx%d P2 %dx%d\n", rq->src_part.JL, rq->src_part.IL,
          rq->src_part.JS, rq->src_part.IS);
#endif

  return rq;
}

void nanorq_free(nanorq *rq) {
  /*u*/byte num_sbn = nanorq_blocks(rq);
  if (rq) {
    for (/*u*/byte sbn = 0; sbn < num_sbn; sbn++) {
      nanorq_encode_cleanup(rq, sbn);
      nanorq_decode_cleanup(rq, sbn);
    }
    free(rq);
  }
}

/*u*/long nanorq_oti_common(nanorq *rq) {
  /*u*/long ret = 0;
  ret = rq->common.F << 24; /* transfer length */
  ret |= rq->common.T;      /* symbol size */

  return ret;
}

/*u*/int nanorq_oti_scheme_specific(nanorq *rq) {
  /*u*/int ret = 0;
  rq->scheme.Z %= (UINT8_MAX + 1);
  rq->scheme.N %= (UINT16_MAX + 1);
  ret = rq->scheme.Z << 24; /* number of source blocks */
  ret |= rq->scheme.N << 8; /* number of sub-blocks */
  ret |= rq->common.Al;     /* symbol alignment */

  return ret;
}

/*u*/int nanorq_fid(/*u*/byte sbn, /*u*/int esi) {
  /*u*/int ret = (/*u*/int)(sbn) << 24;
  ret += esi % (/*u*/int)(1 << 24);
  return ret;
}

/*u*/long nanorq_transfer_length(nanorq *rq) { return rq->common.F; }

/*u*/short nanorq_symbol_size(nanorq *rq) { return rq->common.T; }

nanorq *nanorq_decoder_new(/*u*/long common, /*u*/int scheme) {
  /*u*/long F = common >> 24;
  /*u*/short T = common & 0xffff;

  nanorq *rq = NULL;

  if (F > NANORQ_MAX_TRANSFER)
    return NULL;

  rq = calloc(1, sizeof(nanorq));

  rq->common.F = F;
  rq->common.T = T;

  rq->scheme.Z = (scheme >> 24) & 0xff;
  rq->scheme.N = (scheme >> 8) & 0xffff;
  rq->common.Al = scheme & 0xff;
  rq->scheme.Kt = div_ceil(rq->common.F, rq->common.T);

  if (rq->scheme.Z == 0)
    rq->scheme.Z = (UINT8_MAX + 1);
  if (rq->scheme.N == 0)
    rq->scheme.N = (UINT16_MAX + 1);

  if (rq->common.T < rq->common.Al || rq->common.T % rq->common.Al != 0 ||
      div_ceil(div_ceil(rq->common.F, rq->common.T), rq->scheme.Z) > K_max) {
    free(rq);
    return NULL;
  }

  rq->src_part = fill_partition(rq->scheme.Kt, rq->scheme.Z);
  rq->sub_part = fill_partition(rq->common.T / rq->common.Al, rq->scheme.N);

#ifdef NANORQ_DEBUG
  fprintf(stderr, "T: %06d SS: XXXXXX AL: %d WS: XXXXXX\n", rq->common.T,
          rq->common.Al);
  fprintf(stderr, "Z: %06d N : %06d       Kt: %06lu\n", rq->scheme.Z,
          rq->scheme.N, rq->scheme.Kt);

  fprintf(stderr, "P1 %dx%d P2 %dx%d\n", rq->src_part.JL, rq->src_part.IL,
          rq->src_part.JS, rq->src_part.IS);
#endif

  return rq;
}

/*
+   num(0) JL
+   num(1) JS
+   size(0) IL
+   size(1) IS
*/

/*u*/short nanorq_block_symbols(nanorq *rq, /*u*/byte sbn) {
  if (sbn < rq->src_part.JL)
    return rq->src_part.IL;
  if (sbn - rq->src_part.JL < rq->src_part.JS)
    return rq->src_part.IS;
  return 0;
}

/*u*/int nanorq_encoder_max_repair(nanorq *rq, /*u*/byte sbn) {
  return (/*u*/int)((1 << 20) - nanorq_block_symbols(rq, sbn));
}

/*u*/byte nanorq_blocks(nanorq *rq) {
  return (/*u*/byte)(rq->src_part.JL + rq->src_part.JS);
}

/*u*/long nanorq_encode(nanorq *rq, /*vp*/ByteBuffer data, /*u*/int esi, /*u*/byte sbn,
                       struct ioctx *io) {
  /*u*/long written = 0;

  struct encoder_core *enc = nanorq_block_encoder(rq, sbn);
  if (enc == NULL)
    return 0;

  if (esi < enc->num_symbols) {
    struct source_block blk = get_source_block(rq, sbn, enc->symbol_size);
    /*u*/ByteBuffer dst = ((/*u*/ByteBuffer )data);
    for (/*u*/short i = 0; i < enc->symbol_size;) {
      int offset = get_symbol_offset(&blk, i, enc->num_symbols, esi);
      /*u*/short sublen = (i < blk.part_tot) ? blk.part.IL : blk.part.IS;
      /*u*/short stride = sublen * rq->common.Al;
      /*u*/byte buf[stride];
      i += sublen;

      int got = 0;
      if (io->seek(io, offset)) {
        got = io->read(io, buf, stride);
      }
      for (int byte = 0; byte < got; byte++) {
        *dst = buf[byte];
        dst++;
        written++;
      }
      for (int byte = got; byte < stride; byte++) {
        *dst = 0;
        dst++;
        written++;
      }
    }
  } else {
    // esi is for repair symbol
    struct pparams *prm = &enc->prm;
    if (enc->symbolmat.rows == 0) {
      bool generated = nanorq_generate_symbols(rq, sbn, io);
      if (!generated)
        return 0;
    }

    /*u*/int isi = esi + (prm->K_padded - enc->num_symbols);
    octmat tmp = precode_matrix_encode(prm, &enc->symbolmat, isi);
    /*u*/ByteBuffer dst = ((/*u*/ByteBuffer )data);
    /*u*/ByteBuffer octet = om_P(tmp);
    for (/*u*/short i = 0; i < enc->symbol_size; i++) {
      for (int byte = 0; byte < rq->common.Al; byte++) {
        *dst = (octet == NULL) ? 0 : *(octet++);
        dst++;
        written++;
      }
    }
    om_destroy(&tmp);
  }
  return written;
}

void nanorq_encode_cleanup(nanorq *rq, /*u*/byte sbn) {
  if (rq->encoders[sbn]) {
    struct encoder_core *enc = rq->encoders[sbn];
    om_destroy(&enc->symbolmat);
    free(enc);
    rq->encoders[sbn] = NULL;
  }
}

static struct decoder_core *nanorq_block_decoder(nanorq *rq, /*u*/byte sbn) {
  /*u*/short num_symbols = nanorq_block_symbols(rq, sbn);
  /*u*/short symbol_size = rq->common.T / rq->common.Al;

  if (rq->decoders[sbn])
    return rq->decoders[sbn];

  if (num_symbols == 0 || symbol_size == 0)
    return NULL;

  struct decoder_core *dec = calloc(1, sizeof(struct decoder_core));
  dec->sbn = sbn;
  dec->num_symbols = num_symbols;
  dec->symbol_size = symbol_size;
  dec->prm = params_init(num_symbols);
  dec->mask = bitmask_new(num_symbols);
  om_resize(&dec->symbolmat, num_symbols, symbol_size * rq->common.Al);

  rq->decoders[sbn] = dec;
  return dec;
}

bool nanorq_decoder_add_symbol(nanorq *rq, /*vp*/ByteBuffer data, /*u*/int fid) {

  /*u*/byte sbn = fid >> 24;
  /*u*/int esi = (fid & 0x00ffffff);

  struct decoder_core *dec = nanorq_block_decoder(rq, sbn);

  if (dec == NULL)
    return false;

  /*u*/short cols = dec->symbolmat.cols;

  if (esi >= (1 << 20))
    return false;

  if (bitmask_gaps(dec->mask, dec->num_symbols) == 0) {
    return true; // no gaps! no repair needed.
  }

  if (bitmask_check(dec->mask, esi))
    return true; // already got this esi

  if (esi < dec->num_symbols) {
    memcpy(om_R(dec->symbolmat, esi), data, cols);
  } else {
    struct repair_sym rs = {esi, OM_INITIAL};
    om_resize(&rs.row, 1, cols);
    memcpy(om_R(rs.row, 0), data, cols);
    kv_push(struct repair_sym, dec->repair_bin, rs);
  }
  bitmask_set(dec->mask, esi);

  return true;
}

/*u*/int nanorq_num_missing(nanorq *rq, /*u*/byte sbn) {
  /*u*/short num_symbols = nanorq_block_symbols(rq, sbn);
  struct decoder_core *dec = nanorq_block_decoder(rq, sbn);
  if (dec == NULL)
    return 0;

  return bitmask_gaps(dec->mask, num_symbols);
}

/*u*/int nanorq_num_repair(nanorq *rq, /*u*/byte sbn) {
  struct decoder_core *dec = nanorq_block_decoder(rq, sbn);
  if (dec == NULL)
    return 0;

  return kv_size(dec->repair_bin);
}

/*u*/long nanorq_decode_block(nanorq *rq, struct ioctx *io, /*u*/byte sbn) {
  /*u*/long written = 0;

  struct decoder_core *dec = nanorq_block_decoder(rq, sbn);
  struct pparams *prm = &dec->prm;
  if (dec == NULL)
    return 0;

  bool success =
      precode_matrix_decode(prm, &dec->symbolmat, &dec->repair_bin, dec->mask);
  if (!success) {
    return 0;
  }

  /*u*/short max_esi = dec->symbolmat.rows;
  /*u*/short row = 0, col = 0;
  struct source_block blk = get_source_block(rq, sbn, dec->symbol_size);
  for (; row < max_esi; row++) {
    col = 0;
    for (/*u*/short i = 0; i < dec->symbol_size;) {
      int offset = get_symbol_offset(&blk, i, max_esi, row);
      /*u*/short sublen = (i < blk.part_tot) ? blk.part.IL : blk.part.IS;
      /*u*/short stride = sublen * rq->common.Al;
      i += sublen;

      if (io->seek(io, offset)) {
        /*u*/short len = stride;
        if (offset >= rq->common.F)
          continue;
        if ((offset + stride) >= rq->common.F) {
          len = (rq->common.F - offset);
        }
        written += io->write(io, om_R(dec->symbolmat, row) + col, len);
        col += stride;
      }
    }
  }

  return written;
}

void nanorq_decode_cleanup(nanorq *rq, /*u*/byte sbn) {
  if (rq->decoders[sbn]) {
    struct decoder_core *dec = rq->decoders[sbn];
    om_destroy(&dec->symbolmat);
    if (kv_size(dec->repair_bin) > 0) {
      for (/*u*/short rs = 0; rs < kv_size(dec->repair_bin); rs++) {
        om_destroy(&(kv_A(dec->repair_bin, rs).row));
      }
      kv_destroy(dec->repair_bin);
    }
    bitmask_free(dec->mask);
    free(dec);
    rq->decoders[sbn] = NULL;
  }
}
