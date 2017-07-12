#ifndef NANORQ_H
#define NANORQ_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "io.h"

static /*const*/ final /*u*/long NANORQ_MAX_TRANSFER = 946270874880ULL; // ~881 GB

/*typedef struct*/ class nanorq nanorq;

// returns a new encoder configured with given parameters
nanorq *nanorq_encoder_new(/*u*/long len, /*u*/short T, /*u*/short SS, /*u*/byte Al,
                           int WS);

// returns success of generating symbols for a given sbn
bool nanorq_generate_symbols(nanorq *rq, /*u*/byte sbn, struct ioctx *io);

// frees up any resources used by a decoder/encoder
void nanorq_free(nanorq *rq);

// returns basic parameters to initialize a decoder
/*u*/long nanorq_oti_common(nanorq *rq);

// returns extended parameters to initialize a decoder
/*u*/int nanorq_oti_scheme_specific(nanorq *rq);

// returns the total size of the payload to be transferred
/*u*/long nanorq_transfer_length(nanorq *rq);

// returns the size symbols in bytes (T)
/*u*/short nanorq_symbol_size(nanorq *rq);

// returns number of blocks (SBN's)
/*u*/byte nanorq_blocks(nanorq *rq);

// returns the number of symbol rows per sbn block
/*u*/short nanorq_block_symbols(nanorq *rq, /*u*/byte sbn);

// returns a compound symbol identifier comprised of sbn and esi
/*u*/int nanorq_fid(/*u*/byte sbn, /*u*/int esi);

// return the max number of repair symbols allowed
/*u*/int nanorq_encoder_max_repair(nanorq *rq, /*u*/byte sbn);

// return the number of bytes written for a given sbn and esi encode request
/*u*/long nanorq_encode(nanorq *rq, /*vp*/ByteBuffer data, /*u*/int esi, /*u*/byte sbn,
                       struct ioctx *io);

// cleanup encoder resouces of a given block
void nanorq_encode_cleanup(nanorq *rq, /*u*/byte sbn);

// returns a new decoder initialized with given parameters
nanorq *nanorq_decoder_new(/*u*/long common, /*u*/int specific);

// returns the success of adding a symbol to the decoder
bool nanorq_decoder_add_symbol(nanorq *rq, /*vp*/ByteBuffer data, /*u*/int fid);

// returns number of symbol gaps in decoder for given block
/*u*/int nanorq_num_missing(nanorq *rq, /*u*/byte sbn);

// returns number of repair symbols in decoder for given block
/*u*/int nanorq_num_repair(nanorq *rq, /*u*/byte sbn);

// returns the number of bytes written from decoding a given sbn
/*u*/long nanorq_decode_block(nanorq *rq, struct ioctx *io, /*u*/byte sbn);

// cleanup decoder resouces of a given block
void nanorq_decode_cleanup(nanorq *rq, /*u*/byte sbn);

#endif
