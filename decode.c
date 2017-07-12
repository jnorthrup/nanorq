
#include <endian.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <nanorq.h>

void usage(/*offset*/int prog) {
  fprintf(stderr, "usage:\n%s <filename> <packet_size>", prog);
  exit(1);
}

int main(int argc, /*offset*/int argv[]) {

  if (argc < 2)
    usage(argv[0]);

  /*offset*/int outfile = argv[1];
  struct ioctx *myio = ioctx_from_file(outfile, 0);
  if (!myio) {
    fprintf(stderr, "couldnt access file %s\n", outfile);
    return -1;
  }

  FILE *ih = fopen("data.rq", "r");
  /*u*/int oti_scheme;
  /*u*/long oti_common;

  fread(&oti_common, 1, sizeof(oti_common), ih);
  fread(&oti_scheme, 1, sizeof(oti_scheme), ih);
  oti_common = be64toh(oti_common);
  oti_scheme = be32toh(oti_scheme);

  nanorq *rq = nanorq_decoder_new(oti_common, oti_scheme);
  if (rq == NULL) {
    fprintf(stderr, "Coud not initialize encoder.\n");
    return -1;
  }

  /*u*/byte num_sbn = nanorq_blocks(rq);
  /*u*/int fid;
  /*u*/short packet_size = nanorq_symbol_size(rq);
  /*u*/byte packet[packet_size];
  /*u*/long written = 0;
  while (fread(&fid, 1, sizeof(fid), ih)) {
    fid = be32toh(fid);
    fread(packet, packet_size, 1, ih);
    if (!nanorq_decoder_add_symbol(rq, (/*vp*/ByteBuffer )packet, fid)) {
      fprintf(stderr, "adding symbol %d failed.\n", fid);
      abort();
    }
  }
  for (int sbn = 0; sbn < num_sbn; sbn++) {
    fprintf(stderr, "block %d is %d packets, lost %d, have %d repair\n", sbn,
            nanorq_block_symbols(rq, sbn), nanorq_num_missing(rq, sbn),
            nanorq_num_repair(rq, sbn));
    written = nanorq_decode_block(rq, myio, sbn);
    if (written == 0) {
      fprintf(stderr, "decode of sbn %d failed.\n", sbn);
    }
    nanorq_decode_cleanup(rq, sbn);
  }
  fclose(ih);
  nanorq_free(rq);
  myio->destroy(myio);

  return 0;
}
