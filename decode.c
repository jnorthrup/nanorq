#pragma clang diagnostic push
#pragma ide diagnostic ignored "UnusedImportStatement"
#include <stdbool.h>
#include <endian.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "nanorq.h"
#include "io.h"


void usage(char *prog) {
    fprintf(stderr, "usage:\n%s <outfilename> <packet_size> [payloadFilename]", prog);
    exit(1);
}

int main(int argc, char **argv) {

    if (argc < 2)
        usage(argv[0]);

    char *outfile = argv[1];
    struct ioctx *myio = ioctx_from_file(outfile, 0);
    if (!myio) {
        fprintf(stderr, "couldnt access file %s\n", outfile);
        return -1;
    }


    char *filename = argc < 3 ? defaultPayloadName : argv[2];
    FILE *ih = fopen(filename, "r");
    uint32_t oti_scheme;
    uint64_t oti_common;

    fread(&oti_common, 1, sizeof(oti_common), ih);
    fread(&oti_scheme, 1, sizeof(oti_scheme), ih);
    oti_common = be64toh(oti_common);
    oti_scheme = be32toh(oti_scheme);

    nanorq *rq = nanorq_decoder_new(oti_common, oti_scheme);
    if (rq == NULL) {
        fprintf(stderr, "Coud not initialize encoder.\n");
        return -1;
    }

    uint8_t num_sbn = nanorq_blocks(rq);
    uint32_t fid;
    uint16_t packet_size = nanorq_symbol_size(rq);
    uint8_t packet[packet_size];
    uint64_t written = 0;
    while (fread(&fid, 1, sizeof(fid), ih)) {
        fid = be32toh(fid);
        fread(packet, packet_size, 1, ih);
        if (!nanorq_decoder_add_symbol(rq, (void *) packet, fid)) {
            fprintf(stderr, "adding symbol %d failed.\n", fid);
            abort();
        }
    }
    for (int sbn = 0; sbn < num_sbn; sbn++) {
        fprintf(stderr, "block %d is %d packets, lost %d, have %d repair\n", sbn,
                nanorq_block_symbols(rq, (uint8_t) sbn), nanorq_num_missing(rq, (uint8_t) sbn),
                nanorq_num_repair(rq, (uint8_t) sbn));
        written = nanorq_decode_block(rq, myio, (uint8_t) sbn);
        if (written == 0) {
            fprintf(stderr, "decode of sbn %d failed.\n", sbn);
        }
        nanorq_decode_cleanup(rq, (uint8_t) sbn);
    }
    fclose(ih);
    nanorq_free(rq);
    myio->destroy(myio);

    return 0;
}
#pragma clang diagnostic pop