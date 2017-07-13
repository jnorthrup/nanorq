#pragma clang diagnostic push
#pragma ide diagnostic ignored "UnusedImportStatement"

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>


#include <time.h>
#include <memory.h>
#include "nanorq.h"

#include "io.h"


void dump_esi(nanorq *rq, struct ioctx *myio, FILE *oh, uint8_t sbn,
              uint32_t esi) {
    uint32_t fid = htobe32(nanorq_fid(sbn, esi));
    uint16_t packet_size = nanorq_symbol_size(rq);
    uint8_t data[packet_size];
    uint64_t written = nanorq_encode(rq, (void *) data, esi, sbn, myio);

    if (written != packet_size) {
        fprintf(stderr, "failed to encode packet data for sbn %d esi %d.", sbn,
                esi);
        abort();
    } else {
        fwrite(&fid, 1, sizeof(fid), oh);
        fwrite(data, 1, packet_size, oh);
    }
}

void dump_block(nanorq *rq, struct ioctx *myio, FILE *oh, uint8_t sbn, double expected_loss, int OVERHEAD) {

    int overhead = OVERHEAD;

    uint32_t num_esi = nanorq_block_symbols(rq, sbn);
    int num_dropped = 0, num_rep = 0;
    for (uint32_t esi = 0; esi < num_esi; esi++) {
        float dropped = (float) rand() / (float) RAND_MAX * (float) 100.0;
        float drop_prob = (float) expected_loss;
        if (dropped < drop_prob) {
            num_dropped++;
        } else {
            dump_esi(rq, myio, oh, sbn, esi);
        }
    }
    for (uint32_t esi = num_esi; esi < num_esi + num_dropped + overhead; esi++) {
        dump_esi(rq, myio, oh, sbn, esi);
        num_rep++;
    }
    nanorq_encode_cleanup(rq, sbn);
    fprintf(stderr, "block %d is %d packets, dropped %d, created %d repair\n",
            sbn, num_esi, num_dropped, num_rep);
}

void usage(char *prog) {
    fprintf(stderr,
            "usage:\n%s <srcFilename> [packet_size|600] [payloadFilename|data.rq] [expected_loss|6.0] [overhead|5]",
            prog);
    exit(1);
}

int main(int argc, char **argv) {
    if (argc < 2)
        usage(argv[0]);

    char *srcFilename = argv[1];
    char *packetSize = argc < 3 ? "600" : argv[2];
    char *filename = argc < 4 ? defaultPayloadName : argv[3];
    double expected_loss;
    if (argc < 5) {
        expected_loss = 6.0;
    } else {
        char *dummy;
        expected_loss = strtod(argv[4], &dummy);
    }
    int overhead;
    if (argc < 6) {
        overhead = (int) 5;
    } else {
        char *dummy;
        overhead = (int) strtol(argv[5], &dummy, 10);
    }

    struct ioctx *myio = NULL;
    {
        myio = ioctx_from_file(srcFilename, 1);
        if (!myio) {
            fprintf(stderr, "couldnt access file %s\n", srcFilename);
            return -1;
        }
    }

    nanorq *rq = NULL;
    {
        size_t filesize = myio->size(myio);

        // determine chunks, symbol size, memory usage from size
        uint16_t packet_size = (uint16_t) strtol(packetSize, NULL, 10); // T
        uint8_t align = 4;
        uint16_t ss = (uint16_t) (packet_size / 2);
        uint32_t ws = (uint32_t) (packet_size * 100);

        srand((unsigned int) time(0));

        rq = nanorq_encoder_new(filesize, packet_size, ss, align, ws);

        if (rq == NULL) {
            fprintf(stderr, "Coud not initialize encoder.\n");
            return -1;
        }
    }

    uint8_t num_sbn = 0;
    {
        num_sbn = nanorq_blocks(rq);
        for (uint8_t b = 0; b < num_sbn; b++) {
            nanorq_generate_symbols(rq, b, myio);
        }
    }

    {
        uint64_t oti_common = htobe64(nanorq_oti_common(rq));
        uint32_t oti_scheme = htobe32(nanorq_oti_scheme_specific(rq));
        FILE *oh = fopen(filename, "w+");
        fwrite(&oti_common, 1, sizeof(oti_common), oh);
        fwrite(&oti_scheme, 1, sizeof(oti_scheme), oh);
        for (uint8_t sbn = 0; sbn < num_sbn; sbn++) {
            dump_block(rq, myio, oh, sbn, expected_loss, overhead);
        }
        fclose(oh);
    }

    nanorq_free(rq);
    myio->destroy(myio);

    return 0;
}

#pragma clang diagnostic pop