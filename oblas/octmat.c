#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <stdint.h>
#include <memory.h>

#include "octmat.h"
#include "octtables.h"

#include "oblas.h"

void om_resize(octmat *v, uint16_t r, uint16_t c) {
    v->rows = r;
    v->cols = c;

    void *aligned = NULL;
    v->cols_al = (uint16_t) ((c / OCTMAT_ALIGN + ((c % OCTMAT_ALIGN) ? 1 : 0)) * OCTMAT_ALIGN);

    if (posix_memalign(&aligned, OCTMAT_ALIGN, r * v->cols_al) != 0) {
        exit(ENOMEM);
    }
    ozero(aligned, 0, v->cols_al * r);
    v->data = (uint8_t *) aligned;
}

void om_copy(octmat *v1, octmat *v0) {
    v1->rows = v0->rows;
    v1->cols = v0->cols;
    v1->cols_al = v0->cols_al;

    if (!v1->data) {
        void *aligned = NULL;
        if (posix_memalign(&aligned, OCTMAT_ALIGN, v0->rows * v0->cols_al) != 0)
            exit(ENOMEM);
        v1->data = (uint8_t *) aligned;
    }
    memcpy(v1->data, v0->data, v0->rows * v0->cols_al);
}

void om_destroy(octmat *v) {
    v->rows = 0;
    v->cols = 0;
    v->cols_al = 0;
    free(v->data);
    v->data = NULL;
}

void om_print(FILE *stream, octmat m) {
    fprintf(stream, "[%ux%u]\n", m.rows, m.cols);
    for (int i = 0; i < m.rows; i++) {
        fprintf(stream, "|%3d", om_A(m, i, 0));
        for (int j = 1; j < m.cols; j++) {
            fprintf(stream, ", %3d", om_A(m, i, j));
        }
        fprintf(stream, "  |\n");
    }
}
