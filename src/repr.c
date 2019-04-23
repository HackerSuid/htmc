#include "repr.h"

/* user wants their 2D structure, but it isn't
literally stored that way. it's just a large
enough 1D array of ints to store the bits. */
repr_t *
new_repr(uint32_t r, uint32_t c)
{
    repr_t *repr = (repr_t *)calloc(
        sizeof(repr_t), sizeof(repr_t));

    repr->repr = (uint32_t *)calloc(
        INT_LEN(r, c), SZ);
    repr->rows = r;
    repr->cols = c;

    return repr;
}

void
print_repr(repr_t *r)
{
    uint32_t b = r->rows*r->cols;
    printf("%u bits in %u ints\n",
        b, (uint32_t)(INT_LEN(r->rows, r->cols)));
    int i;
    for (i=0; i<b; i++) {
        if (i>0 && i%SZ==0)
            printf("\n");
        printf("%u",
            (r->repr[i/SZ] & ((1 << (i%SZ)))) >> (i%SZ)
        );
    }
    printf("\n");
}

