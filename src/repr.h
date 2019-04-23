/* Interface for creating and modifying binary
arrays */
#ifndef REPR_H_
#define REPR_H_ 1

#include <stdint.h>

/* bits in an int */
#define SZ (8*sizeof(uint32_t))

#define INT_EXTRA(r, c) \
    (r)*(c) % SZ
#define INT_LEN(r, c) \
    (r)*(c)/SZ + (INT_EXTRA(r, c)? 1 : 0)

/* index of bit in int array. rep must be a pointer to repr_t */
#define BIT_IDX(rep, r ,c) \
    (((r)*(rep)->cols+(c))/SZ)
/* position of bit in int at index */
#define BIT_POS(rep, r, c) \
    (((r)*(rep)->cols+(c))%SZ)
/* perform operation if it is valid, otherwise print a
warning. */
#define VALID_MODIFY(operation, rep, r, c) \
    __builtin_expect(((r)<0||(r)>(rep)->rows-1) ||  \
                     ((c)<0||(c)>(rep)->cols-1), 0) ? \
        fprintf(stderr, "[WARN] Attempt to modify repr bit " \
            "(%u, %u), outside range (%u, %u).\n",  \
            (r), (c), (rep)->rows, (rep)->cols) : \
        (operation)

/* set bit in representation to 1. rep must be a
pointer to a repr_t to allow modifying the memory. */
#define SET_REPR_BIT(rep, r, c) \
    VALID_MODIFY( \
        (rep)->repr[BIT_IDX(rep, r, c)] |= 1<<BIT_POS(rep, r, c), \
        rep, r, c \
    )

#define CLR_REPR_BIT(rep, r, c) \
    VALID_MODIFY( \
        (rep)->repr[BIT_IDX(rep, r, c)] &= ~(1<<BIT_POS(rep, r, c)), \
        rep, r, c \
    )

/* Fast, unchecked versions */
#define SET_REPR_BIT_FAST(rep, r, c) \
        ( (rep)->repr[BIT_IDX(rep, r, c)] |= 1<<BIT_POS(rep, r, c) )

#define CLR_REPR_BIT_FAST(rep, r, c) \
        ( (rep)->repr[BIT_IDX(rep, r, c)] &= ~(1<<BIT_POS(rep, r, c)) )

#define TEST_REPR_BIT_FAST(rep, r, c) \
    ( (rep)->repr[BIT_IDX(rep, r, c)] & 1<<BIT_POS(rep, r, c) )


/* raw input patterns are represented by a binary
array of bits. */
typedef struct
{
    uint32_t rows, cols;
    uint32_t *repr;
} repr_t;


repr_t*
new_repr(uint32_t r, uint32_t c);

void
print_repr(repr_t *rep);

/* data structure containing the various patterns within
the input provided by the external encoder */
typedef struct
{
    repr_t *sensory_pattern, *location_pattern;
} input_patterns;

/* declaration of a pointer type to the encoder callback
function. */
typedef input_patterns* (*codec_cb)(void);

#endif

