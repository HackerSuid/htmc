#ifndef CODEC_H_
#define CODEC_H_ 1

/* researched more compact ways of storing the sdr bits
   but bitfields and bitwise operations are apparently
   even worse performance than bytes or words due to
   unpacking/deconstruction of the data. */

/* the sdr is two-dimensional, to make "horizontal/lateral"
   connections and "vertical" columns a bit more
   straightforward. */

typedef char** sdr_t;

typedef struct
{
    unsigned int height;
    unsigned int width;
} pattern_sz;

typedef struct
{
    sdr_t sensory_pattern;
    sdr_t location_pattern;
    pattern_sz sensory_sz;
    pattern_sz location_sz;
} input_patterns;

#endif

