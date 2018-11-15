/* Prototypes and definitions for HTM implementation in
ISO 1999 Standard C. */

#ifndef HTM_H_
# define HTM_H_ 1

/* trigraph sequence, neat trick */
??=define DEFAULT_CONF_PATH "/etc/htmc.conf"

/* C99 fixed-width data types for increased integer
portability */
#include <stdint.h>

/* inform C++ callers that this is C code */
#ifdef __cplusplus
extern "C" {
#endif

/* I have researched more compact ways of storing the sdr
bits, but bitfields and bitwise operations are apparently
worse performance than bytes or words due to
unpacking/deconstruction of the data. */

/* the sdr is two-dimensional, to make "horizontal" and
"lateral" synapses and "vertical" columns a bit more
straightforward. */
typedef char** sdr_t;

/* data structure that describes the size of a pattern in
the HTM, which may or may not be an sdr */
typedef struct
{
    uint32_t height, width;
} pattern_sz;

/* data structure containing the various patterns within
the input provided by the external encoder */
typedef struct
{
    sdr_t sensory_pattern, location_pattern;
    pattern_sz sensory_sz, location_sz;
} input_patterns;

/* declaration of a pointer type to the encoder callback
function. */
typedef input_patterns* (*codec_cb)(void);

/* parsed configuration parameters for the layer 6 HTM
object. */
struct layer6_conf
{
    uint32_t num_gcms, num_cells_per_gcm;
};

/* parsed configuration parameters for the layer 4 HTM
object. */
struct layer4_conf
{
    uint32_t height, width;
    uint16_t cells_per_col;
    char sensorimotor;
    uint32_t loc_patt_sz;
    uint16_t loc_patt_bits;

    struct columns_conf
    {
        float rec_field_sz;
        float local_activity;
        float column_complexity;
        char high_tier;
        uint32_t activity_cycle_window;
    } colconf;

};

/* parsed configuration parameters for the HTM object. */
struct htm_conf
{
    char *target;
    char allow_boosting;
    struct layer6_conf layer6conf;
    struct layer4_conf layer4conf;
};

/* initialize the htmc library: parses the XML configuration
file, stores the node data in an htm_conf structure, and
sets the encoder callback. */
extern int32_t
init_htm (codec_cb cb);

/* calls the HTM learning and inference algorithms on
the encoded input patterns. */
extern int32_t
process_subcortical_input (void);

extern struct layer*
get_layer4 (void);

extern input_patterns*
get_htm_input_patterns (void);

/* HTM "layer" functions & data structures */
struct layer
{
    struct minicolumn ***minicolumns;
    uint32_t height, width, inhibition_radius;
};


/* HTM "minicolumn" functions & data structures */
struct minicolumn
{
    uint32_t overlap;
    /* 1 byte mask: LSB represents activity at current timestep t,
       LSB+1 t-1, etc. */
    unsigned char active_mask;
    float boost;
    struct cell *cells;
    struct synapse *proximal_dendrite_segment;
    uint32_t num_synapses;
    uint32_t input_xcent;
    uint32_t input_ycent;
    struct minicolumn **neighbors;
};

extern unsigned char
mc_active_at (struct minicolumn *mc, uint32_t t);

#ifdef __cplusplus
}
#endif

#endif

