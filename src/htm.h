/* Prototypes and definitions for HTM implementation. */

/* htmc requires some POSIX and GNU extensions, which are
not part of ISO C, such as
 - strdup()
 - secure_getenv() */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#ifndef _DEFAULT_SOURCE
#define _DEFAULT_SOURCE
#endif

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifndef HTM_H_
#define HTM_H_ 1

#define DEFAULT_CONF_PATH "/etc/htmc.conf"

/* C99 fixed-width data types */
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* researched more compact ways of storing the sdr bits
   but bitfields and bitwise operations are apparently
   even worse performance than bytes or words due to
   unpacking/deconstruction of the data. */

/* the sdr is two-dimensional, to make "horizontal/later
l"
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

struct layer6_conf
{
    uint32_t num_gcms, num_cells_gcm;
};

struct layer4_conf
{
    uint32_t height, width;
    short cells_per_col;
    char sensorimotor;
    uint32_t loc_patt_sz;
    short loc_patt_bits;

    struct columns_conf
    {
        float rec_field_sz;
        float local_activity;
        float column_complexity;
        char high_tier;
        uint32_t activity_cycle_window;
    } colconf;

};

struct htm_conf
{
    char *target;
    char allow_boosting;
    struct layer6_conf layer6conf;
    struct layer4_conf layer4conf;
};

typedef input_patterns* (*codec_cb)(void);

/* initialize the htmc library. parses the XML configuration
   file and stores the node data in an htm_conf structure. */
int32_t init_htm (codec_cb cb);

/* htm learning and inference algorithms implemented
   procedurally */
int32_t process_subcortical_input (void);

struct layer* get_layer4 (void);
input_patterns* get_htm_input_patterns (void);

/* layer functions & data structures */
struct layer
{
    struct minicolumn ***minicolumns;
    uint32_t height;
    uint32_t width;
    uint32_t inhibition_radius;
};


/* minicolumn functions & data structures */
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

unsigned char
mc_active_at (struct minicolumn *mc, uint32_t t);

#ifdef __cplusplus
}
#endif

#endif

