/* Prototypes and definitions for HTM implementation. */

#ifndef HTM_H_
#define HTM_H_ 1

#define DEFAULT_CONF_PATH "/etc/htmc.conf"

#include "codec.h"

#ifdef __cplusplus
extern "C" {
#endif

struct layer6_conf
{
    unsigned long num_gcms;
    unsigned long num_cells_gcm;
};

struct layer4_conf
{
    unsigned long height;
    unsigned long width;
    short cells_per_col;
    char sensorimotor;
    unsigned long loc_patt_sz;
    short loc_patt_bits;

    struct columns_conf
    {
        float rec_field_sz;
        float local_activity;
        float column_complexity;
        char high_tier;
        unsigned long activity_cycle_window;
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
int init_htm (codec_cb cb);

/* htm learning and inference algorithms implemented
   procedurally */
int process_subcortical_input (void);

struct layer* get_layer4 (void);
input_patterns* get_htm_input_patterns (void);

/* layer functions & data structures */
struct layer
{
    struct minicolumn ***minicolumns;
    unsigned int height;
    unsigned int width;
    unsigned int inhibition_radius;
};


/* minicolumn functions & data structures */
struct minicolumn
{
    unsigned int overlap;
    /* 1 byte mask: LSB represents activity at current timestep t,
       LSB+1 t-1, etc. */
    unsigned char active_mask;
    float boost;
    struct cell *cells;
    struct synapse *proximal_dendrite_segment;
    unsigned int num_synapses;
    unsigned int input_xcent;
    unsigned int input_ycent;
    struct minicolumn **neighbors;
};

unsigned char mc_active_at(struct minicolumn *mc, unsigned int t);

#ifdef __cplusplus
}
#endif

#endif

