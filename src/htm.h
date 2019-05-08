/* Public prototypes and definitions interface for HTM implementation
in ISO 1999 Standard C (C99). */

#ifndef HTM_H_
# define HTM_H_ 1

#define DEFAULT_CONF_PATH "/etc/htmc.conf"

/* C99 fixed-width data types for increased type
portability */
#include <stdint.h>

/* Renaming trick to make sure users include the htm
interface. This will cause linker errors if they don't. The
benefit is that function prototypes will be properly enforced
by the compiler. Without including the interface, C allows
the user to misuse it. NOTE: I think this will increase compile
times by some amount due to the extra token replacements. */
#define init_htm INT_init_htm
#define process_subcortical_input INT_process_subcortical_input
#define get_layer4 INT_get_layer4
#define get_htm_input_patterns INT_get_htm_input_patterns
#define mc_active_at INT_mc_active_at

/* inform C++ callers that this is C code */
#ifdef __cplusplus
extern "C" {
#endif

#include "repr.h"

/* initialize the htmc library: parses the XML configuration
file, and sets the encoder callback. */
extern int32_t
init_htm (codec_cb cb);

/* calls the HTM learning and inference algorithms on
the encoded input patterns. */
extern int32_t
run_cortical_algorithm (void);

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

#define RST_MC_SP_INDICATOR(mc) \
    ( (mc)->active_mask <<= 2 )
/* bitwise AND of active and SP-process bits */
#define MC_ACTIVE_AT(mc, t) \
    ( ( (mc)->active_mask & (1<<((t)+1)) ) && \
      ( (mc)->active_mask & 1 ) )

extern unsigned char
mc_active_at (struct minicolumn *mc, uint32_t t);

#ifdef __cplusplus
}
#endif

#endif

