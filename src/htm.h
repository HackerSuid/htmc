/* Prototypes and definitions for HTM implementation. */

#ifndef HTM_H_
#define HTM_H_ 1

#define DEFAULT_CONF_PATH "/etc/htmc.conf"

#include "codec.h"

typedef input_patterns* (*codec_cb)(void);

/* initialize the htmc library. parses the XML configuration
   file and stores the node data in an htm_conf structure. */
int init_htm (codec_cb cb);

/* htm learning and inference algorithms implemented
   procedurally */
int process_subcortical_input ();

#endif

