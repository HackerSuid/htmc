#ifndef PARSE_CONF_H_
#define PARSE_CONF_H_ 1

#include <libxml/tree.h>

#include "htm.h"

/* update the conf data pointers with addresses of new
   sublayer conf struct */
/*TODO delete this function later
void update_sublayer_conf_attrs(
    struct sublayer_conf *conf,
    unsigned long conf_nodes
);*/

/* walk through the Htm node linked list and count the
   sublayer nodes. */
unsigned long cnt_htm_sublayer_nodes (xmlNodePtr);

typedef enum c_type
{
    BOOLEAN=1,
    STRING=2,
    ULONG=3,
    FLOAT=4,
    SHORT=5
} c_type;

typedef struct xml_element
{
    char *name;
    void *conf_data;
    c_type type;
    char required;
} xml_el;

int parse_htm_conf (void);

#endif

