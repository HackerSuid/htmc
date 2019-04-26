#include "htm.h"
#include "conf.h"
#include "parse_conf.h"
#include "utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <libxml/tree.h>

#define CONF_NODE(OBJ, FIELD, TYPE, REQ) \
    { #FIELD, &OBJ.FIELD, TYPE, REQ }

#define HTMCONF_NODE(FIELD, TYPE, REQ) \
    CONF_NODE(htmconf, FIELD, TYPE, REQ)

#define L6CONF_NODE(FIELD, TYPE, REQ) \
    CONF_NODE(htmconf.layer6conf, FIELD, TYPE, REQ)

#define L4CONF_NODE(FIELD, TYPE, REQ) \
    CONF_NODE(htmconf.layer4conf, FIELD, TYPE, REQ)

#define COLCONF_NODE(FIELD, TYPE, REQ) \
    CONF_NODE(htmconf.layer4conf.colconf, FIELD, TYPE, REQ)

static int32_t
set_conf_node_attr (xmlNodePtr node, xml_el attr);

struct htm_conf htmconf;

xml_el htm_conf_attrs[] =
{
/*
    {
        "target",
        &htmconf.target,
        STRING, 1
    },
*/
    HTMCONF_NODE(target, STRING, 1),
    HTMCONF_NODE(allow_boosting, BOOLEAN, 0)
};

xml_el layer6_conf_attrs[] =
{
    L6CONF_NODE(num_gcms, ULONG, 1),
    L6CONF_NODE(num_cells_per_gcm, ULONG, 1)
};

xml_el layer4_conf_attrs[] =
{
    L4CONF_NODE(height, ULONG, 1),
    L4CONF_NODE(width, ULONG, 1),
    L4CONF_NODE(cells_per_col, SHORT, 1),
    L4CONF_NODE(loc_patt_sz, ULONG, 1),
    L4CONF_NODE(loc_patt_bits, SHORT, 1)
};

xml_el columns_conf_attrs[] =
{
    COLCONF_NODE(rec_field_sz, FLOAT, 1),
    COLCONF_NODE(local_activity, FLOAT, 1),
    COLCONF_NODE(column_complexity, FLOAT, 1),
    COLCONF_NODE(high_tier, BOOLEAN, 1),
    COLCONF_NODE(activity_cycle_window, ULONG, 1)
};

int parse_htm_conf (void)
{
    char *conf_path=DEFAULT_CONF_PATH;
    char *conf_env=NULL;
    xmlDocPtr doc;
    xmlNodePtr node, colnode;
    xmlChar *xmlstr=NULL;
    unsigned long attr_num=0;
    unsigned long htm_conf_nodes=0;
    unsigned long layer6_conf_nodes=0;
    unsigned long layer4_conf_nodes=0;
    unsigned long col_conf_nodes=0;
    uint32_t c;

    if ((conf_env=(char *)getenv("HTM_CONF_PATH"))!=NULL)
        conf_path=conf_env;
    if (access(conf_path, R_OK)<0) {
        perror("Could not open config file for reading");
        goto fail_jmp;
    }

    INFO("[*] Parsing configuration: %s\n", conf_path);

    /* ignore whitespace "text" nodes between nodes */
    xmlKeepBlanksDefault(0);
    if ((doc=xmlParseFile(conf_path))==NULL) {
        ERR("Failed to parse %s\n", conf_path);
        goto fail_jmp;
    }
    /* parse htm root node attributes */
    if ((node = xmlDocGetRootElement(doc))==NULL) {
        ERR("Config file %s seems empty\n",
            conf_path);
        xmlFreeDoc(doc);
        goto fail_jmp;
    }
    if (xmlStrcmp(node->name, (const xmlChar *)"Htm")) {
        ERR("Htm config node missing/misplaced\n");
        goto fail_jmp;
    }

    /* htmconf is an uninitialized global, so it will be
       stored in .bss, preset with zeroes by the loader.
       so i don't need to do it myself. */
    /* memset(&htmconf, 0, sizeof(struct htm_conf)); */

    /* iterate through htm conf attrs to parse config nodes */
    htm_conf_nodes = sizeof(htm_conf_attrs)/sizeof(xml_el);
    for (c=0; c<htm_conf_nodes; c++) {
        if (set_conf_node_attr(node, htm_conf_attrs[c])) {
            ERR("Conf parsing failed at attribute %s\n",
            htm_conf_attrs[c].name);
            goto fail_jmp;
        }
    }
    
    /* parse layer4 node attributes */
    node = node->xmlChildrenNode;
    if (xmlStrcmp(node->name, (const xmlChar *)"Layer4")) {
        ERR("Layer4 config node seems missing/misplaced\n");
        goto fail_jmp;
    }

    /* iterate through layer4 conf attrs to parse config nodes */
    layer4_conf_nodes = sizeof(layer4_conf_attrs)/sizeof(xml_el);
    for (c=0; c<layer4_conf_nodes; c++) {
        if (set_conf_node_attr(node, layer4_conf_attrs[c])) {
            ERR("Conf parsing failed at attribute %s\n",
            layer4_conf_attrs[c].name);
            goto fail_jmp;
        }
    }

    /* parse columns node attributes */
    if ((colnode = node->xmlChildrenNode)==NULL) {
        ERR("Config file is incomplete.\n");
        goto fail_jmp;
    }
    if (xmlStrcmp(colnode->name, (const xmlChar *)"Minicolumns")) {
        ERR("Minicolumns config node seems missing/misplaced\n");
        goto fail_jmp;
    }
    
    /* iterate through columns conf attrs to parse config nodes */
    col_conf_nodes = sizeof(columns_conf_attrs)/sizeof(xml_el);
    for (c=0; c<col_conf_nodes; c++) {
        if (set_conf_node_attr(colnode, columns_conf_attrs[c])) {
            ERR("Conf parsing failed at attribute %s\n",
            columns_conf_attrs[c].name);
            goto fail_jmp;
        }
    }

    if (xmlstr) xmlFree(xmlstr);
    if (doc) xmlFreeDoc(doc);
    return 0;

    fail_jmp:
        if (xmlstr) xmlFree(xmlstr);
        if (doc) xmlFreeDoc(doc);
        return 1;
}

void free_sublayer_confs()
{
    /* free the htmconf layerconf linked list */
}

unsigned long cnt_xml_node_attributes(xmlNodePtr node)
{
    xmlAttr *attr=NULL;
    unsigned long attr_num=0;

    attr = node->properties;
    while (attr) {
        attr=attr->next;
        attr_num++;
    }

    return attr_num;
}

unsigned long cnt_htm_sublayer_nodes(xmlNodePtr htmnode)
{
    unsigned long nsublayers=0;
    xmlNodePtr tmpnode = htmnode->xmlChildrenNode;

    while (tmpnode) {
        nsublayers++;
        tmpnode = tmpnode->next;
    }

    return nsublayers;
}

static int32_t
set_conf_node_attr (xmlNodePtr node, xml_el attr)
{
    xmlChar *xmlstr=NULL;

    /* get the string representation of the node value */
    xmlstr = xmlGetProp(node, (const xmlChar *)attr.name);
    if (!xmlstr) {
        if (attr.required==1) {
            ERR("Config requires '%s' attribute.\n",
                attr.name);
            return 1;
        }
        return 0;
    }
    /* parse node value given expected format */
    switch (attr.type) {
        case BOOLEAN:
            /* non-false == true */
            if (xmlStrcmp(xmlstr, (const xmlChar *)"false"))
                *(char *)attr.conf_data = 1;
            break;
        case STRING:
            attr.conf_data =
                strdup((char *)xmlstr);
            break;
        case ULONG:
            *(long *)attr.conf_data =
                strtoul((char *)xmlstr, NULL, 10);
            break;
        case FLOAT:
            *(float *)attr.conf_data =
                atof((char *)xmlstr);
            break;
        case SHORT:
            *(short *)attr.conf_data =
                (short)atoi((char *)xmlstr);
            break;
        default:
            ERR("Unsupported type\n");
            xmlFree(xmlstr);
            return 1;
    }

    xmlFree(xmlstr);
    return 0;
}

