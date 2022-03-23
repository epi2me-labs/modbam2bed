#ifndef _MODBAMBED_ARGS_H
#define _MODBAMBED_ARGS_H

#include <stdbool.h>

#include "common.h"

typedef struct arguments {
    const char** bam;
    char* ref;
    char* region;
    char* read_group;
    char tag_name[2];
    int tag_value;
    mod_base mod_base;
    bool mask;
    bool cpg;
    bool chh;
    bool chg;
    bool extended;
    int threads;
    int lowthreshold;
    int highthreshold;
    char* prefix;
    int hts_maxcnt;
} arguments_t;

arguments_t parse_arguments(int argc, char** argv);

#endif
