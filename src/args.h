#ifndef _MODBAMBED_ARGS_H
#define _MODBAMBED_ARGS_H

#include <stdbool.h>

typedef struct mod_base {
    char *name;
    char *abbrev;
    char base;
    char code;
} mod_base;

static const size_t n_mod_bases = 10;
static const mod_base mod_bases[] = {
    {"5-methylcytosine", "5mC", 'C', 'm'},
    {"5-hydroxymethylcytosine", "5hmC", 'C', 'h'},
    {"5-formylcytosine", "5fC", 'C', 'f'},
    {"5-carboxylcytosine", "5caC", 'C', 'c'},
    {"5-hydroxymethyluracil", "5hmU", 'T', 'g'},
    {"5-formyluracil", "5fU", 'T', 'e'},
    {"5-carboxyluracil", "5caU", 'T', 'b'},
    {"6-methyladenine", "6mA", 'A', 'a'},
    {"8-Oxoguanine", "8oxoG", 'G', 'o'},
    {"Xanthosine", "Xao", 'T', 'n'},
};
static const mod_base default_mod_base = {"5-methylcytosine", "5mC", 'C', 'm'};

typedef struct arguments {
    const char** bam;
    char* ref;
    char* region;
    char* read_group;
    char tag_name[2];
    int tag_value;
    mod_base mod_base;
    bool cpg;
    bool extended;
    int threads;
    int lowthreshold;
    int highthreshold;
} arguments_t;

arguments_t parse_arguments(int argc, char** argv);

#endif
