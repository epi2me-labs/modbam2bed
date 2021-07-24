#ifndef _MODBAMBED_ARGS_H
#define _MODBAMBED_ARGS_H

typedef struct arguments {
    char *bam;
    char *ref;
    char* region;
    char* read_group;
    float lowthreshold;
    float highthreshold;
} arguments_t;

arguments_t parse_arguments(int argc, char** argv);

#endif
