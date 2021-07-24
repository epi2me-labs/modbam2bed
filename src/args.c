#include <stdlib.h>
#include <argp.h>

#include "args.h"

const char *argp_program_version = "0.0.4";
const char *argp_program_bug_address = "chris.wright@nanoporetech.com";
static char doc[] = 
  "modbambed -- summarised a BAM with modified base tags to bed-methyl.\
  \vPositions absent from the methylation tags are assumed to be canonical. \
  Postions with modified probability between upper and lower thresholds \
  are removed from the counting process.";
static char args_doc[] = "[options] <reads.bam> <reference.fasta>";
static struct argp_option options[] = {
    {"canon_threshold", 'a', "CANON_THRESHOLD", 0,
        "Bases with mod. probability < CANON_THRESHOLD are counted as canonical."},
    {"mod_threshold", 'b', "MOD_THRESHOLD", 0,
        "Bases with mod. probability > MOD_THRESHOLD are counted as modified."},
    {"region", 'r', "chr:start-end", 0,
        "Genomic region to process."},
    {"read_group", 'g', "READ_GROUP", 0,
        "Only process reads from given read group."},
    { 0 }
};


static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    arguments_t *arguments = state->input;
    float thresh;
    switch (key) {
        case 'a':
            thresh = atof(arg);
            if (thresh < 0 || thresh > 1.0) {
                argp_error (state, "Threshold parameter must be in (0,1), got %s", arg);
            }
            arguments->lowthreshold = (int)(thresh * 255);
            break;
        case 'b':
            thresh = atof(arg);
            if (thresh < 0 || thresh > 1.0) {
                argp_error (state, "Threshold parameter must be in (0,1), got %s", arg);
            }
            arguments->highthreshold = (int)(thresh * 255);
            break;
        case 'r':
            arguments->region = arg;
            break;
        case 'g':
            arguments->read_group = arg;
            break;
        case ARGP_KEY_NO_ARGS:
            argp_usage (state);
            break;
        case ARGP_KEY_ARG:
            switch (state->arg_num) {
                case 0:
                    arguments->bam = arg;
                    break;
                case 1:
                    arguments->ref = arg;
                    break;
                default:
                    argp_usage (state);
            }
            break;
        case ARGP_KEY_END:
            if (state->arg_num < 2)
                argp_usage (state);
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};

arguments_t parse_arguments(int argc, char** argv) {
    arguments_t args;
    args.lowthreshold = (int)(0.33 * 255);
    args.highthreshold = (int)(0.66 * 255);
    args.bam = NULL;
    args.ref = NULL;
    args.region = NULL;
    args.read_group = NULL;
    argp_parse(&argp, argc, argv, 0, 0, &args);
    return args;
}
