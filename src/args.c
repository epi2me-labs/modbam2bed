#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <argp.h>

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "args.h"

const char *argp_program_version = "0.2.2";
const char *argp_program_bug_address = "chris.wright@nanoporetech.com";
static char doc[] = 
 "modbam2bed -- summarise a BAM with modified base tags to bedMethyl.\
 \vPositions absent from the methylation tags are assumed to be canonical. \
 Positions with modified probability between upper and lower thresholds\
 are removed from the counting process. Column 5 (\"score\") of the output\
 is calculated as the proportion of bases called as the canonical or modified\
 reference base with respect to the number of spanning reads, scaled to a\
 maximum of 1000. Column 11 is the percentage of reference base calls identified\
 as being modified.";
static char args_doc[] = " <reads.bam> <reference.fasta> > ";
static struct argp_option options[] = {
    {0, 0, 0, 0,
        "General options:"},
    {"region", 'r', "chr:start-end", 0,
        "Genomic region to process."},
    {"extended", 'e', 0, 0,
        "Output extended bedMethyl including counts of canonical, modified, and filtered bases (in that order)."},
    {"mod_base", 'm', "BASE", 0,
        "Modified base of interest, one of: 5mC, 5hmC, 5fC, 5caC, 5hmU, 5fU, 5caU, 6mA, 5oxoG, Xao."},
    {"threads", 't', "THREADS", 0,
        "Number of threads for BAM processing."},
    {0, 0, 0, 0,
        "Base filtering options:"},
    {"canon_threshold", 'a', "THRESHOLD", 0,
        "Bases with mod. probability < THRESHOLD are counted as canonical.", 2},
    {"mod_threshold", 'b', "THRESHOLD", 0,
        "Bases with mod. probability > THRESHOLD are counted as modified.", 2},
    {"cpg", 'c', 0, 0,
        "Output records filtered to CpG sited.", 2},
    {0, 0, 0, 0,
        "Read filtering options:"},
    {"read_group", 'g', "RG", 0,
        "Only process reads from given read group.", 3},
    {"tag_name", 0x100, "TN", 0,
        "Only process reads with a given tag (see --tag_value).", 3},
    {"tag_value", 0x200, "VAL", 0,
        "Only process reads with a given tag value.", 3},
    {"haplotype", 0x300, "VAL", 0,
        "Only process reads from a given haplotype. Equivalent to --tag_name HP --tag_value VAL.", 3},
    { 0 }
};

bool file_exists(char* filename) {
    struct stat st;
    return (stat(filename, &st) == 0);
}

static int tag_items = 0;
static bool tag_given = false;
static bool hp_given = false;
static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    arguments_t *arguments = state->input;
    float thresh;
    bool found = false;
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
        case 'm':
            for (size_t i = 0; i < n_mod_bases; ++i) {
                if (!strcmp(mod_bases[i].abbrev, arg)) {
                    arguments->mod_base = mod_bases[i];
                    found = true;
                    break;
                }
            }
            if (!found) {
                argp_error (state, "Unrecognised modified base type: '%s'.", arg);
            }
            break;
        case 'r':
            arguments->region = arg;
            break;
        case 'c':
            arguments->cpg = true;
            break;
        case 'e':
            arguments->extended = true;
            break;
        case 'g':
            arguments->read_group = arg;
            break;
        case 0x100:
            if (strlen(arg) > 2) {
                argp_error(state, "Tag name should be a two-letter code, received: '%s'.", arg);
            }
            memcpy(arguments->tag_name, arg, 2 *sizeof(char));
            tag_items += 1;
            tag_given = true;
            break;
        case 0x200:
            arguments->tag_value = atoi(arg);
            tag_items += 1;
            tag_given = true;
            break;
        case 0x300:
            memcpy(arguments->tag_name, "HP", 2 * sizeof(char));
            arguments->tag_value = atoi(arg);
            tag_items += 2;
            hp_given = true;
            break;
        case 't':
            arguments->threads = atoi(arg);
            break;
        case ARGP_KEY_NO_ARGS:
            argp_usage (state);
            break;
        case ARGP_KEY_ARG:
            switch (state->arg_num) {
                case 0:
                    arguments->bam = arg;
                    if (!file_exists(arg)) {
                        argp_error(state, "Cannot access BAM input file: '%s'.", arg);
                    }
                    htsFile *fp = hts_open(arg, "rb");
                    hts_idx_t *idx = sam_index_load(fp, arg);
                    sam_hdr_t *hdr = sam_hdr_read(fp);
                    if (hdr == 0 || idx == 0 || fp == 0) {
                        argp_error(state, "Failed to read .bam file '%s'.", arg);
                    }
                    hts_close(fp); hts_idx_destroy(idx); sam_hdr_destroy(hdr);
                    break;
                case 1:
                    arguments->ref = arg;
                    if (!file_exists(arg)) {
                        argp_error(state, "Cannot access reference input file: '%s'.", arg);
                    }
                    faidx_t *fai = fai_load(arg);
                    if (fai == NULL) {
                        argp_error(state, "Cannot read .fasta(.gz) file: '%s'.", arg);
                    }
                    fai_destroy(fai);
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
    args.mod_base = default_mod_base;
    args.lowthreshold = (int)(0.33 * 255);
    args.highthreshold = (int)(0.66 * 255);
    args.bam = NULL;
    args.ref = NULL;
    args.region = NULL;
    args.read_group = NULL;
    args.tag_name[0] = '\0';
    args.tag_value = -1;
    args.cpg = false;
    args.extended = false;
    args.threads = 1;
    argp_parse(&argp, argc, argv, 0, 0, &args);
    // allow CpG only for C!
    if(args.cpg) {
        if (args.mod_base.base != 'C') {
            fprintf(stderr, "ERROR: Option '--cpg' can only be used with cytosine modifications.");
            exit(1);
        }; 
    }
    if (tag_items % 2 > 0) {
        fprintf(stderr, "ERROR: Both or neither of --tag_name and --tag_value must be given.\n");
        exit(1);
    }
    if (tag_given && hp_given) {
        fprintf(stderr, "ERROR: If --haplotype is given neither of --tag_name or --tag_value should be provided.\n");
        exit(1);
    }
    if (args.highthreshold < args.lowthreshold) {
        fprintf(stderr, "ERROR: --highthreshold must be larger than --lowthreshold\n");
        exit(1);
    }
    return args;
}
