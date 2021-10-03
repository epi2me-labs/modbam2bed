#ifndef _MODBAMBED_COUNTS_H
#define _MODBAMBED_COUNTS_H

#include <stdbool.h>

#include "args.h"

// medaka-style feature data
typedef struct _plp_data {
    size_t buffer_cols;
    size_t n_cols;
    char *rname;
    size_t *matrix;
    size_t *major;
} _plp_data;
typedef _plp_data *plp_data;


// medaka-style base encoding - augmented with (a) modified base counts
static const char plp_bases[] = "acgtACGTdDmMfF";  // f for "filtered"
static const size_t featlen = 14; // len of the above
static const size_t fwd_del = 9;  // position of D
static const size_t rev_del = 8;  // position of d
static const size_t fwd_mod = 11; // position of M
static const size_t rev_mod = 10; // position of m
static const size_t fwd_filt = 13; // position of F
static const size_t rev_filt = 12; // position of f

// convert 16bit IUPAC (+16 for strand) to plp_bases index
static const int num2countbase[32] = {
 -1,  4,  5, -1,  6, -1, -1, -1,
  7, -1, -1, -1, -1, -1, -1, -1,
 -1,  0,  1, -1,  2, -1, -1, -1,
  3, -1, -1, -1, -1, -1, -1, -1,
};


/** Constructs a pileup data structure.
 *
 *  @param buffer_cols maximum number of pileup columns.
 *  @param rname reference name.
 *  @see destroy_plp_data
 *  @returns a plp_data pointer.
 *
 *  The return value can be freed with destroy_plp_data.
 *
 */
plp_data create_plp_data(size_t buffer_cols, const char *rname);


/** Destroys a pileup data structure.
 *
 *  @param data the object to cleanup.
 *  @returns void.
 *
 */
void destroy_plp_data(plp_data data);


/** Prints a pileup data structure.
 *
 *  @param pileup a pileup counts structure.
 *  @returns void
 *
 */
void print_pileup_data(plp_data pileup);


/** Prints a pileup data structure as bedmethyl file
 *
 *  @param pileup a pileup counts structure.
 *  @param ref reference sequence.
 *  @param rstart starting reference coordinate corresponding to ref.
 *  @param extended whether to include counts of canonical, modified and filtered bases.
 *  @param feature name to use for feature column of BED (e.g. 5mC).
 *  @param canon_base canonical base to match.
 *  @param cpg filter output to only CpG sites.
 *  @returns void
 *
 */
void print_bedmethyl(plp_data pileup, char *ref, int rstart, bool extended, char *feature, char canon_base, bool cpg);


/** Generates base counts from a region of a bam.
 *
 *  @param bam_file input aligment file.
 *  @param chr bam target name.
 *  @param start start position of chr to consider.
 *  @param end end position of chr to consider.
 *  @param read_group by which to filter alignments.
 *  @param tag_name by which to filter alignments.
 *  @param tag_value associated with tag_name
 *  @param lowthreshold highest probability to call base as canonical.
 *  @param highthreshold lowest probablity to call base as modified.
 *  @param mod_base BAM code for modified base to report. (e.g. h for 5hmC).
 *  @returns a pileup data pointer.
 *
 *  The return value can be freed with destroy_plp_data.
 *
 */
plp_data calculate_pileup(
    const char **bam_file, const char *chr, int start, int end,
    const char *read_group, const char tag_name[2], const int tag_value,
    int lowthreshold, int highthreshold, char mod_base);


/* Process and print a single region using a threadpool
 *
 * @param args program arguments.
 * @param chr reference sequence to process.
 * @param start reference coordinate to process (0-based).
 * @param end reference coordiate to process (exclusive).
 * @param ref reference sequence.
 *
 */
void process_region(arguments_t args, const char *chr, int start, int end, char *ref);

#endif
