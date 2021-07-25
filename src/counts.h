#ifndef _MODBAMBED_COUNTS_H
#define _MODBAMBED_COUNTS_H

#include <stdbool.h>

// medaka-style feature data
typedef struct _plp_data {
    size_t buffer_cols;
    size_t n_cols;
    char* rname;
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
plp_data create_plp_data(size_t buffer_cols, const char* rname);


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
 *  @returns void
 *
 */
void print_bedmethyl(plp_data pileup, char *ref, int rstart, bool extended);


/** Generates base counts from a region of a bam.
 *
 *  @param bam_file input aligment file.
 *  @param chr bam target name.
 *  @param start start position of chr to consider.
 *  @param end end position of chr to consider.
 *  @param lowthreshold highest probability to call base as canonical.
 *  @param highthreshold lowest probablity to call base as modified.
 *  @returns a pileup data pointer.
 *
 *  The return value can be freed with destroy_plp_data.
 *
 */
plp_data calculate_pileup(
    const char *bam_file, const char *chr, int start, int end,
    const char *read_group, int lowthreshold, int highthreshold);


#endif
