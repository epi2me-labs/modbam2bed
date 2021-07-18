#ifndef _MEDAKA_COUNTS_H
#define _MEDAKA_COUNTS_H

// medaka-style feature data
typedef struct _plp_data {
    size_t buffer_cols;
    size_t n_cols;
    char* rname;
    size_t *matrix;
    size_t *major;
} _plp_data;
typedef _plp_data *plp_data;


// medaka-style base encoding
static const char plp_bases[] = "acgtACGTdDmM";
static const size_t featlen = 12; // len of the above
static const size_t fwd_del = 9;  // position of D
static const size_t rev_del = 8;  // position of d


static const size_t fwd_mod = 10; // position of M
static const size_t rev_mod = 11; // position of m

// bam tag used for datatypes
static const char datatype_tag[] = "DT";

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
plp_data create_plp_data(size_t buffer_cols, char* rname);



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
 *  @returns void
 *
 */
void print_bedmethyl(plp_data pileup, char *ref, int rstart);


/** Generates medaka-style feature data in a region of a bam.
 *
 *  @param region 1-based region string.
 *  @param bam_file input aligment file.
 *  @param ref reference sequence corresponding to region.
 *  @param rstart the 0-based start of the reference.
 *  @returns a pileup data pointer.
 *
 *  The return value can be freed with destroy_plp_data.
 *
 */
plp_data calculate_pileup(
        const char *region, const char *bam_file, const char *read_group, char *ref, int rstart);


#endif
