#ifndef _MODBAMBED_COUNTS_H
#define _MODBAMBED_COUNTS_H

#include <stdbool.h>
#include <limits.h>

#include "common.h"

static const int _INT_MAX = INT_MAX;

// medaka-style feature data
typedef struct _plp_data {
    size_t buffer_cols;
    size_t n_cols;
    char *rname;
    size_t *matrix;
    size_t *major;
} _plp_data;
typedef _plp_data *plp_data;

typedef struct bed_buffer {
    int pos;
    bool isrev;
    size_t depth, cd, md, fd, xd, od;
} bed_buffer;

// files open for writing outputs
// this buf_size is silly, but its to work around CFFI sillyness
static const size_t _buf_size = 2;
typedef struct _output_files {
    bool multi;
    bool take_all;
    bool accumulated;
    bool cpg;
    bool chh;
    bool chg;
    FILE *fcpg;
    FILE *fchh;
    FILE *fchg;
    FILE *fcpg_acc;
    FILE *fchh_acc;
    FILE *fchg_acc;
    size_t buf_size;
    bed_buffer out_buffer[2];
    size_t motif_offsets[2];
    FILE* motif_acc_files[2];
} _output_files;
typedef _output_files *output_files;


output_files open_bed_files(char* prefix, bool cpg, bool chh, bool chg, bool accumulated);
void close_bed_files(output_files);
// reset state of buffers (to handle loci split by thread blocks)
void init_output_buffers(output_files bed_files);
void flush_output_buffers(output_files bed_files, const char* chr, bool extended, char* feature);


// Check sequences for motifs
// CpG
bool extern inline is_cpg_fwd(size_t rpos, int rlen, char* ref);
bool extern inline is_cpg_rev(size_t rpos, int rlen, char* ref);
// CHN
bool extern inline _is_chn_fwd(size_t rpos, int rlen, char* ref);
bool extern inline _is_chn_rev(size_t rpos, int rlen, char* ref);
// CHH
bool extern inline is_chh_fwd(size_t rpos, int rlen, char* ref);
bool extern inline is_chh_rev(size_t rpos, int rlen, char* ref);
// CHG
bool extern inline is_chg_fwd(size_t rpos, int rlen, char* ref);
bool extern inline is_chg_rev(size_t rpos, int rlen, char* ref);

// medaka-style base encoding - augmented with (a) modified base counts
static const char plp_bases[] = "acgtACGTdDmMfoOfFxX";  // o: "other mod", f:"filtered", x:"no call"

enum plp_index {
    rev_A, rev_C, rev_G, rev_T,
    fwd_A, fwd_C, fwd_G, fwd_T, 
    rev_del, fwd_del,
    rev_mod, fwd_mod,
    rev_other, fwd_other,
    rev_filt, fwd_filt,
    rev_nocall, fwd_nocall,
    featlen
};
static const size_t fwdbases[] = 
    {fwd_A, fwd_C, fwd_G, fwd_T, fwd_del, fwd_mod, fwd_other, fwd_filt, fwd_nocall}; 
static const size_t revbases[] = 
    {rev_A, rev_C, rev_G, rev_T, rev_del, rev_mod, rev_other, rev_filt, rev_nocall};
static const size_t numbases = featlen / 2;

// convert 16bit IUPAC (+16 for strand) to plp_bases index
// e.g. G=4 => fwd_G => plp_bases[6]
static const int num2countbase[32] = {
      -1, fwd_A, fwd_C,  -1, fwd_G,  -1, -1, -1,
   fwd_T,    -1,    -1,  -1,    -1,  -1, -1, -1,
      -1, rev_A, rev_C,  -1, rev_G,  -1, -1, -1,
   rev_T,    -1,    -1,  -1,    -1,  -1, -1, -1,
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
 *  @param bed_files output file handles (and filters).
 *  @returns void
 *
 */
void print_bedmethyl(
    plp_data pileup, char *ref, int rstart, bool extended,
    char *feature, char canon_base, output_files bed_files);


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
 *  @param mod_base a mod_base instance
 *  @param combine combine all modified bases corresponding to same canonical base as mb
 *  @param max_depth maximum depth of pileup.
 *  @returns a pileup data pointer.
 *
 *  The return value can be freed with destroy_plp_data.
 *
 */
plp_data calculate_pileup(
    const set_fsets *fsets, const char *chr, int start, int end,
    const char *read_group, const char tag_name[2], const int tag_value,
    int threshold, mod_base mb, bool combine, int max_depth);

#endif
