#ifndef _MODBAMBED_BAMITER_H
#define _MODBAMBED_BAMITER_H

#include <stdbool.h>
#include "htslib/sam.h"

// a wrapper to keep track of a bam file, its index and header
typedef struct {
    htsFile *fp;
    hts_idx_t *idx;
    sam_hdr_t *hdr;
} bam_fset;

// a collection of such things
typedef struct set_fsets {
    bam_fset **fsets;
    size_t n;
} set_fsets;

// a pool of file sets for multithreaded processing
typedef struct {
    set_fsets **a;
    int cap;
    int top;
    pthread_mutex_t mu;
} fspool;

// parameters for bam iteration
typedef struct {
    htsFile *fp;
    hts_idx_t *idx;
    sam_hdr_t *hdr;
    hts_itr_t *iter;
    int min_mapQ;
    char tag_name[2];
    int tag_value;
    bool keep_missing;
    const char *read_group;
} mplp_data;


// Initialise BAM file, index and header structures
bam_fset* create_bam_fset(const char* fname, const char* ref_file);

// Destroy BAM file, index and header structures
void destroy_bam_fset(bam_fset* fset);

// Initialise multiple BAM filesets
set_fsets *create_filesets(const char **bams, const char* ref_file);

// Destroy multiple BAM filesets
void destroy_filesets(set_fsets *s);


// Create a pool of filesets for use in multithreaded processing
fspool *fspool_create(const char **bam_files, int nworkers, const char* ref_file);

// Destroy a pool of filesets
void fspool_destroy(fspool *p);

// Acquire a fileset from the pool
set_fsets *fspool_acquire(fspool *p);

// Release a fileset back to the pool
void fspool_release(fspool *p, set_fsets *fs);


/** Set up a bam file for reading (filtered) records.
 *
 *  @param bam_fset A BAM fileset from create_bam_fset
 *  @param chr bam target name.
 *  @param start start position of chr to consider.
 *  @param end end position of chr to consider.
 *  @param read_group by which to filter alignments.
 *  @param tag_name by which to filter alignments.
 *  @param tag_value associated with tag_name.
 *  @param min_mapQ minimum mapping quality of reads.
 *
 *  The return value can be freed with destroy_bam_iter_data.
 *
 */
mplp_data *create_bam_iter_data(
    const bam_fset* fset, const char *chr, int start, int end,
    const char *read_group, const char tag_name[2], const int tag_value,
    const int min_mapQ);


/** Clean up auxiliary bam reading data.
 *
 *  @param data auxiliary structure to clean.
 *
 */
void destroy_bam_iter_data(mplp_data *data);


/** Read a bam record.
 *
 *  @param data an mplp_data encoding the bam file to read with filter options.
 *  @param b output pointer.
 *
 */
int read_bam(void *data, bam1_t *b);


/** Create an map of query position to reference position
 *
 *  @param b alignment record
 *
 *  The length of the returned array is b->core->l_qlen.
 */
int *qpos2rpos(bam1_t *b);

#endif
