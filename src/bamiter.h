#ifndef _MODBAMBED_BAMITER_H
#define _MODBAMBED_BAMITER_H

#include <stdbool.h>
#include "htslib/sam.h"

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


typedef struct {
    htsFile *fp;
    hts_idx_t *idx;
    sam_hdr_t *hdr;
} bam_fset;

typedef struct set_fsets {
    bam_fset **fsets;
    size_t n;
} set_fsets;


// Initialise BAM file, index and header structures
bam_fset* create_bam_fset(const char* fname);

// Destory BAM file, index and header structures
void destroy_bam_fset(bam_fset* fset);

// Initialise multiple BAM filesets
set_fsets *create_filesets(const char **bams);

// Destroy multiple BAM filesets
void destroy_filesets(set_fsets *s);


/** Set up a bam file for reading (filtered) records.
 *
 *  @param bam_fset A BAM fileset from create_bam_fset
 *  @param chr bam target name.
 *  @param start start position of chr to consider.
 *  @param end end position of chr to consider.
 *  @param read_group by which to filter alignments.
 *  @param tag_name by which to filter alignments.
 *  @param tag_value associated with tag_name.
 *
 *  The return value can be freed with destroy_bam_iter_data.
 *
 */
mplp_data *create_bam_iter_data(
    const bam_fset* fset, const char *chr, int start, int end,
    const char *read_group, const char tag_name[2], const int tag_value);

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
