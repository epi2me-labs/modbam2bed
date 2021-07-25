#define _GNU_SOURCE
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"

#include "bamiter.h"
#include "common.h"
#include "counts.h"
#include "args.h"

#define bam1_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)
#define bam1_seqi(s, i) (bam_seqi((s), (i)))
#define bam_nt16_rev_table seq_nt16_str
#define bam_nt16_table seq_nt16_table


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
plp_data create_plp_data(size_t buffer_cols, const char *rname) {
    plp_data data = xalloc(1, sizeof(_plp_data), "plp_data");
    data->buffer_cols = buffer_cols;
    data->n_cols = 0;
    data->matrix = xalloc(featlen * buffer_cols, sizeof(size_t), "matrix");
    data->major = xalloc(buffer_cols, sizeof(size_t), "major");
    data->rname = xalloc(strlen(rname) + 1, sizeof(char), "chr");
    strcpy(data->rname, rname);
    return data;
}


/** Destroys a pileup data structure.
 *
 *  @param data the object to cleanup.
 *  @returns void.
 *
 */
void destroy_plp_data(plp_data data) {
    free(data->matrix); free(data->major); free(data->rname); free(data);
}


/** Prints a pileup data structure.
 *
 *  @param pileup a pileup structure.
 *  @returns void
 *
 */
void print_pileup_data(plp_data pileup){
    fprintf(stdout, "chrom\tpos\t");
    for (size_t j = 0; j < featlen; ++j){
        fprintf(stdout, "%c\t", plp_bases[j]);
    }
    fprintf(stdout, "depth\n");
    for (size_t j = 0; j < pileup->n_cols; ++j) {
        int s = 0;
        fprintf(stdout, "%s\t%zu\t", pileup->rname, pileup->major[j]);
        for (size_t i = 0; i < featlen; ++i){
            size_t c = pileup->matrix[j * featlen + i];
            s += c;
            fprintf(stdout, "%zu\t", c);
        }
        fprintf(stdout, "%d\n", s);
    }
}


/** Prints a pileup data structure as bedmethyl file
 *
 *  @param pileup a pileup counts structure.
 *  @param ref reference sequence.
 *  @param rstart starting reference coordinate corresponding to ref.
 *  @param extended whether to include counts of canonical, modified and filtered bases.
 *  @returns void
 *
 */
void print_bedmethyl(plp_data pileup, char *ref, int rstart, bool extended){
    // ecoli1  100718  100719  .       4       +       100718  100719  0,0,0   3       0
    
    // this is a bit naff, we should introspect these indices, or have them
    // as data in the header.
    const size_t numbases = 7;
    size_t fwdbases[] = {4,5,6,7,9,11,13};
    size_t revbases[] = {0,1,2,3,8,10,12};
    size_t ci, mi, fi;
    size_t *bases;
    bool isrev;
    for (size_t i = 0; i < pileup->n_cols; ++i) {
        size_t pos = pileup->major[i];
        // check if this is a position we care about, this could be better by passing
        // in list of relevant contexts:
        //     eg. cCagg, ccaGg, cCtgg, cCtgg
        // where the upper case base corresponds to pos
        char rbase = ref[pos - rstart];
        if (rbase == 'C'){
            isrev = 0; mi = fwd_mod; fi = fwd_filt; ci = 5;
            bases = fwdbases;
        } else if (rbase == 'G') {
            // G on rev strand is C in reads
            isrev = 1; mi = rev_mod; fi = rev_filt; ci = 2;
            bases = revbases;
        }
        else {
            continue;
        }
        // calculate depth on strand
        size_t depth = 0;
        for (size_t j = 0; j < numbases; ++j) {
            depth += pileup->matrix[i * featlen + bases[j]];
        }
        // https://www.encodeproject.org/data-standards/wgbs/
        // column 11: "Percentage of reads that show methylation at this position in the genome"
        //  - Seems to disregard possibility of non-C canonical calls
        // lets calculate this as proportion of meth:non-meth C
        size_t cd = pileup->matrix[i * featlen + ci];
        size_t md = pileup->matrix[i * featlen + mi];
        size_t fd = pileup->matrix[i * featlen + fi];
        size_t tot = cd + md;
        float meth = tot == 0 ? 0 : (100.0f * md) / tot;
        // column 5: "Score from 0-1000. Capped number of reads"
        // lets go with proportion of (mod or canon):(mod or canon or filtered)
        size_t score = depth == 0 ? 0 : (1000 * tot) / depth;

        fprintf(stdout,
            "%s\t%zu\t%zu\t"
            "5mC\t%zu\t%c\t"
            "%zu\t%zu\t0,0,0\t%zu\t%.2f",
            pileup->rname, pos, pos + 1,
            score, "+-"[isrev],
            pos, pos + 1, depth, meth
        );
        if (extended) {
            fprintf(stdout, "\t%zu\t%zu\t%zu\n", cd, md, fd);
        } else {
            fprintf(stdout, "\n");
        }
    }
}


// Control client data for pileup: in this case the mod base data
int pileup_cd_create(void *data, const bam1_t *b, bam_pileup_cd *cd) {
    hts_base_mod_state *m = hts_base_mod_state_alloc();
    bam_parse_basemod(b, m); cd->p = m;
    return 0;
}

int pileup_cd_destroy(void *data, const bam1_t *b, bam_pileup_cd *cd) {
    hts_base_mod_state_free(cd->p);
    return 0;
}


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
        const char *read_group, int lowthreshold, int highthreshold) {

    // open bam etc.
    htsFile *fp = hts_open(bam_file, "rb");
    hts_idx_t *idx = sam_index_load(fp, bam_file);
    sam_hdr_t *hdr = sam_hdr_read(fp);
    if (hdr == 0 || idx == 0 || fp == 0) {
        hts_close(fp); hts_idx_destroy(idx); sam_hdr_destroy(hdr);
        fprintf(stderr, "Failed to read .bam file '%s'.", bam_file);
        return NULL;
    }

    // find the target index for query below
    int mytid = -1;
    for (int i=0; i < hdr->n_targets; ++i) {
        if(!strcmp(hdr->target_name[i], chr)) {
            mytid = i;
            break;
        }
    }
    if (mytid == -1) {
        hts_close(fp); hts_idx_destroy(idx); sam_hdr_destroy(hdr);
        fprintf(stderr, "Failed to find reference sequence '%s' in bam '%s'.", chr, bam_file);
        return NULL;

    }

    // setup bam interator
    mplp_data *data = xalloc(1, sizeof(mplp_data), "pileup init data");
    data->fp = fp; data->hdr = hdr;
    data->iter = bam_itr_queryi(idx, mytid, start, end);
    data->min_mapQ = 1; data->read_group = read_group;

    bam_mplp_t mplp = bam_mplp_init(1, read_bam, (void **)& data);
    const bam_pileup1_t **plp = xalloc(1, sizeof(bam_pileup1_t *), "pileup");
    int ret, pos, tid, n_plp;

    bam_mplp_constructor(mplp, pileup_cd_create);
    bam_mplp_destructor(mplp, pileup_cd_destroy);

    // allocate output, not doing insertions here, so know maximum width
    plp_data pileup = create_plp_data(end - start, chr);

    // get counts
    int n_cols = 0;  // number of processed columns (not all ref positions included)
    size_t major_col = 0;
    int last_pos = start - 1;
    while ((ret=bam_mplp_auto(mplp, &tid, &pos, &n_plp, plp) > 0)) {
        const char *c_name = data->hdr->target_name[tid];
        if (strcmp(c_name, chr) != 0) continue;
        if (pos < start) continue;
        if (pos >= end) break;

        // pileup can skip positions with no coverage
        if (pos - last_pos > 1) fprintf(stderr, "WARNING: no reads span: %s:%d-%d\n", chr, last_pos+1, pos);
        
        n_cols++;
        last_pos = pos;
        pileup->major[n_cols] = pos;  // dont need insert columns for this

        // loop through all reads at this position
        for (int i = 0; i < n_plp; ++i) {
            const bam_pileup1_t *p = plp[0] + i;
            if (p->is_refskip) continue;

            int base_i;
            if (p->is_del) {
                // deletions are interesting for counting depth
                base_i = bam_is_rev(p->b) ? rev_del : fwd_del;
            } else { // handle pos and any following ins
                int base_j = bam1_seqi(bam1_seq(p->b), p->qpos);

                // Simple mod detection
                size_t n_mods = 256;
                hts_base_mod_state *mod_state = p->cd.p;
                hts_base_mod mod[n_mods];
                int nm = bam_mods_at_qpos(p->b, p->qpos, mod_state, mod, n_mods);
                if (nm < 0 ) continue;  // ignore reads which give error
                if (nm > 0 && mod[0].modified_base == 'm') {
                    // just assume mC for now
                    //fprintf(stderr, "Modified %c to %c at %d\n",
                    //    mod[0].canonical_base, mod[0].modified_base, pos);
                    // make decision between mod and unmod.
                    //float q = -10 * log10(1 - ((mod[0].qual + 0.5) / 256)) + 0.5;
                    if (mod[0].qual > highthreshold) {
                        base_i = bam_is_rev(p->b) ? rev_mod : fwd_mod;
                    } else if (mod[0].qual < lowthreshold) {
                        // canonical
                        if bam_is_rev(p->b) base_j += 16;
                        base_i = num2countbase[base_j];
                    } else {
                        // filter out
                        base_i = bam_is_rev(p->b) ? rev_filt : fwd_filt;
                    }
                } else {
                    // no mod call - assume this means canonical
                    if bam_is_rev(p->b) base_j += 16;
                    base_i = num2countbase[base_j];
                }
            }
            if (base_i != -1) {  // not an ambiguity code
                pileup->matrix[major_col + base_i] += 1;
            }
        }
        major_col += featlen;
    }
    if ((end - 1) - last_pos > 1) fprintf(stderr, "WARNING: no reads span: %s:%d-%d\n", chr, last_pos, end);
    pileup->n_cols = n_cols;

    bam_itr_destroy(data->iter);
    bam_mplp_destroy(mplp);
    free(data);
    free(plp);
    hts_close(fp);
    hts_idx_destroy(idx);
    sam_hdr_destroy(hdr);

    return pileup;
}


// Process a region
void process_region(arguments_t args, const char *chr, int start, int end, char *ref) {
    fprintf(stderr, "Processing: %s:%d-%d\n", chr, start, end);
    plp_data pileup = calculate_pileup(
        args.bam, chr, start, end,
        args.read_group, args.lowthreshold, args.highthreshold);
    if (pileup == NULL) return;
    print_bedmethyl(pileup, ref, start, args.extended);
    destroy_plp_data(pileup);
}


// Run program
int main(int argc, char *argv[]) {
    arguments_t args = parse_arguments(argc, argv);
    if (args.highthreshold < args.lowthreshold) {
        fprintf(stderr, "ERROR: --highthreshold must be larger than --lowthreshold\n");
        exit(1);
    }

    // load ref sequence
    faidx_t *fai = fai_load(args.ref);
    if (fai == NULL) exit(1);
    if (args.region == NULL) {
        // process all regions
        int nseq = faidx_nseq(fai);
        for (int i = 0; i < nseq; ++i) {
            const char *chr = faidx_iseq(fai, i);
            int len = faidx_seq_len(fai, chr);
            int alen;
            char *ref = faidx_fetch_seq(fai, chr, 0, len, &alen);
            fprintf(stderr, "Fetched %s, %i %i\n", chr, len, alen);
            process_region(args, chr, 0, len, ref);
            free((void*) chr);
            free(ref);
        }
        exit(0);
    } else {
        // process given region
        int len;
        char *ref = fai_fetch(fai, args.region, &len);
        if (len < 0) {
            exit(1);
        }
        int start, end;
        char *chr = xalloc(strlen(args.region) + 1, sizeof(char), "chr");
        strcpy(chr, args.region);
        char *reg_chr = (char *) hts_parse_reg(chr, &start, &end);
        // start and end now zero-based end exclusive
        if (reg_chr) {
            *reg_chr = '\0';  // sets chr to be terminated at correct point
        } else {
            fprintf(stderr, "Failed to parse region: '%s'.\n", args.region);
        }
        process_region(args, chr, start, end, ref);

        free(chr);
        free(ref);
    }
    exit(0);
}
