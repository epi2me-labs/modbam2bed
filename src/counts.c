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
#include <pthread.h>
#include <unistd.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
//#include "threadpool.h"
#include "htslib/thread_pool.h"

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
    //fprintf(stderr, buffer_cols); 
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
 *  @param feature name to use for feature column of BED (e.g. 5mC).
 *  @param canon_base canonical base to match.
 *  @param cpg filter output to only CpG sites.
 *  @returns void
 *
 */
void print_bedmethyl(plp_data pileup, char *ref, int rstart, bool extended, char* feature, char canon_base, bool cpg){
    // ecoli1  100718  100719  .       4       +       100718  100719  0,0,0   3       0
    
    // this is a bit naff, we should introspect these indices, or have them
    // as data in the header.
    const size_t numbases = 7;
    size_t fwdbases[] = {4,5,6,7,9,11,13};
    size_t revbases[] = {0,1,2,3,8,10,12};
    size_t ci, mi, fi;
    size_t *bases;
    bool isrev;
    char rc_canon_base = ' ';
    size_t cif, cir;
    if (canon_base == 'A') {cif=4; cir=3; rc_canon_base = 'T';}
    else if (canon_base == 'C') {cif=5; cir=2; rc_canon_base = 'G';}
    else if (canon_base == 'G') {cif=6; cir=1; rc_canon_base = 'C';}
    else if (canon_base == 'T') {cif=7; cir=0; rc_canon_base = 'A';}
    else {fprintf(stderr, "Unrecognised canonical base: '%c'\n", canon_base); exit(1);}
    if (canon_base != 'C' && cpg) {
        fprintf(stderr, "CpG filtering cannot be used when canonical base is not 'C'.\n");
        exit(1);
    }

    int rlen = strlen(ref);

    for (size_t i = 0; i < pileup->n_cols; ++i) {
        size_t pos = pileup->major[i];
        size_t rpos = pos - rstart;
        char rbase = ref[rpos];
        if (rbase == canon_base){
            if (cpg && rpos < rlen - 1  && ref[rpos + 1] != 'G') {
                continue;
            }
            isrev = 0; mi = fwd_mod; fi = fwd_filt; ci = cif;
            bases = fwdbases;
        } else if (rbase == rc_canon_base) {
            // e.g. G on rev strand is C in reads
            if (cpg && rpos != 0 && ref[rpos - 1] != 'C') {
                continue;
            }
            isrev = 1; mi = rev_mod; fi = rev_filt; ci = cir;
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
        if (depth == 0) continue;
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
            "%s\t%zu\t%c\t"
            "%zu\t%zu\t0,0,0\t%zu\t%.2f",
            pileup->rname, pos, pos + 1,
            feature, score, "+-"[isrev],
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
 *  @param mod_base BAM code for modified base to report. (e.g. h for 5hmC).
 *  @returns a pileup data pointer.
 *
 *  The return value can be freed with destroy_plp_data.
 *
 */
plp_data calculate_pileup(
        const char *bam_file, const char *chr, int start, int end,
        const char *read_group, int lowthreshold, int highthreshold, char mod_base) {

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
    while ((ret=bam_mplp_auto(mplp, &tid, &pos, &n_plp, plp) > 0)) {
        const char *c_name = data->hdr->target_name[tid];
        if (strcmp(c_name, chr) != 0) continue;
        if (pos < start) continue;
        if (pos >= end) break;

        pileup->major[n_cols] = pos;  // dont need insert columns for this

        // loop through all reads at this position
        for (int i = 0; i < n_plp; ++i) {
            const bam_pileup1_t *p = plp[0] + i;
            if (p->is_refskip) continue;

            int base_i = -1;
            if (p->is_del) {
                // deletions are interesting for counting depth
                base_i = bam_is_rev(p->b) ? rev_del : fwd_del;
            } else { // handle pos and any following ins
                int base_j = bam1_seqi(bam1_seq(p->b), p->qpos);

                // Simple mod detection
                size_t n_mods = 256;
                hts_base_mod_state *mod_state = p->cd.p;
                hts_base_mod allmod[n_mods];
                int nm = bam_mods_at_qpos(p->b, p->qpos, mod_state, allmod, n_mods);
                if (nm < 0 ) continue;  // ignore reads which give error
                if (nm > 0) {
                    hts_base_mod mod;
                    for (int k = 0; k < nm && k < n_mods; ++k) {
                        if (allmod[k].modified_base == mod_base) {
                            mod = allmod[k];
                            // we found our mod
                            //fprintf(stderr, "Modified %c to %c at %d\n",
                            //    mod.canonical_base, mod.modified_base, pos);
                            // make decision between mod and unmod.
                            //float q = -10 * log10(1 - ((mod[0].qual + 0.5) / 256)) + 0.5;
                            if (mod.qual > highthreshold) {
                                base_i = bam_is_rev(p->b) ? rev_mod : fwd_mod;
                            } else if (mod.qual < lowthreshold) {
                                // canonical
                                if bam_is_rev(p->b) base_j += 16;
                                base_i = num2countbase[base_j];
                            } else {
                                // filter out
                                base_i = bam_is_rev(p->b) ? rev_filt : fwd_filt;
                            }
                            break;
                        }
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
        n_cols++;
    }
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


typedef struct twarg {
    arguments_t args;
    const char *chr;
    int start;
    int end;
} twarg;

void *pileup_worker(void *arg) {
    twarg j = *(twarg *)arg;
    plp_data pileup = calculate_pileup(
        j.args.bam, j.chr, j.start, j.end,
        j.args.read_group, j.args.lowthreshold, j.args.highthreshold, j.args.mod_base.code);
    free(arg);
    return pileup;
}

// Process and print a single region using a threadpool
void process_region_threads(arguments_t args, const char *chr, int start, int end, char *ref) {
    // create thread pool
    hts_tpool *p = hts_tpool_init(args.threads);
    hts_tpool_process *q = hts_tpool_process_init(p, 2 * args.threads, 0);
    hts_tpool_result *r;
    const int width = 1000000;

    fprintf(stderr, "Processing: %s:%d-%d\n", chr, start, end);
    int nregs = 1 + (end - start) / width; float done = 0;
    for (int rstart = start; rstart < end; rstart += width) {
        twarg *tw_args = xalloc(1, sizeof(*tw_args), "thread worker args");  // freed in worker
        tw_args->args = args;
        tw_args->chr = chr; tw_args->start = rstart; tw_args->end=min(rstart + width, end);
        int blk;
        do {
            blk = hts_tpool_dispatch2(p, q, pileup_worker, tw_args, 1);
            if ((r = hts_tpool_next_result(q))) {
                plp_data res = (plp_data)hts_tpool_result_data(r);
                if (res != NULL) {
                    print_bedmethyl(
                        res, ref, 0,
                        args.extended, args.mod_base.abbrev, args.mod_base.base, args.cpg);
                    destroy_plp_data(res);
                    done++;
                    fprintf(stderr, "\r%.1f %%", 100*done/nregs);
                }
                hts_tpool_delete_result(r, 0);
            }
        } while (blk == -1);
    }

    // wait for jobs, then collect.
    hts_tpool_process_flush(q);
    while ((r = hts_tpool_next_result(q))) {
        plp_data res = (plp_data)hts_tpool_result_data(r);
        if (res != NULL) {
            print_bedmethyl(
                res, ref, 0,
                args.extended, args.mod_base.abbrev, args.mod_base.base, args.cpg);
            destroy_plp_data(res);
            done++;
            fprintf(stderr, "\r%.1f %%", 100*done/nregs);
        }
        hts_tpool_delete_result(r, 0);
    }
    fprintf(stderr, "\r100 %%  ");
    fprintf(stderr, "\n");
    // clean up pool
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
}


// Process and print a single region
void process_region(arguments_t args, const char *chr, int start, int end, char *ref) {
    plp_data pileup = calculate_pileup(
        args.bam, chr, start, end,
        args.read_group, args.lowthreshold, args.highthreshold, args.mod_base.code);
    if (pileup == NULL) return;
    print_bedmethyl(pileup, ref, start, args.extended, args.mod_base.abbrev, args.mod_base.base, args.cpg);
    destroy_plp_data(pileup);
}


// Run program
int main(int argc, char *argv[]) {
    arguments_t args = parse_arguments(argc, argv);
    if (args.highthreshold < args.lowthreshold) {
        fprintf(stderr, "ERROR: --highthreshold must be larger than --lowthreshold\n");
        exit(1);
    }
    fprintf(
        stderr, "Analysing: %s (%s, %c>%c)\n",
        args.mod_base.name, args.mod_base.abbrev, args.mod_base.base, args.mod_base.code);

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
            process_region_threads(args, chr, 0, len, ref);
            free((void*) chr);
            free(ref);
        }
        exit(0);
    } else {
        // process given region
        // simplify things for later (motif matching) on by fetching whole chr
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
        int len;
        char *ref = fai_fetch(fai, chr, &len);
        if (len < 0) {
            exit(1);
        }
        end = min(end, len);
        process_region_threads(args, chr, start, end, ref);

        free(chr);
        free(ref);
    }
    fai_destroy(fai);
    exit(0);
}
