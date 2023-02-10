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
#include <unistd.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/thread_pool.h"

#include "bamiter.h"
#include "common.h"
#include "counts.h"

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


output_files open_bed_files(char* prefix, bool cpg, bool chh, bool chg, bool accumulated) {
    output_files files = xalloc(1, sizeof(_output_files), "output_files");
    // default to stdout for zero or one filters
    files->multi = (int)cpg + chh + chg > 1;
    files->take_all = (int)cpg + chh + chg == 0;
    files->accumulated = accumulated;
    files->fcpg = stdout;
    files->fchh = stdout;
    files->fchg = stdout;
    files->fcpg_acc = NULL;
    files->fchh_acc = NULL;
    files->fchg_acc = NULL;
    files->cpg = cpg;
    files->chh = chh;
    files->chg = chg;
    // use distinct files if more than one filter
    if (files->multi) {
        char* fname = xalloc(strlen(prefix) + 9, sizeof(char), "fname");
        if (cpg) {
            strcpy(fname, prefix); strcat(fname, ".cpg.bed");
            files->fcpg = fopen(fname, "w");
        }
        if (chh) {
            strcpy(fname, prefix); strcat(fname, ".chh.bed");
            files->fchh = fopen(fname, "w");
        }
        if (chg) {
            strcpy(fname, prefix); strcat(fname, ".chg.bed");
            files->fchg = fopen(fname, "w");
        }
        free(fname);
    }

    if (files->accumulated) {
        char* fname_acc = xalloc(strlen(prefix) + 13, sizeof(char), "fname");
        if (cpg) {
            strcpy(fname_acc, prefix); strcat(fname_acc, ".cpg.acc.bed");
            files->fcpg_acc = fopen(fname_acc, "w");
        }
        if (chg) {
            strcpy(fname_acc, prefix); strcat(fname_acc, ".chg.acc.bed");
            files->fchg_acc = fopen(fname_acc, "w");
        }
        free(fname_acc);
    }

    // store these in an array for later ease
    // [CpG, CHG]
    init_output_buffers(files);
    files->buf_size = _buf_size;
    files->motif_offsets[0] = 1;
    files->motif_offsets[1] = 2;
    files->motif_acc_files[0] = files->fcpg_acc;
    files->motif_acc_files[1] = files->fchg_acc;
    return files;
}

void close_bed_files(output_files files) {
   if (files->fcpg != stdout) { fclose(files->fcpg); }
   if (files->fchh != stdout) { fclose(files->fchh); }
   if (files->fchg != stdout) { fclose(files->fchg); }
   if (files->fcpg_acc != NULL) { fclose(files->fcpg_acc); }
   if (files->fchh_acc != NULL) { fclose(files->fchh_acc); }
   if (files->fchg_acc != NULL) { fclose(files->fchg_acc); }
   free(files);
}


// Check sequences for motifs

// CpG
bool extern inline is_cpg_fwd(size_t rpos, int rlen, char* ref){
    return rpos < rlen - 1 && ref[rpos] == 'C' && ref[rpos + 1] == 'G';
}
bool extern inline is_cpg_rev(size_t rpos, int rlen, char* ref){
    return rpos != 0 && ref[rpos] == 'G' && ref[rpos - 1] == 'C';
}
// CHN
bool extern inline _is_chn_fwd(size_t rpos, int rlen, char* ref) {
    bool is_chn = false;
    if (rpos < rlen - 2 && ref[rpos] == 'C') {
        char b = ref[rpos + 1];
        // these are all not G
        is_chn = (b == 'A' || b == 'C' || b == 'T' || b == 'M' || b == 'W' || b == 'Y' || b == 'H');
    }
    return is_chn;
}
bool extern inline _is_chn_rev(size_t rpos, int rlen, char* ref) {
    bool is_chn = false;
    if (rpos > 1 && ref[rpos] == 'G') {
        char b = ref[rpos - 1];
        // these are all not C
        is_chn = (b == 'A' || b == 'G' || b == 'T' || b == 'R' || b == 'W' || b == 'K' || b == 'D');
    }
    return is_chn;
}
// CHH
bool extern inline is_chh_fwd(size_t rpos, int rlen, char* ref) {
    bool is_chh = _is_chn_fwd(rpos, rlen, ref);
    if (is_chh) { 
        char b = ref[rpos + 2];
        // these are all not G
        is_chh = (b == 'A' || b == 'C' || b == 'T' || b == 'M' || b == 'W' || b == 'Y' || b == 'H');
    }
    return is_chh;
}
bool extern inline is_chh_rev(size_t rpos, int rlen, char* ref) {
    bool is_chh = _is_chn_rev(rpos, rlen, ref);
    if (is_chh) {
        char b = ref[rpos - 2];
        // these are all not C
        is_chh = (b == 'A' || b == 'G' || b == 'T' || b == 'R' || b == 'W' || b == 'K' || b == 'D');
    }
    return is_chh;
}
// CHG
bool extern inline is_chg_fwd(size_t rpos, int rlen, char* ref) {
    bool is_chg = _is_chn_fwd(rpos, rlen, ref);
    if (is_chg) {
        is_chg = ref[rpos + 2] == 'G';
    }
    return is_chg;
}
bool extern inline is_chg_rev(size_t rpos, int rlen, char* ref) {
    bool is_chg = _is_chn_rev(rpos, rlen, ref);
    if (is_chg) {
        is_chg = ref[rpos - 2] == 'C';
    }
    return is_chg;
}


void inline print_record(
        FILE* fout, const char* rname, size_t start, size_t end,
        char* feature, char orient, size_t depth,
        bool extended, size_t cd, size_t md, size_t fd, size_t xd, size_t od) {
    // https://www.encodeproject.org/data-standards/wgbs/
    // column 11: "Percentage of reads that show methylation at this position in the genome"
    //  - Seems to disregard possibility of non-C canonical calls
    // lets calculate this as proportion of meth:non-meth C
    size_t tot = cd + md + od;
    float meth = tot == 0 ? nanf("") : (100.0f * md) / tot;
    // column 5: "Score from 0-1000. Capped number of reads"
    // lets go with proportion of (mod or canon):(mod or canon or filtered)
    size_t score = depth == 0 ? nanf("") : (1000 * tot) / depth;

    // TODO: don't print when nan?
    fprintf(fout,
        "%s\t%zu\t%zu\t"
        "%s\t%zu\t%c\t"
        "%zu\t%zu\t0,0,0\t%zu\t%.2f",
        rname, start, end,
        feature, score, orient,
        start, end, depth, meth);
    if (extended) {
        fprintf(fout, "\t%zu\t%zu\t%zu\t%zu\t%zu\n", cd, md, fd, xd, od);
    } else {
        fprintf(fout, "\n");
    }
}


void init_output_buffers(output_files bed_files) {
    // information regarding motif offset pairing
    for (size_t i=0; i < bed_files->buf_size; ++i) {
        bed_files->out_buffer[i] = (bed_buffer){-1, false, 0, 0, 0, 0};
    }
}

void flush_output_buffers(output_files bed_files, const char* chr, bool extended, char* feature) {
    // flush accumulation buffers
    if (bed_files->accumulated) {
        for(size_t ibuf=0; ibuf < bed_files->buf_size; ++ibuf) {
            bed_buffer buf = bed_files->out_buffer[ibuf];
            FILE* fout = bed_files->motif_acc_files[ibuf];
            if (buf.pos != -1 && fout != NULL) {
                print_record(
                    fout, chr, buf.pos, buf.pos + 1, feature, "+-"[buf.isrev],
                    buf.depth, extended, buf.cd, buf.md, buf.fd, buf.xd, buf.od);
            }
        }
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
 *  @param output_files file handles and output options.
 *  @param out_buffer state for strand accumulation (modified on output).
 *  @returns void
 *
 */
void print_bedmethyl(
        plp_data pileup, char *ref, int rstart, bool extended,
        char* feature, char canon_base, output_files bed_files) {
    // ecoli1  100718  100719  .       4       +       100718  100719  0,0,0   3       0
    
    // this is a bit naff, we should introspect these indices, or have them
    // as data in the header.
    size_t ci, mi, fi, xi, oi;
    size_t *bases;
    bool isrev;
    char rc_canon_base = ' ';
    size_t cif, cir;

    // TODO: if canon_base were passed as an htslib int this would be cleaner
    if      (canon_base == 'A') {cif=fwd_A; cir=rev_T; rc_canon_base = 'T';}
    else if (canon_base == 'C') {cif=fwd_C; cir=rev_G; rc_canon_base = 'G';}
    else if (canon_base == 'G') {cif=fwd_G; cir=rev_C; rc_canon_base = 'C';}
    else if (canon_base == 'T') {cif=fwd_T; cir=rev_A; rc_canon_base = 'A';}
    else {fprintf(stderr, "ERROR: Unrecognised canonical base: '%c'\n", canon_base); exit(1);}

    int rlen = strlen(ref);

    for (size_t i = 0; i < pileup->n_cols; ++i) {
        size_t pos = pileup->major[i];
        size_t rpos = pos - rstart;
        char rbase = ref[rpos];
        bool is_cpg = false;
        bool is_chh = false;
        bool is_chg = false;
        if (rbase == canon_base) {
            if (!bed_files->take_all) {
                if (!(
                       (bed_files->cpg && (is_cpg = is_cpg_fwd(rpos, rlen, ref)))
                    || (bed_files->chh && (is_chh = is_chh_fwd(rpos, rlen, ref)))
                    || (bed_files->chg && (is_chg = is_chg_fwd(rpos, rlen, ref)))
                    ) ) { continue; }
            }
            isrev = 0; mi = fwd_mod; fi = fwd_filt; xi = fwd_nocall; oi = fwd_other; ci = cif;
            bases = (size_t *) fwdbases;
        } else if (rbase == rc_canon_base) {
            if (!bed_files->take_all) {
                if (!(
                       (bed_files->cpg && (is_cpg = is_cpg_rev(rpos, rlen, ref)))
                    || (bed_files->chh && (is_chh = is_chh_rev(rpos, rlen, ref)))
                    || (bed_files->chg && (is_chg = is_chg_rev(rpos, rlen, ref)))
                    ) ) { continue; }
            }
            isrev = 1; mi = rev_mod; fi = rev_filt; xi = rev_nocall; oi = rev_other; ci = cir;
            bases = (size_t *)revbases;
        }
        else {
            continue;
        }
        // calculate depth on strand
        size_t depth = 0;
        for (size_t j = 0; j < numbases; ++j) {
            depth += pileup->matrix[i * featlen + bases[j]];
        }
        size_t cd = pileup->matrix[i * featlen + ci];
        size_t md = pileup->matrix[i * featlen + mi];
        size_t fd = pileup->matrix[i * featlen + fi];
        size_t xd = pileup->matrix[i * featlen + xi];
        size_t od = pileup->matrix[i * featlen + oi];

        // choose output for this locus, the motifs are mutually exclusive so
        // no need to loop
        FILE* fout = stdout;
        if (bed_files->multi) {
            if (is_cpg) { fout = bed_files->fcpg; }
            else if (is_chh) { fout = bed_files->fchh; }
            else if (is_chg) { fout = bed_files->fchg; }
        }
        print_record(
            fout, pileup->rname, pos, pos + 1, feature, "+-"[isrev],
            depth, extended, cd, md, fd, xd, od);

        // strand accumulated
        if (bed_files->accumulated && (is_cpg || is_chg)) {
            size_t ibuf, motif_offset;
            bool do_output;
            if (is_cpg) {
                ibuf = 0; do_output = bed_files->cpg;
            } else { // chg
                ibuf = 1; do_output = bed_files->chh;
            }
            motif_offset = bed_files->motif_offsets[ibuf];
            fout = bed_files->motif_acc_files[ibuf];
            if (do_output) {
                assert(fout != NULL);
                bed_buffer buf = bed_files->out_buffer[ibuf];
                if (buf.pos == -1) {
                    bed_files->out_buffer[ibuf] = (bed_buffer){pos, isrev, depth, cd, md, fd, xd, od};
                } else if (pos - buf.pos == motif_offset ) { // paired
                    assert(buf.isrev != isrev); // shouldn't happen, they can't be same
                    buf.depth += depth;
                    buf.cd += cd;
                    buf.md += md;
                    buf.fd += fd;
                    buf.xd += xd;
                    buf.od += od;
                    print_record(
                        fout, pileup->rname, buf.pos, buf.pos + motif_offset + 1, feature, '.',
                        buf.depth, extended, buf.cd, buf.md, buf.fd, buf.xd, buf.od);
                    bed_files->out_buffer[ibuf] = (bed_buffer){-1, false, 0, 0, 0, 0, 0, 0};
                } else { // unrelated
                    print_record(
                        fout, pileup->rname, buf.pos, buf.pos + 1, feature, "+-"[buf.isrev],
                        buf.depth, extended, buf.cd, buf.md, buf.fd, buf.xd, buf.od);
                    bed_files->out_buffer[ibuf] = (bed_buffer){pos, isrev, depth, cd, md, fd, xd, od};
                }
            }
        }

    } // position loop

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


// TODO: this is taken from sam.c, its here so we can introspec some things
//       for which there's no public interface. A little spicey to redefine
//       this, but we do what we can.
//       https://github.com/samtools/htslib/issues/1550
#define MAX_BASE_MOD 256
struct hts_base_mod_state {
    int type[MAX_BASE_MOD];     // char or minus-CHEBI
    int canonical[MAX_BASE_MOD];// canonical base, as seqi (1,2,4,8,15)
    char strand[MAX_BASE_MOD];  // strand of modification; + or -
    int MMcount[MAX_BASE_MOD];  // no. canonical bases left until next mod
    char *MM[MAX_BASE_MOD];     // next pos delta (string)
    char *MMend[MAX_BASE_MOD];  // end of pos-delta string
    uint8_t *ML[MAX_BASE_MOD];  // next qual
    int MLstride[MAX_BASE_MOD]; // bytes between quals for this type
    int implicit[MAX_BASE_MOD]; // treat unlisted positions as non-modified?
    int seq_pos;                // current position along sequence
    int nmods;                  // used array size (0 to MAX_BASE_MOD-1).
};


// Query if a specific MM subtag is present
bool query_mod_subtag(hts_base_mod_state *state, int qtype, int qcanonical, char qstrand, int qimplicit) {
    bool found = false;
    for (size_t i=0; i<state->nmods; ++i) {
        if ((state->type[i] == qtype || state->type[i] == -qtype)
                && state->canonical[i] == qcanonical
                // although strand is typed char and documented as + or -, its actually 0/1
                && "+-"[state->strand[i]] == qstrand
                && state->implicit[i] == qimplicit) {
            found = true;
            break;
        }
    }
    return found;
}

/** Generates base counts from a region of a bam.
 *
 *  @param bam_file input aligment file.
 *  @param chr bam target name.
 *  @param start start position of chr to consider.
 *  @param end end position of chr to consider.
 *  @param read_group by which to filter alignments.
 *  @param tag_name by which to filter alignments.
 *  @param tag_value associated with tag_name
 *  @param threshold probability filter for excluding calls from counts.
 *  @param mb BAM code for modified base to report. (e.g. h for 5hmC), or a ChEBI code.
 *  @param combine combine all modified bases corresponding to same canonical base as mb
 *  @returns a pileup data pointer.
 *
 *  The return value can be freed with destroy_plp_data.
 *
 */
plp_data calculate_pileup(
        const set_fsets *fsets, const char *chr, int start, int end,
        const char *read_group, const char tag_name[2], const int tag_value,
        int threshold, mod_base mb, bool combine, int max_depth) {

    static bool shown_second_strand_warning = false;

    // counting mod calls other than the one asked for
    int rev_in_family = rev_other;
    int fwd_in_family = fwd_other;
    if (combine) { rev_in_family = rev_mod; fwd_in_family = fwd_mod; }

    // setup bam reading
    size_t nfile = fsets->n;
    mplp_data **data = xalloc(fsets->n, sizeof(mplp_data*), "bam files");
    for (size_t i = 0; i < nfile; ++i) {
        data[i] = create_bam_iter_data(
            fsets->fsets[i], chr, start, end, read_group, tag_name, tag_value);
        if (data[i] == NULL) {
            // TODO: clean-up all j<i data[i], and free data
            return NULL;
        }
    }

    bam_mplp_t mplp = bam_mplp_init(nfile, read_bam, (void **)data);
    int *n_plp = xalloc(nfile, sizeof(int), "bam read cover");
    const bam_pileup1_t **plp = xalloc(nfile, sizeof(bam_pileup1_t *), "pileup");
    int ret, pos, tid;

    bam_mplp_constructor(mplp, pileup_cd_create);
    bam_mplp_destructor(mplp, pileup_cd_destroy);
    bam_mplp_set_maxcnt(mplp, max_depth);

    // allocate output, not doing insertions here, so know maximum width
    plp_data pileup = create_plp_data(end - start, chr);

    // get counts
    int n_cols = 0;  // number of processed columns (not all ref positions included)
    size_t major_col = 0;
    while ((ret=bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0)) {
        const char *c_name = data[0]->hdr->target_name[tid];
        if (strcmp(c_name, chr) != 0) continue;
        if (pos < start) continue;
        if (pos >= end) break;

        pileup->major[n_cols] = pos;  // dont need insert columns for this

        // go through all files, and all reads in each
        for (size_t file = 0; file < nfile; ++file) {
            for (int i = 0; i < n_plp[file]; ++i) {
                const bam_pileup1_t *p = plp[file] + i;
                if (p->is_refskip) continue;

                // ONT calls are "query based", this means an attempt at a mod call is
                // made only if the first-pass canon basecall was the base of interest.
                // They are NOT "reference based": a mod call being attempted when the
                // query position aligns to a reference position containing the
                // of-interest base. (Actually reading between the lines of the spec
                // discussions, there was an implied assumption that mod calls are
                // always query based).
                //
                // There are two modes:
                //  i) "." - implicit = 1; Unlisted positions are assumed canonical
                // ii) "?" - implicit = 0; Nothing can really be said about unlisted
                //
                // Case i) is trivial and easy to handle: no mod calls, assume canonical.
                // This is like just not having a tag at all. If the above found no mods,
                // any query base (ACGT) is assumed canonical
                //
                // Case ii) is a bit more icky for us. Before deciding canon/no-call we
                // need to know if there was even a tag present, e.g. C+m for 5mC. For
                // canon base types other than that relating to our mod base, we make
                // no claims about modification status: all forms are lumped together.
                //
                // For the most part ONT callers output `?` and have a call for every
                // of-interest base. There are two cases where this isn't true:
                //  i) Guppy elided some low prob calls (as in the `.` mode)
                // ii) callers which specialise to CpG (so don't have an entry for every C)
                //
                // To complicate things further we can have tags such as "G-m" indicating
                // methylation on the second strand of the sequenced read. Such tags ought
                // not to occur without a corresponding "C+m" tag: in a simple case this
                // would imply a caller had called methylation on the strand that wasn't
                // sequenced but not on the strand that was sequenced. A more realistic
                // situation would be making calls only on the second strands of duplex reads.
                // 
                // Here we simplify our lives by restricting to the case of skipping any
                // such second strand tags, for the reasons above but also primarily
                // because ideally the second strand tag should be jointly interpreted
                // with the first strand tag:
                //    to detect hemimethylation
                //    understand and correctly report depth
                //    made hard by them being on different positions

                int base_i = -1;  // index into counts matrix
                int base_j = bam1_seqi(bam1_seq(p->b), p->qpos);
                if (p->is_del) {
                    // deletions are interesting for counting depth
                    base_i = bam_is_rev(p->b) ? rev_del : fwd_del;
                } else if (!(
                        (base_j == mb.base_i && !bam_is_rev(p->b))
                        || (seqi_rc[base_j] == mb.base_i && bam_is_rev(p->b)))) {
                    // e.g. if query we're looking for 5mC and qbase in {A,T}
                    //      we'll just count a plain A/T
                    // NOTE: this test assumes only first strand subtags (e.g. C+m, not C-m)
                    base_i = num2countbase[bam_is_rev(p->b) ? base_j + 16: base_j];
                } else {		
                    // We have the correct query base for the orientation of the alignment
                    // so now look for modified bases.
                    size_t n_mods = 256;
                    hts_base_mod_state *mod_state = p->cd.p;
                    hts_base_mod allmod[n_mods];
                    int nm = bam_mods_at_qpos(p->b, p->qpos, mod_state, allmod, n_mods);
                    if (nm < 0 ) continue;  // ignore reads which give error
                    hts_base_mod mod;
                    int our_mod = -1;
                    int best_mod = -1;
                    int best_score = 0;
                    int canon_score = MAX_QUAL;  // we subtract from this below
                    if (nm > 0) {
                        for (int k = 0; k < nm && k < n_mods; ++k) {
                            mod = allmod[k];
                            if (mod.strand == 1) {  // second strand tag
                                if (!shown_second_strand_warning) {
                                    fprintf(stderr, "WARNING: Skipping second strand tag.");
                                    shown_second_strand_warning = true;
                                }
                                continue;
                            }
                            // our mod
                            if (mb.code == mod.modified_base || mb.code == -mod.modified_base) {
                                our_mod = k;
                            }
                            // any mod in the family
                            if (mod.canonical_base == mb.base) {
                                if (mod.qual > best_score) { best_mod = k; best_score = mod.qual; }
                                canon_score -= mod.qual;
                            }
                        }
                    }
                    
                    // Now analyse scores. Note: ignoring the old lowthreshold here.
                    if (best_mod != -1) {
                        // we found some mods, lets not worry about funny mixes
                        // of calls and no calls i.e. were assuming we have a call
                        // for all the mods present (implicit non-mod doesn't matter here therefore).
                        if (canon_score > threshold) { // implied canon score 
                            base_i = num2countbase[bam_is_rev(p->b) ? base_j + 16 : base_j];
                        }
                        else if (best_mod == our_mod) { // the mod requested
                            base_i = (best_score > threshold) ?
                                (bam_is_rev(p->b) ? rev_mod : fwd_mod) :
                                (bam_is_rev(p->b) ? rev_filt : fwd_filt);
                        }
                        else { // some other mod in the family
                            base_i = (best_score > threshold) ?
                                (bam_is_rev(p->b) ? rev_in_family : fwd_in_family) :  // either mod or other depending on combine
                                (bam_is_rev(p->b) ? rev_filt : fwd_filt);
                        }
                    }
                    else {
                        // we didn't find any mods in the family
                        // In the case of explicit `?`
                        // tag we should not assume canonical, otherwise we can.
                        // NOTE: we don't look for second strand `-` tags.
                        //       or a mess of `?` and `.` for alternative mods
                        if (query_mod_subtag(mod_state, mb.code, mb.base_i, '+', 0)) {
                            // we had an explicit tag, but no call for this position
                            base_i = bam_is_rev(p->b) ? rev_nocall : fwd_nocall;
                        }
                        else {
                            // for everything else theres canonical
                            base_i = num2countbase[bam_is_rev(p->b) ? base_j + 16 : base_j];
                        }
                    }
                }
                if (base_i != -1) {  // not an ambiguity code
                    pileup->matrix[major_col + base_i] += 1;
                } // read loop
            } // file loop
        }
        major_col += featlen;
        n_cols++;
    }
    pileup->n_cols = n_cols;

    free(plp);
    free(n_plp);
    bam_mplp_destroy(mplp);
    for (size_t i = 0; i < nfile; ++i) {
        destroy_bam_iter_data(data[i]);
    }
    free(data);

    return pileup;
}


