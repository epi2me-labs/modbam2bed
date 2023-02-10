// modbam2bed program

#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <time.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/thread_pool.h"

#include "bamiter.h"
#include "common.h"
#include "counts.h"
#include "args.h"


typedef struct twarg {
    arguments_t args;
    const char *chr;
    int start;
    int end;
} twarg;


void *pileup_worker(void *arg) {
    twarg j = *(twarg *)arg;
    set_fsets *files = create_filesets(j.args.bam);
    if (files == NULL) { free(arg); return NULL; }
    plp_data pileup = calculate_pileup(
        files, j.chr, j.start, j.end,
        j.args.read_group, j.args.tag_name, j.args.tag_value,
        j.args.threshold, j.args.mod_base, j.args.combine,
        j.args.hts_maxcnt);
    destroy_filesets(files);
    free(arg);
    return pileup;
}


/* Process and print a single region using a threadpool
 *
 * @param args program arguments.
 * @param chr reference sequence to process.
 * @param start reference coordinate to process (0-based).
 * @param end reference coordiate to process (exclusive).
 * @param ref reference sequence.
 *
 */
#ifdef NOTHREADS
void process_region(arguments_t args, const char *chr, int start, int end, char *ref, output_files bed_files) {
    fprintf(stderr, "Processing: %s:%d-%d\n", chr, start, end);
    set_fsets* files = create_filesets(j.args.bam);
    if (files == NULL) return;
    plp_data pileup = calculate_pileup(
        args.bam, chr, start, end,
        args.read_group, args.tag_name, args.tag_value,
        args.threshold, args.mod_base, args.combine,
        args.hts_maxcnt);
    if (pileup == NULL) return;

    init_output_buffers(bed_files);
    if (args.pileup) {
        print_pileup_data(pileup);
    } else {
        print_bedmethyl(pileup, ref, 0, args.extended, args.mod_base.abbrev, args.mod_base.base, bed_files);
    }
    flush_output_buffers(bed_files, chr, args.extended, args.mod_base.abbrev);
    destroy_plp_data(pileup);
}
#else
void process_region(arguments_t args, const char *chr, int start, int end, char *ref, output_files bed_files) {
    fprintf(stderr, "Processing: %s:%d-%d\n", chr, start, end);
    // create thread pool
    hts_tpool *p = hts_tpool_init(args.threads);
    hts_tpool_process *q = hts_tpool_process_init(p, 2 * args.threads, 0);
    hts_tpool_result *r;
    const int width = 1000000;

    init_output_buffers(bed_files);
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
                    if (args.pileup) {
                        print_pileup_data(res);
                    } else {
                        print_bedmethyl(
                            res, ref, 0,
                            args.extended, args.mod_base.abbrev, args.mod_base.base, bed_files);
                    }
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
            if (args.pileup) {
                print_pileup_data(res);
            } else {
                print_bedmethyl(
                    res, ref, 0,
                    args.extended, args.mod_base.abbrev, args.mod_base.base, bed_files);
            }
            destroy_plp_data(res);
            done++;
            fprintf(stderr, "\r%.1f %%", 100*done/nregs);
        }
        hts_tpool_delete_result(r, 0);
    }

    // finalise any remaining singleton strands
    flush_output_buffers(bed_files, chr, args.extended, args.mod_base.abbrev);

    fprintf(stderr, "\r100 %%  ");
    fprintf(stderr, "\n");
    // clean up pool
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
}
#endif


int main(int argc, char *argv[]) {
    clock_t begin = clock();
    arguments_t args = parse_arguments(argc, argv);
    fprintf(
        stderr, "Analysing: %s (%s, %c>%c)\n",
        args.mod_base.name, args.mod_base.abbrev, args.mod_base.base, args.mod_base.code);
#ifdef NOTHREADS
    if (args.threads != 1) {
        fprintf(
            stderr,
            "--threads set to %d, but threading not supported by this build.\n", args.threads);
    }
#endif

    // large basecaller runs can produce more files than a single
    // process can open, check this ahead of time.
#ifndef WASM
    struct rlimit reslimit;
    int nfile = 0; for (; args.bam[nfile]; nfile++);
    if (getrlimit(RLIMIT_NOFILE, &reslimit) == 0) {
        if (nfile * args.threads > reslimit.rlim_cur - 100) {
            fprintf(stderr,
                "ERROR: Too many BAM files provided (%i). Try running "
                "samtools merge on subsets of files to produce fewer files", nfile);
            exit(EXIT_FAILURE);
        }
    }
#endif

    // open output files, sort out filter options
    output_files bed_files = open_bed_files(
        args.prefix, args.cpg, args.chh, args.chg, args.accumulated);

    // load ref sequence
    faidx_t *fai = fai_load(args.ref);
    if (fai == NULL) {
        fprintf(stderr,
            "ERROR: Failed to parse reference file\n");
        exit(EXIT_FAILURE);
    }
    if (args.region == NULL) {
        // process all regions
        int nseq = faidx_nseq(fai);
        for (int i = 0; i < nseq; ++i) {
            const char *chr = faidx_iseq(fai, i);
            int len = faidx_seq_len(fai, chr);
            int alen;
            char *ref = faidx_fetch_seq(fai, chr, 0, len, &alen);
            if (!args.mask) {
                for (size_t i=0; i<alen; ++i){ ref[i] = toupper(ref[i]); }
            }
            fprintf(stderr, "Fetched %s, %i %i\n", chr, len, alen);
            process_region(args, chr, 0, len, ref, bed_files);
            free(ref);
        }
    } else {
        // process given region
        int start, end;
        char *chr = xalloc(strlen(args.region) + 1, sizeof(char), "chr");
        strcpy(chr, args.region);
        char *reg_chr = (char *) hts_parse_reg(chr, &start, &end);
        // start and end now zero-based end exclusive
        if (reg_chr) {
            *reg_chr = '\0';  // sets chr to be terminated at correct point
        } else {
            fprintf(stderr, "ERROR: Failed to parse region: '%s'.\n", args.region);
            exit(EXIT_FAILURE);
        }
        // simplify things for later (motif matching) on by fetching whole chr
        int len;
        char *ref = fai_fetch(fai, chr, &len);
        if (len < 0) {
            fprintf(stderr, "ERROR: Failed to fetch reference region: '%s'.\n", args.region);
            exit(EXIT_FAILURE);
        }
        if (!args.mask) {
            for (size_t i=0; i<len; ++i){ ref[i] = toupper(ref[i]); }
        }
        end = min(end, len);
        process_region(args, chr, start, end, ref, bed_files);

        free(chr);
        free(ref);
    }
    close_bed_files(bed_files);
    fai_destroy(fai);
    clock_t end = clock();
    fprintf(stderr, "Total time: %fs\n", (double)(end - begin) / CLOCKS_PER_SEC);
    exit(EXIT_SUCCESS);
}
