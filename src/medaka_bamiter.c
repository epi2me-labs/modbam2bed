#include <errno.h>
#include <string.h>

#include "medaka_bamiter.h"

// iterator for reading bam
int read_bam(void *data, bam1_t *b) {
    mplp_data *aux = (mplp_data*) data;
    bool have_rg = (aux->read_group != NULL);
    uint8_t *rg;
    char *rg_val;
    int ret;
    while (1) {
        ret = aux->iter ? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
        if (ret<0) break;
        // only take primary alignments
        if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FDUP)) continue;
        // filter by mapping quality
        if ((int)b->core.qual < aux->min_mapQ) continue;
        // filter by RG (read group):
        if (have_rg) {
            rg = bam_aux_get((const bam1_t*) b, "RG");
            if (rg == NULL) continue;  // missing
            rg_val = bam_aux2Z(rg);
            if (errno == EINVAL) continue;  // bad parse
            if (strcmp(aux->read_group, rg_val) != 0) continue;  // not wanted
        }
        break;
    }
    return ret;
}
