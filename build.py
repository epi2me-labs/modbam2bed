import itertools
import os

from cffi import FFI

libraries=['m', 'z', 'lzma', 'bz2', 'pthread', 'curl', 'crypto']
library_dirs=[]
src_dir='src'

ffibuilder = FFI()
ffibuilder.set_source("libmodbampy",
    r"""
    #include "htslib/sam.h"
    #include "args.h"
    #include "bamiter.h"
    #include "common.h"
    #include "counts.h"

    """,
    libraries=libraries,
    library_dirs=library_dirs,
    include_dirs=[src_dir, 'htslib'],
    #sources=[
    #    os.path.join(src_dir, x)
    #    for x in ('common.c', 'bamiter.c', 'counts.c', 'args.c')],
    extra_compile_args=['-std=c99', '-msse3', '-O3'],
    extra_objects=[
        'pymod.a',
        os.path.join('htslib', 'libhts.a')]
)

cdef = ["""
    // START: custom header

    // export free
    void free(void *ptr);

    // basic bam opening/handling
    typedef struct bam1_core_t {uint16_t flag; uint32_t l_qseq; ...;} bam1_core_t;
    typedef struct bam1_t {bam1_core_t core; uint8_t *data; ...;} bam1_t;
    typedef struct mplp_data {...;} mplp_data;
    mplp_data *create_bam_iter_data(
        const char *bam_file, const char *chr, int start, int end,
        const char *read_group, const char tag_name[2], const int tag_value);
    void destroy_bam_iter_data(mplp_data *data);
    // iterate a file
    int read_bam(void *data, bam1_t *b);
    // cigar parsing
    int *qpos2rpos(bam1_t *b);

    // retrieving mod data
    typedef struct hts_base_mod_state hts_base_mod_state;
    hts_base_mod_state *hts_base_mod_state_alloc();
    void hts_base_mod_state_free(hts_base_mod_state *state);
    int bam_parse_basemod(const bam1_t *b, hts_base_mod_state *state);

    typedef struct hts_base_mod {
        int modified_base;
        int canonical_base;
        int strand;
        int qual;
    } hts_base_mod;
    int bam_next_basemod(
        const bam1_t *b, hts_base_mod_state *state,
        hts_base_mod *mods, int n_mods, int *pos);

    // END: custom header
"""]

# add in some things from headers, removing directives
for header in ('src/args.h', 'src/counts.h'):
    with open(header, 'r') as fh:
        cdef.append("// START: {}".format(header))
        cdef.append(
            ''.join(x for x in fh.readlines() if not x.startswith('#')))
        cdef.append("// END: {}".format(header))

ffibuilder.cdef('\n\n'.join(cdef))


if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
