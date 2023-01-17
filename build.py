import itertools
import os

from cffi import FFI

dir_path = os.path.dirname(os.path.realpath(__file__))
src_dir='src'
libraries=['m', 'lzma', 'bz2', 'pthread', 'curl', 'crypto']
library_dirs=[]
print("WITHDEFLATE:", os.getenv('WITHDEFLATE'))
if os.getenv('WITHDEFLATE') == "1":
    print("Using deflate")
    libraries.append('deflate')
    library_dirs.append(os.path.join(dir_path, 'libdeflate'))

ffibuilder = FFI()
ffibuilder.set_source("libmodbampy",
    r"""
    #include "htslib/sam.h"
    #include "bamiter.h"
    #include "common.h"
    #include "counts.h"

    """,
    libraries=libraries,
    library_dirs=library_dirs,
    include_dirs=[src_dir, 'htslib'],
    extra_compile_args=['-std=c99', '-msse3', '-O3'],
    extra_objects=[
        'pymod.a',
        os.path.join('htslib', 'libhts.a')]
)

cdef = ["""
    // START: custom header

    // export free
    void free(void *ptr);

    typedef int64_t hts_pos_t;

    // basic bam opening/handling
    typedef struct bam1_core_t {
        hts_pos_t pos;
        int32_t tid;
        uint16_t bin; // NB: invalid on 64-bit pos
        uint8_t qual;
        uint8_t l_extranul;
        uint16_t flag;
        uint16_t l_qname;
        uint32_t n_cigar;
        int32_t l_qseq;
        int32_t mtid;
        hts_pos_t mpos;
        hts_pos_t isize;
    } bam1_core_t;


    typedef struct bam1_t {
        bam1_core_t core;
        uint64_t id;
        uint8_t *data;
        int l_data;
        uint32_t m_data;
        uint32_t mempolicy:2, :30 /* Reserved */;
    } bam1_t;

    bam1_t *bam_init1();
    void bam_destroy1(bam1_t *b);
    typedef struct mplp_data {...;} mplp_data;

    // opening bam with idx and hdr info
    typedef struct { ...; } bam_fset;
    bam_fset* create_bam_fset(char* fname);
    void destroy_bam_fset(bam_fset* fset);
    typedef struct set_fsets {
        bam_fset **fsets;
        size_t n;
    } set_fsets;
    set_fsets *create_filesets(const char **bams);
    void destroy_filesets(set_fsets *s);

    mplp_data *create_bam_iter_data(
        const bam_fset* fset, const char *chr, int start, int end,
        const char *read_group, const char tag_name[2], const int tag_value);
    void destroy_bam_iter_data(mplp_data *data);
    // iterate a file
    int read_bam(void *data, bam1_t *b);
    // cigar parsing
    int *qpos2rpos(bam1_t *b);

    // things from htslib
    hts_pos_t bam_endpos(const bam1_t *b);

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

    // from common.h needed in functions in counts.h
    //typedef struct mod_base {...;} mod_base;

    // END: custom header
"""]

# add in some things from headers, removing directives
for header in ('src/common.h', 'src/counts.h'):
    with open(header, 'r') as fh:
        cdef.append("// START: {}".format(header))
        cdef.append(
            ''.join(
                x for x in fh.readlines()
                if not (x.startswith('#') or x.startswith("static inline int"))))
        cdef.append("// END: {}".format(header))

ffibuilder.cdef('\n\n'.join(cdef))


if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
