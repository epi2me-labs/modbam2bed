import argparse

import numpy as np
import libmodbampy

__version__ = "0.0.1"
ffi = libmodbampy.ffi
libbam = libmodbampy.lib


def _tidy_args(read_group, tag_name, tag_value):
    if read_group is None:
        read_group = ffi.NULL
    else:
        read_group = ffi.new("char[]", read_group.encode())
    if tag_name is None:
        tag_name = ffi.new("char[2]", "".encode())
        tag_value = 0
        keep_missing = False
    elif len(tag_name) != 2:
        raise ValueError("'tag_name' must be a length-2 string.")
    else:
        tag_name = ffi.new("char[2]", tag_name.encode())
    return read_group, tag_name, tag_value



class ModBam:
    """A minimal class to iterate over a bam."""

    def __init__(self, bam_file, chrom, start, end,
            read_group=None, tag_name=None, tag_value=None):

        read_group, tag_name, tag_value = _tidy_args(
            read_group, tag_name, tag_value)

        self._bam1_t = ffi.new("bam1_t *")
        self._data = ffi.gc(
            libbam.create_bam_iter_data(
                bam_file.encode(), chrom.encode(), start, end,
                read_group, tag_name, tag_value),
            libbam.destroy_bam_iter_data)
        self._mod_state = ffi.gc(
            libbam.hts_base_mod_state_alloc(),
            libbam.hts_base_mod_state_free)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass

    def reads(self):
        rtn = 1
        while rtn > 0:
            rtn = libbam.read_bam(self._data, self._bam1_t)
            yield ModRead(self._bam1_t, self._mod_state)


class ModRead:
    """Proxy for a bam alignment.

    The class is not intended to be instantiated by users.
    """

    def __init__(self, bam1_t, mod_state):
        self._bam1_t = bam1_t
        self._mod_state = mod_state
        self.read_id = ffi.string((ffi.cast("char*", self._bam1_t.data))).decode()
        self.flags = self._bam1_t.core.flag
        self.is_reverse = self._bam1_t.core.flag & 16 > 0
        self.strand = '+-'[self.is_reverse]

    def alignment(self):
        return ffi.gc(libbam.qpos2rpos(self._bam1_t), libbam.free)

    def mod_sites(self):
        """Iterator over all modified bases in read."""
        mods = ffi.new("hts_base_mod[256]")
        pos = ffi.new("int *")
        align = self.alignment()
        libbam.bam_parse_basemod(self._bam1_t, self._mod_state);
        n = 1
        while n > 0:
            n = libbam.bam_next_basemod(self._bam1_t, self._mod_state, mods, 256, pos)
            rpos = align[pos[0]]
            if n > 0:
                for i in range(n):
                    m = mods[i]
                    # note m.strand refers to the strand recorded in the Mm tag.
                    yield self.read_id, rpos, pos[0], self.strand, m.strand, chr(m.canonical_base), chr(m.modified_base), m.qual


def pileup(
        bam, chrom, start, end,
        read_group=None, tag_name=None, tag_value=None,
        low_threshold=0.33, high_threshold=0.66, mod_base="m"):

    for thresh in (low_threshold, high_threshold):
        if thresh < 0.0 or thresh > 1.0:
            raise ValueError("Thresholds should be in (0,1).")
    low_threshold, high_threshold = (
        int(x * 255.0) for x in (low_threshold, high_threshold))
    read_group, tag_name, tag_value = _tidy_args(
        read_group, tag_name, tag_value)
    plp_data = libbam.calculate_pileup(
        bam.encode(), chrom.encode(), start, end, read_group, tag_name, tag_value,
        low_threshold, high_threshold, mod_base.encode())

    # copy data to numpy, we could be more clever here an wrap
    #   the pointer in a subclass of ndarray to track its lifetime
    #   and avoid the explicit copy
    n_rows = libbam.featlen
    size_sizet = np.dtype(np.uintp).itemsize
    np_counts = np.frombuffer(ffi.buffer(
        plp_data.matrix, size_sizet * plp_data.n_cols * n_rows),
        dtype=np.uintp
    ).reshape(plp_data.n_cols, n_rows).copy()
    np_positions = np.frombuffer(
        ffi.buffer(plp_data.major, size_sizet * plp_data.n_cols),
        dtype=np.uintp).copy()
    libbam.destroy_plp_data(plp_data)
    return np_positions, np_counts


def main():
    parser = argparse.ArgumentParser(description="Modified base demo program.")
    parser.add_argument("bam", help="Indexed .bam file.")
    parser.add_argument("chrom", help="Chromosome for which to fetch read")
    parser.add_argument("start", type=int, help="Reference start coordinate.")
    parser.add_argument("end", type=int, help="Reference end coordinate.")
    parser.add_argument("--pileup", action="store_true",
        help="Create pileup counts rather than per-read modified base data")
    parser.add_argument("--mod_base", default="m", help="Modified base to count during pileup.")
    args = parser.parse_args()

    if args.pileup:
        codes = ffi.string(libbam.plp_bases).decode()
        print("pos\t", end="")
        print("\t".join(x for x in codes))
        positions, counts = pileup(args.bam, args.chrom, args.start, args.end, mod_base=args.mod_base)
        for p, row in zip(positions, counts):
            print(p, end='\t')
            print("\t".join(str(x) for x in row))

    else:
        with ModBam(args.bam, args.chrom, args.start, args.end) as bam:
            for read in bam.reads():
                for pos_mod in read.mod_sites():
                    print(*pos_mod)
