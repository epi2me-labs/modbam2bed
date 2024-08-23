"""Functionality for interacting with modified base tags in BAM files."""

import argparse
import collections

import numpy as np

import libmodbampy

# remember to bump version in src/version.h too
__version__ = "0.10.0"
ffi = libmodbampy.ffi
libbam = libmodbampy.lib

MAX_MODS = 256  # from htslib

ModInfo = collections.namedtuple(
    'ModInfo', (
        'query_name', 'rpos', 'qpos', 'strand', 'mstrand',
        'cbase', 'mbase', 'qual'))


def _tidy_args(read_group, tag_name, tag_value):
    """Turn Python variables into CFFI ones."""
    if read_group is None:
        read_group = ffi.NULL
    else:
        read_group = ffi.new("char[]", read_group.encode())
    if tag_name is None:
        tag_name = ffi.new("char[2]", "".encode())
        tag_value = 0
    elif len(tag_name) != 2:
        raise ValueError("'tag_name' must be a length-2 string.")
    else:
        tag_name = ffi.new("char[2]", tag_name.encode())
    return read_group, tag_name, tag_value


class ModBase:
    """Helper to create a mod_base instance.

    :param code: modified base ChEBI code (e.g. "h" or 104)
    :param base: one of {A, C, G, T}
    :param name: long name of modified base (e.g. "5-methylcytosine")
    :param abbrev: short name of modified base (e.g. "5mC")

    Actually just a compatible list is created. Reuses the predefined
    instances from header where possible.
    """

    def __init__(self, code, base=None, name="unknown", abbrev="unknown"):
        """Initialise the instance."""
        self._name = ffi.new("char[]", name.encode())
        self._abbrev = ffi.new("char[]", abbrev.encode())
        self._base = base
        err = TypeError(
            "'base' should be a single character or None")
        if isinstance(self._base, str):
            if len(self._base) != 1:
                raise err
            self._base = base.encode()
            self._base_i = {"A": 1, "C": 2, "G": 4, "T": 8}[base]
        elif self._base is not None:
            raise err

        err = TypeError(
            "'code' should be a single character or an "
            "integer (ChEBI) code.")
        self._code = code
        if isinstance(self._code, str):
            # ffi won't coerce a char to int, so we need to do it
            if len(self._code) != 1:
                raise err
            self._code = ord(self._code)
        elif not isinstance(self._code, int):
            raise err

    @property
    def struct(self):
        """Return a list compatible with C structure."""
        for i in range(libbam.n_mod_bases):
            if libbam.mod_bases[i].code == self._code:
                return libbam.mod_bases[i]

        # make a new mod_base using a code and a canonical base
        if self._base is None:
            raise ValueError(
                f"Modified base type '{self._code}' unknown. Please provide "
                "a value for 'base' to describe the unmodified base.")
        mod_base_type = [
            self._name, self._abbrev,
            self._base, self._base_i, self._code]
        return mod_base_type


class ModBam:
    """A minimal class to iterate over a bam."""

    def __init__(self, bam):
        """Open a BAM file.

        :param bam: BAM file to open.
        """
        self.bam = bam
        self._bam_fset = ffi.gc(
            libbam.create_bam_fset(self.bam.encode()),
            libbam.destroy_bam_fset)

    def __enter__(self):
        """Open context."""
        return self

    def __exit__(self, type, value, traceback):
        """Exit context."""
        pass

    def reads(
            self, chrom, start, end,
            read_group=None, tag_name=None, tag_value=None, min_mapq=0):
        """Iterate over (filtered) alignments in file.

        :param chrom: reference sequence from BAM.
        :param start: reference start coordinate.
        :param end: reference end coordinate.
        :param read group: read group of read to return.
        :param tag_name: read tag to check during read filtering.
        :param tag_value: tag value for reads to keep.
        :param min_mapq: minimum read mapping quality.
        """
        read_group, tag_name, tag_value = _tidy_args(
            read_group, tag_name, tag_value)

        it = libbam.create_bam_iter_data(
                self._bam_fset, chrom.encode(), start, end,
                read_group, tag_name, tag_value, min_mapq)
        if it == ffi.NULL:
            return

        data = ffi.gc(it, libbam.destroy_bam_iter_data)
        mod_state = ffi.gc(
            libbam.hts_base_mod_state_alloc(),
            libbam.hts_base_mod_state_free)

        bam1_t = ffi.gc(libbam.bam_init1(), libbam.bam_destroy1)
        while libbam.read_bam(data, bam1_t) > 0:
            yield ModRead(bam1_t, mod_state)

    def pileup(
            self, chrom, start, end,
            read_group=None, tag_name=None, tag_value=None,
            low_threshold=0.33, high_threshold=0.66, threshold=0.66,
            mod_base="m", max_depth=None, canon_base=None, combine=False,
            min_mapq=0):
        """Create a base count matrix.

        :param chrom: reference sequence from BAM.
        :param start: reference start coordinate.
        :param end: reference end coordinate.
        :param read group: read group of read to return.
        :param tag_name: read tag to check during read filtering.
        :param tag_value: tag value for reads to keep.
        :param threshold: probability filter threshold for excluding
            calls from counts.
        :param mod_base: ChEBI code of modified base to examine.
        :param max_depth: maximum read depth to examine.
        :param canon_base: canonical base corresponding to `mod_base`.
            Required only if `mod_base` is not a modification known to
            the code.
        :param combine: combine (include) all alternative modifications
            with the same parent canonical base.
        :param min_mapq: minimum read mapping quality.
        """
        for thresh in (low_threshold, high_threshold):
            if thresh < 0.0 or thresh > 1.0:
                raise ValueError("Thresholds should be in (0,1).")
        threshold = int(threshold * 255.0)
        # C code currently uses high_threshold as the only threshold
        high_threshold = threshold
        read_group, tag_name, tag_value = _tidy_args(
            read_group, tag_name, tag_value)

        if max_depth is None:
            max_depth = libbam._INT_MAX

        _f = ffi.new("bam_fset *[]", [self._bam_fset])
        fsets = ffi.new("set_fsets *", {"fsets": _f, "n": 1})
        mod_base = ModBase(code=mod_base, base=canon_base)
        plp_data = libbam.calculate_pileup(
            fsets, chrom.encode(), start, end,
            read_group, tag_name, tag_value,
            threshold, mod_base.struct,
            combine, max_depth, min_mapq)
        # TODO: check for NULL

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


class ModRead:
    """Proxy for a bam alignment.

    The class is not intended to be instantiated by users.
    """

    def __init__(self, bam1_t, mod_state, header=None):
        """Create an interface to alignment.

        The input alignment is copied.
        """
        self._bam1_t = ffi.gc(libbam.bam_init1(), libbam.bam_destroy1)
        libbam.bam_copy1(self._bam1_t, bam1_t)
        self._mod_state = mod_state
        self._header = header

    @property
    def flags(self):
        """Return alignment flags."""
        return self._bam1_t.core.flag

    @property
    def is_unmapped(self):
        """Return if read is unmapped."""
        return self._bam1_t.core.flag & 4 > 0

    @property
    def is_reverse(self):
        """Return if alignment is to reverse strand."""
        return self._bam1_t.core.flag & 16 > 0

    @property
    def is_secondary(self):
        """Return if alignment is a secondary alignment."""
        return self._bam1_t.core.flag & 256 > 0

    @property
    def is_supplementary(self):
        """Return is alignment is a supplementary alignment."""
        return self._bam1_t.core.flag & 2048 > 0

    @property
    def mapping_quality(self):
        """Return mapping quality."""
        return self._bam1_t.core.qual

    @property
    def strand(self):
        """Return strand as '+' or '-'."""
        return "+-"[self.is_reverse]

    @property
    def query_name(self):
        """Return query name."""
        return ffi.string(
            (ffi.cast("char*", self._bam1_t.data))).decode()

    @property
    def query_length(self):
        """Return query length as record in BAM. See `query_sequence`."""
        return self._bam1_t.core.l_qseq

    @property
    def query_sequence(self):
        """Return the query sequence as recorded in the BAM.

        Includes soft-clipped bases, does not include hard-clipped bases, and
        may return an error when sequence is not recorded.
        """
        # bam1_seq() define
        # (b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname)
        raise NotImplementedError("query_sequence not implemented")

    @property
    def query_qualities(self):
        """Return the query quality array.

        Includes soft-clipped bases as for `query_sequence`.
        """
        # bam1_qual define
        # ((b)->data + ((b)->core.n_cigar<<2)
        #   + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1))
        raise NotImplementedError("query_qualities not implemented")

    @property
    def reference_name(self):
        """Return the reference name associated with the alignment."""
        if self._bam1_t.core.tid == -1:
            return None
        elif self.header is None:
            raise IndexError(
                "Require header information to retrieve reference_name")
        else:
            raise NotImplementedError(
                "Fetching reference_name not implemented")

    @property
    def reference_start(self):
        """Return the 0-based start position of the alignment."""
        return self._bam1_t.core.pos

    @property
    def reference_end(self):
        """Return the 0-based (exclusive) end position of the alignment."""
        return libbam.bam_endpos(self._bam1_t)

    @property
    def reference_length(self):
        """Return the length of the alignment on the reference."""
        return self.reference_end - self.reference_start

    @property
    def get_aligned_pairs(self):
        """Return aligned query and reference positions."""
        raise NotImplementedError("get_aligned_pairs not implemented")

    @property
    def alignment(self):
        """Create array representing alignment.

        The returned item is of length self.query_length
        """
        if not hasattr(self, "_alignment"):
            self._alignment = ffi.gc(
                libbam.qpos2rpos(self._bam1_t), libbam.free)
        return self._alignment

    @property
    def mod_sites(self):
        """Iterate over all modified bases in read.

        :yields: (read_id, ref. pos., query pos., ref. strand,
            mod. strand, canon. base, mod. base, mod. quality)

        The ref. strand is that recorded in the Mm tag from the bam.
        """
        mods = ffi.new("hts_base_mod[{}]".format(MAX_MODS))
        pos = ffi.new("int *")
        align = self.alignment
        libbam.bam_parse_basemod(self._bam1_t, self._mod_state)
        n = 1
        while n > 0:
            n = libbam.bam_next_basemod(
                self._bam1_t, self._mod_state, mods, MAX_MODS, pos)
            rpos = align[pos[0]]
            if n > 0:
                for i in range(n):
                    m = mods[i]
                    # note m.strand refers to strand recorded in the Mm tag.
                    modbase = m.modified_base
                    if modbase > 0:
                        modbase = chr(modbase)
                    else:
                        modbase = -modbase
                    yield ModInfo(
                        self.query_name, rpos, pos[0], self.strand, m.strand,
                        chr(m.canonical_base), modbase, m.qual)


def main():
    """Test entry point."""
    parser = argparse.ArgumentParser(
        description="Modified base demo program.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "bam", help="Indexed .bam file.")
    parser.add_argument(
        "chrom", help="Chromosome for which to fetch read")
    parser.add_argument(
        "start", type=int,
        help="Reference start coordinate.")
    parser.add_argument(
        "end", type=int,
        help="Reference end coordinate.")
    parser.add_argument(
        "--pileup", action="store_true",
        help="Create pileup counts rather than per-read modified base data")
    parser.add_argument(
        "--mod_base", default="m",
        help="Modified base to count during pileup.")
    parser.add_argument(
        "--low_threshold", type=float, default=0.33,
        help="Lower threshold for filtering.")
    parser.add_argument(
        "--high_threshold", type=float, default=0.66,
        help="High threshold for filtering.")
    args = parser.parse_args()

    with ModBam(args.bam) as bam:
        if args.pileup:
            codes = ffi.string(libbam.plp_bases).decode()
            print("pos\t", end="")
            print("\t".join(x for x in codes))
            positions, counts = bam.pileup(
                args.chrom, args.start, args.end, mod_base=args.mod_base,
                low_threshold=args.low_threshold,
                high_threshold=args.high_threshold)
            for p, row in zip(positions, counts):
                print(p, end='\t')
                print("\t".join(str(x) for x in row))
        else:
            counts = collections.Counter()
            for read in bam.reads(args.chrom, args.start, args.end):
                for pos_mod in read.mod_sites:
                    counts[pos_mod.qual] += 1
            for k in sorted(counts.keys()):
                print(k, counts[k])
