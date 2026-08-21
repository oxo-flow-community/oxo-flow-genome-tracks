#!/usr/bin/env python3
"""Deterministic single-cell fixtures for the oxo-flow-genome-tracks sc path.

Generates (stdlib only — no pysam/samtools needed):

  sc_bams/sc1.bam, sc_bams/sc2.bam   valid BGZF/BAM files, coordinate-sorted,
                                     reads mapped to chr1:1000-2000, each read
                                     carrying a CB cell-barcode tag
  sc_metadata/sc1.tsv, sc_metadata/sc2.tsv
                                     2-column barcode<TAB>group files (no
                                     header), the sinto filterbarcodes input

Fixture semantics (exercise the upstream touch-empty-bam fallback):
  * sc1 has cells in group g1 AND group g2  -> both split BAMs get real reads
  * sc2 has cells in group g1 ONLY          -> sc2/g2.bam exists only as the
    touched-empty file, so the downstream samtools merge still succeeds

Run from anywhere:  python3 test/fixtures/make_sc_fixtures.py
Re-running overwrites the outputs; the script self-verifies by reading the
BAMs back and asserting read counts, CB tags and sort order.
"""

import binascii
import os
import struct
import sys
import zlib

# --- constants (keep in sync with main.oxoflow [config] + [[values]]) --------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SC_BAM_DIR = os.path.join(SCRIPT_DIR, "sc_bams")
SC_META_DIR = os.path.join(SCRIPT_DIR, "sc_metadata")

# canonical BGZF EOF marker (what bgzip writes; 28 bytes)
BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")

# chr1 length must cover the reads (all in the 1000-2000 region)
REF_NAME = b"chr1"
REF_LEN = 200000
SAM_HEADER = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:200000\n"

# per-sample cell barcodes: (barcode, group)
SC_CELLS = {
    "sc1": [
        ("AAAAAAAAAAAAAA01-1", "g1"),
        ("AAAAAAAAAAAAAA02-1", "g1"),
        ("AAAAAAAAAAAAAA03-1", "g2"),
        ("AAAAAAAAAAAAAA04-1", "g2"),
    ],
    "sc2": [
        ("TTTTTTTTTTTTTT01-1", "g1"),
        ("TTTTTTTTTTTTTT02-1", "g1"),
    ],
}

READS_PER_CELL = 2
READ_LEN = 50  # even, so the packed sequence is byte-aligned
READ_START = 1000  # 0-based; reads span chr1:1001-2050 in 1-based coords


def bgzf_block(data: bytes) -> bytes:
    """Wrap a data chunk in a BGZF block (gzip header + BC/BS extra + raw deflate).

    Extra field layout (htslib): subfield tag "BC" (2 bytes) + subfield size (2
    bytes, = 2) + BS (2 bytes). NOTE: htslib stores and reads BS as
    (actual block size - 1) — bgzf_read_block does unpackInt16(&h[16]) + 1
    ("+1 because when writing this number, we used -1") — so a canonical
    28-byte EOF block stores 27. Store size-1 here or htslib shifts the
    CRC/ISIZE trailer by one byte and rejects the block.
    """
    comp = zlib.compressobj(6, zlib.DEFLATED, -15)
    payload = comp.compress(data) + comp.flush()
    block_size = 26 + len(payload)  # 18 header+extra + payload + 8 trailer
    header = b"\x1f\x8b\x08\x04" + b"\x00\x00\x00\x00" + b"\x00\xff"
    header += struct.pack("<H", 6)  # XLEN
    header += b"BC" + struct.pack("<H", 2) + struct.pack("<H", block_size - 1)
    return header + payload + struct.pack(
        "<II", binascii.crc32(data) & 0xFFFFFFFF, len(data)
    )


def bgzf_encode(raw: bytes) -> bytes:
    """Encode raw BAM bytes as BGZF (single block for tiny fixtures) + EOF marker."""
    return bgzf_block(raw) + BGZF_EOF


def reg2bin(beg: int, end: int) -> int:
    """htslib reg2bin — bin of a read spanning [beg, end) (0-based)."""
    end -= 1
    if beg >> 14 == end >> 14:
        return ((1 << 15) - 1) // 7 + (beg >> 14)
    if beg >> 17 == end >> 17:
        return ((1 << 12) - 1) // 7 + (beg >> 17)
    if beg >> 20 == end >> 20:
        return ((1 << 9) - 1) // 7 + (beg >> 20)
    if beg >> 23 == end >> 23:
        return ((1 << 6) - 1) // 7 + (beg >> 23)
    if beg >> 26 == end >> 26:
        return ((1 << 3) - 1) // 7 + (beg >> 26)
    return 0


BASE_CODES = {"A": 1, "C": 2, "G": 4, "T": 8, "N": 15}


def pack_seq(seq: str) -> bytes:
    """Pack a 4-bit-per-base sequence (A=1, C=2, G=4, T=8)."""
    if len(seq) % 2:
        raise ValueError("sequence length must be even")
    out = bytearray()
    for i in range(0, len(seq), 2):
        out.append((BASE_CODES[seq[i]] << 4) | BASE_CODES[seq[i + 1]])
    return bytes(out)


def make_alignment(
    read_id: int, pos: int, qname: bytes, flag: int, seq: str, barcode: str
) -> bytes:
    """Encode one BAM alignment record (1 M cigar, CB barcode tag).

    Layout follows htslib's bam_write1 exactly: a 4-byte l_data prefix
    (= 32 + qname + cigar + seq + qual + aux, padding NOT written to disk —
    htslib synthesizes the qname padding in memory on read), then the 32-byte
    core with the UNPADDED l_read_name, then qname+NUL, cigar, packed seq,
    quals, aux tags. bam_read1 rejects records whose prefix is < 32.
    """
    l_read_name = len(qname) + 1
    n_cigar = 1
    l_seq = len(seq)
    core = struct.pack(
        "<iiBBHHHiiii",
        0,  # refID (chr1)
        pos,
        l_read_name,
        60,  # MAPQ
        reg2bin(pos, pos + l_seq),
        n_cigar,
        flag,
        l_seq,
        -1,  # next_refID
        -1,  # next_pos
        0,  # tlen
    )
    cigar = struct.pack("<i", (l_seq << 4) | 0)  # l_seq M
    aux = b"CB" + b"Z" + barcode.encode("ascii") + b"\x00"
    record = core + qname + b"\x00" + cigar + pack_seq(seq) + (b"\xff" * l_seq) + aux
    return struct.pack("<i", len(record)) + record


def build_bam(reads: list[bytes]) -> bytes:
    """Assemble BAM bytes: magic + header + refs + alignments."""
    out = b"BAM\x01"
    out += struct.pack("<i", len(SAM_HEADER)) + SAM_HEADER
    out += struct.pack("<i", 1)  # n_ref
    out += struct.pack("<i", len(REF_NAME) + 1) + REF_NAME + b"\x00"
    out += struct.pack("<i", REF_LEN)
    for rec in reads:
        out += rec
    return out


def generate_sample(sample: str, seq: str) -> tuple[bytes, bytes]:
    """BAM bytes + metadata TSV bytes for one sc sample (deterministic)."""
    reads = []
    lines = []
    for i, (barcode, group) in enumerate(SC_CELLS[sample]):
        lines.append(f"{barcode}\t{group}")
        for r in range(READS_PER_CELL):
            pos = READ_START + i * 250 + r * 60  # spread reads across the region
            qname = f"R{i}_{r}_{sample}".encode("ascii")
            reads.append(make_alignment(len(reads), pos, qname, 0, seq, barcode))
    # coordinate-sorted (already in order); header declares SO:coordinate
    return build_bam(reads), "\n".join(lines).encode("ascii") + b"\n"


def write_bytes(path: str, data: bytes) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        f.write(data)


# ---------------------------------------------------------------------------
# self-verification: read the BAM back and parse it (BGZF + BAM binary)
# ---------------------------------------------------------------------------

def read_bgzf(path: str) -> bytes:
    raw = open(path, "rb").read()
    assert raw.endswith(BGZF_EOF), "missing BGZF EOF marker"
    body = raw[: -len(BGZF_EOF)]
    data = bytearray()
    view = memoryview(body)
    while view:
        assert view[:3] == b"\x1f\x8b\x08", "bad BGZF magic"
        flags = view[3]
        assert flags & 0x04, "BGZF block without FEXTRA"
        xlen = struct.unpack_from("<H", view, 10)[0]
        total = 12 + xlen
        assert view[12:14] == b"BC", "BGZF block without BC subfield"
        # BS is stored as actual block size - 1 (htslib adds 1 on read)
        stored_bs = struct.unpack_from("<H", view, 16)[0]
        block_size = stored_bs + 1
        assert total + 8 <= block_size <= len(view), "truncated BGZF block"
        crc = struct.unpack_from("<I", view, block_size - 8)[0]
        isize = struct.unpack_from("<I", view, block_size - 4)[0]
        payload = bytes(view[total : block_size - 8])
        chunk = zlib.decompress(payload, -15)
        # BGZF CRC32 is over the uncompressed data
        assert binascii.crc32(chunk) & 0xFFFFFFFF == crc, "CRC32 mismatch"
        assert len(chunk) == isize, "ISIZE mismatch"
        data += chunk
        view = view[block_size:]
    return bytes(data)


def parse_bam(data: bytes) -> tuple[int, int, dict]:
    """Return (n_reads, ref_len, {read_name: [cb tags]}) — asserts as it walks."""
    assert data[:4] == b"BAM\x01", "bad BAM magic"
    off = 4
    (l_text,) = struct.unpack_from("<i", data, off)
    off += 4 + l_text
    (n_ref,) = struct.unpack_from("<i", data, off)
    off += 4
    ref_len = None
    for _ in range(n_ref):
        (l_name,) = struct.unpack_from("<i", data, off)
        off += 4 + l_name
        (l_ref,) = struct.unpack_from("<i", data, off)
        ref_len = l_ref
        off += 4
    n_reads = 0
    tags = {}
    prev_pos = -1
    while off < len(data):
        # each record is prefixed by its l_data (32 + payload); walk by prefix
        # like bam_read1 does
        (block_len,) = struct.unpack_from("<i", data, off)
        assert 32 <= block_len <= len(data) - off - 4, "bad record l_data"
        rec = data[off + 4 : off + 4 + block_len]
        off += 4 + block_len
        (ref_id, pos, l_read_name, _mapq, _bin, n_cigar, _flag, l_seq,
         _nref, _npos, _tlen) = struct.unpack_from("<iiBBHHHiiii", rec, 0)
        r = 32
        qname = bytes(rec[r : r + l_read_name - 1]).decode("ascii")
        r += l_read_name
        r += n_cigar * 4
        r += (l_seq + 1) // 2
        r += l_seq
        assert ref_id == 0 and pos >= prev_pos, "reads not coordinate-sorted"
        prev_pos = pos
        assert rec[r : r + 3] == b"CBZ", "expected CB:Z tag"
        end = rec.index(b"\x00", r + 3)
        cb = bytes(rec[r + 3 : end]).decode("ascii")
        tags.setdefault(qname, []).append(cb)
        n_reads += 1
    return n_reads, ref_len, tags


def main() -> int:
    seq = "A" * READ_LEN
    expected_reads = sum(READS_PER_CELL * len(cells) for cells in SC_CELLS.values())
    total = 0
    for sample, cells in SC_CELLS.items():
        bam, tsv = generate_sample(sample, seq)
        bam_path = os.path.join(SC_BAM_DIR, f"{sample}.bam")
        tsv_path = os.path.join(SC_META_DIR, f"{sample}.tsv")
        write_bytes(bam_path, bgzf_encode(bam))
        write_bytes(tsv_path, tsv)
        n_reads, ref_len, tags = parse_bam(read_bgzf(bam_path), )
        expected = READS_PER_CELL * len(cells)
        assert n_reads == expected, f"{sample}: {n_reads} reads != {expected}"
        assert ref_len == REF_LEN, f"{sample}: ref len {ref_len}"
        assert len(tags) == expected, f"{sample}: {len(tags)} named reads"
        with open(tsv_path, "rb") as f:
            tsv_text = f.read().decode("ascii")
        assert len(tsv_text.splitlines()) == len(cells), f"{sample}: TSV rows"
        assert all(t.split("\t")[1] in {"g1", "g2"} for t in tsv_text.splitlines())
        assert set(c[0] for c in cells) == set(
            tag for tags_ in tags.values() for tag in tags_
        ), f"{sample}: CB tags vs TSV barcodes"
        total += n_reads
        print(f"  {sample}: {n_reads} reads, {len(cells)} cells -> {bam_path}")
    print(f"OK — {total} reads total; fixture files written to {SC_BAM_DIR} and {SC_META_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
