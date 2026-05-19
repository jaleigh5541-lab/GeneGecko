#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 23 10:19:36 2025

@author: JessicaGeneWeaver

Merge DNA sequences using longest common substring (LCS) strategy:
- Forward = direct DNA sequence
- Reverse = reverse complemented DNA sequence
- Merge at longest common substring (LCS) ≥ min_lcs (default 12 bp)
- Translate merged DNA in six frames
- Extract longest ORF (M → * or uninterrupted from M if no stop codon)
- Save both merged DNA and longest AA ORF
- Compare longest ORF with reference amino acid sequences from refseq.xlsx
- Output summary with percent identity
- Save formatted text alignments

Reference Excel file requirements:
    refseq.xlsx with columns:
        "ID"  → sample prefix (before underscore, case-insensitive)
        "Sequence"       → reference amino acid sequence
"""

import os
import re
from typing import Union, Tuple, Dict, List, Any, Optional
from io import StringIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
import pandas as pd

# Feature detection imports
from feature_detection import (
    detect_features,
    map_features_to_gapped_alignment,
    build_feature_position_map,
    DetectedFeature,
    BUILTIN_FEATURES
)


def _read_ref_table(filepath: Union[str, Any]) -> "pd.DataFrame":
    """Read a reference file as a DataFrame, supporting .xlsx and .csv."""
    name = getattr(filepath, "name", str(filepath)).lower()
    try:
        if name.endswith(".csv"):
            return pd.read_csv(filepath)
        return pd.read_excel(filepath)
    except Exception as e:
        raise ValueError(f"Failed to read reference file: {e}")


# ───────────────────────────────────────────────────────────────
#  Sequence utilities
# ───────────────────────────────────────────────────────────────

def read_seq_file(file_or_path: Union[str, Any]) -> str:
    """
    Read DNA sequence from a file path or file-like object.

    Handles both FASTA format files and raw sequence files. Non-nucleotide
    characters are converted to 'N'.

    Args:
        file_or_path: File path (str) or file-like object with .read() method

    Returns:
        DNA sequence as uppercase string

    Raises:
        ValueError: If no valid sequence is found in the file
    """
    if hasattr(file_or_path, "read"):
        # file-like object - read content once
        content = file_or_path.read()
        if isinstance(content, bytes):
            content = content.decode('utf-8-sig')
        raw_lines = [line.strip() for line in content.splitlines() if line.strip()]
        content_str = content
    else:
        with open(file_or_path, "r", encoding="utf-8-sig") as f:
            content_str = f.read()
        raw_lines = [line.strip() for line in content_str.splitlines() if line.strip()]

    # -------------------------
    # CASE 1: FASTA header present
    # -------------------------
    if any(line.startswith(">") for line in raw_lines):
        try:
            records = list(SeqIO.parse(StringIO(content_str), "fasta"))
            if records:
                seq = str(records[0].seq).upper()
                if seq:
                    return seq
        except Exception:
            pass  # fallback to raw parsing

    # -------------------------
    # CASE 2: RAW SEQUENCE FILE
    # -------------------------
    seq_lines = []
    for line in raw_lines:
        # only keep characters that could be nucleotides
        clean = re.sub(r"[^ACGTNacgtn]", "N", line)
        if any(c in "ACGTNacgtn" for c in clean):
            seq_lines.append(clean.upper())

    seq = "".join(seq_lines)

    if not seq:
        file_repr = getattr(file_or_path, 'name', str(file_or_path))
        raise ValueError(f"No usable sequence found in {file_repr}")

    return seq


def longest_common_substring(s1: str, s2: str) -> Tuple[str, int, int]:
    """
    Find the longest common substring between two strings using dynamic programming.

    Args:
        s1: First string
        s2: Second string

    Returns:
        Tuple of (longest_common_substring, start_position_in_s1, end_position_in_s1)
    """
    m, n = len(s1), len(s2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    max_len = 0
    end_pos_s1 = 0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s1[i - 1] == s2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
                if dp[i][j] > max_len:
                    max_len = dp[i][j]
                    end_pos_s1 = i
    lcs = s1[end_pos_s1 - max_len:end_pos_s1] if max_len > 0 else ""
    return lcs, end_pos_s1 - max_len, end_pos_s1


def merge_by_lcs_dna(fwd_dna: str, rev_dna: str, min_lcs: int = 12) -> Dict[str, Any]:
    """
    Merge forward and reverse DNA sequences at the longest common substring (LCS).

    Finds the LCS between forward and reverse complement sequences and merges
    them by trimming the forward sequence at the LCS end and the reverse
    sequence at the LCS start.

    Args:
        fwd_dna: Forward DNA sequence
        rev_dna: Reverse complement DNA sequence
        min_lcs: Minimum LCS length required for merging (default: 12 bp)

    Returns:
        Dictionary with keys:
            - merged: Merged DNA sequence or None if LCS too short
            - lcs: The longest common substring found
            - fwd_trimmed_len: Length of forward sequence used
            - rev_trimmed_len: Length of reverse sequence used

    Raises:
        ValueError: If input sequences are empty
    """
    if not fwd_dna or not rev_dna:
        raise ValueError("Forward and reverse sequences cannot be empty")

    lcs, fwd_start, fwd_end = longest_common_substring(fwd_dna, rev_dna)
    if not lcs or len(lcs) < min_lcs:
        return {"merged": None, "lcs": "", "fwd_trimmed_len": 0, "rev_trimmed_len": 0}
    
    rev_pos = rev_dna.find(lcs)
    trimmed_fwd = fwd_dna[:fwd_end]
    trimmed_rev = rev_dna[rev_pos + len(lcs):]
    merged = trimmed_fwd + trimmed_rev
    
    return {
        "merged": merged,
        "lcs": lcs,
        "fwd_trimmed_len": len(trimmed_fwd),
        "rev_trimmed_len": len(trimmed_rev)
    }


def translate_six_frames(seq: str) -> Dict[str, str]:
    """
    Translate DNA sequence in all six reading frames.

    Args:
        seq: DNA sequence string

    Returns:
        Dictionary mapping frame names (+1, +2, +3, -1, -2, -3) to
        translated amino acid sequences
    """
    frames = {}
    dna_seq = Seq(seq)
    rev_seq = dna_seq.reverse_complement()
    for i in range(3):
        fwd_subseq = dna_seq[i:]
        trim_len = len(fwd_subseq) - len(fwd_subseq) % 3
        frames[f"+{i+1}"] = str(fwd_subseq[:trim_len].translate())
        rev_subseq = rev_seq[i:]
        trim_len = len(rev_subseq) - len(rev_subseq) % 3
        frames[f"-{i+1}"] = str(rev_subseq[:trim_len].translate())
    return frames


def longest_ORF(aa_seq: str) -> str:
    """
    Extract the longest open reading frame (ORF) from an amino acid sequence.

    An ORF is defined as starting with M (methionine) and ending at * (stop codon)
    or continuing to the end of the sequence if no stop codon is found.

    Args:
        aa_seq: Amino acid sequence string

    Returns:
        Longest ORF sequence as string (empty if no M found)
    """
    orfs = re.findall(r'M[^*]*', aa_seq)
    if orfs:
        return max(orfs, key=len)
    m_index = aa_seq.find("M")
    if m_index != -1:
        return aa_seq[m_index:]
    return ""


def find_all_orfs(aa_seq: str, min_length: int = 71) -> List[Tuple[str, int]]:
    """Non-overlapping ORFs (M→* or M→end) ≥ min_length in a translated sequence.
    Returns list of (orf_aa_string, start_index_in_aa_seq)."""
    orfs, pos = [], 0
    while pos < len(aa_seq):
        m = aa_seq.find('M', pos)
        if m == -1:
            break
        stop = aa_seq.find('*', m)
        orf = aa_seq[m:stop] if stop != -1 else aa_seq[m:]
        if len(orf) >= min_length:
            orfs.append((orf, m))
        pos = stop + 1 if stop != -1 else len(aa_seq)
    return orfs


def build_ungapped_to_gapped_map(gapped_seq: str) -> List[int]:
    """Returns list where result[i] = gapped position of the i-th non-gap character."""
    return [gi for gi, c in enumerate(gapped_seq) if c != '-']


def _find_orf_dna_span(dna: str, orf_aa: str) -> Tuple[Optional[int], Optional[int]]:
    """Return (dna_start, dna_end) in dna that encodes orf_aa (forward frames only)."""
    dna_upper = dna.upper()
    for offset in range(3):
        fragment = dna_upper[offset:]
        n = (len(fragment) // 3) * 3
        if n == 0:
            continue
        aa = str(Seq(fragment[:n]).translate())
        pos = aa.find(orf_aa)
        if pos >= 0:
            return offset + pos * 3, offset + (pos + len(orf_aa)) * 3
    return None, None


def inspect_single_read(
    raw_dna: str,
    label: str,
    ref_dna: Optional[str] = None,
    ref_aa: Optional[str] = None,
    use_dna_alignment: bool = True,
    query_circular: bool = False,
) -> Dict[str, Any]:
    """Align a single raw read to the reference and enumerate its ORFs.

    Returns a dict suitable for inclusion in result["raw_reads"].
    """
    alignment_html = None
    alignment_html_with_orfs = None
    identity = 0.0
    orf_count = 0
    best_orf_length = 0
    all_orfs_list = []
    aa_seq = ""

    try:
        if use_dna_alignment and ref_dna:
            alignment_html, identity, aa_seq, _, orf_data = align_dna_then_translate(
                raw_dna, ref_dna, label, query_circular=query_circular
            )
            orf_count = orf_data.get("orf_count", 0)
            alignment_html_with_orfs = orf_data.get("alignment_with_orfs")
        elif ref_aa:
            aa_seq, _ = get_longest_orf_six_frames(raw_dna)
            if aa_seq:
                alignment_html, identity = visualize_orf_vs_reference(
                    aa_seq, ref_aa, label, query_circular=query_circular
                )
        else:
            aa_seq, _ = get_longest_orf_six_frames(raw_dna)

        if aa_seq:
            orfs = find_all_orfs(aa_seq)
            all_orfs_list = [
                {"sequence": orf, "length": len(orf), "start": start}
                for orf, start in orfs
            ]
            best_orf_length = max((o["length"] for o in all_orfs_list), default=0)
            if not use_dna_alignment:
                orf_count = len(all_orfs_list)
    except Exception:
        pass

    return {
        "label": label,
        "dna": raw_dna,
        "length": len(raw_dna),
        "alignment_html": alignment_html,
        "alignment_html_with_orfs": alignment_html_with_orfs,
        "identity": identity,
        "orf_count": orf_count,
        "best_orf_length": best_orf_length,
        "all_orfs": all_orfs_list,
    }


def build_investigate_html(
    raw_reads: List[Dict[str, Any]],
    ref_dna: str,
    lcs: Optional[str] = None,
    block_size: int = 60,
    sample_name: str = "",
    include_orfs: bool = False,
) -> str:
    """Pileup-style HTML: one reference row, all reads aligned below it.

    Each read is projected into reference coordinates so the reference is
    displayed only once. LCS positions (if provided) are outlined in amber
    across all read rows to show the overlap region.
    """
    if not ref_dna or not raw_reads:
        return "<p>No alignment data available.</p>"

    ref_len = len(ref_dna)
    read_data = []
    fwd_orf_aa_row: Dict[int, str] = {}  # populated from forward read, reused for reverse

    for rd_idx, rd in enumerate(raw_reads):
        dna = rd.get("dna", "")
        if not dna:
            continue

        alns = pairwise2.align.localms(dna, ref_dna, 2, -1, -2, -0.5)
        if not alns:
            continue

        aln = alns[0]
        seqA_full = aln[0]
        seqB_full = aln[1]
        begin = int(aln[3])
        end   = int(aln[4])

        seqA_reg = seqA_full[begin:end]
        seqB_reg = seqB_full[begin:end]

        # Trim to region where reference has sequence (no leading/trailing ref gaps)
        r_start = next((i for i, c in enumerate(seqB_reg) if c != "-"), 0)
        r_end   = len(seqB_reg) - next((i for i, c in enumerate(reversed(seqB_reg)) if c != "-"), 0)
        seqA_trim = seqA_reg[r_start:r_end]
        seqB_trim = seqB_reg[r_start:r_end]

        # 0-based position in ref_dna where the trimmed alignment begins
        ref_offset = (len(seqB_full[:begin].replace("-", "")) +
                      len(seqB_reg[:r_start].replace("-", "")))

        # Build ref-coordinate array: read_arr[ref_pos] = read char or None
        read_arr = [None] * ref_len
        ref_p = ref_offset
        for qc, rc in zip(seqA_trim, seqB_trim):
            if ref_p >= ref_len:
                break
            if rc != "-":
                read_arr[ref_p] = qc  # "-" means deletion in read
                ref_p += 1
            # rc == "-" → insertion in read; no ref position consumed

        ref_end_pos = ref_p

        # Map LCS positions (in read-DNA coordinates) to ref coordinates
        lcs_ref_positions: set = set()
        if lcs and dna:
            lcs_hit = dna.upper().find(lcs.upper())
            if lcs_hit >= 0:
                lcs_span = (lcs_hit, lcs_hit + len(lcs))
                dna_idx = (len(seqA_full[:begin].replace("-", "")) +
                           len(seqA_reg[:r_start].replace("-", "")))
                ref_p2 = ref_offset
                for qc, rc in zip(seqA_trim, seqB_trim):
                    if rc != "-":
                        if qc != "-":
                            if lcs_span[0] <= dna_idx < lcs_span[1]:
                                lcs_ref_positions.add(ref_p2)
                            dna_idx += 1
                        ref_p2 += 1
                    else:
                        if qc != "-":
                            dna_idx += 1

        # Build dna_pos → ref_pos map to locate ORFs in reference coordinates
        dna_base = (len(seqA_full[:begin].replace("-", "")) +
                    len(seqA_reg[:r_start].replace("-", "")))
        dna_to_ref: Dict[int, int] = {}
        d_idx = dna_base
        r_p3 = ref_offset
        for qc, rc in zip(seqA_trim, seqB_trim):
            if rc != "-":
                if qc != "-":
                    dna_to_ref[d_idx] = r_p3
                    d_idx += 1
                r_p3 += 1
            else:
                if qc != "-":
                    d_idx += 1  # insertion in read, no ref position

        # orf_aa_row: ref_pos → display char for ORF annotation
        # Middle codon nt gets the AA letter; flanking nts get '▓'.
        orf_aa_row: Dict[int, str] = {}

        if rd_idx == 0:
            # Forward read: search forward frames via stored all_orfs
            for orf_info in rd.get("all_orfs", []):
                orf_aa = orf_info.get("sequence", "")
                if not orf_aa:
                    continue
                ds, de = _find_orf_dna_span(dna, orf_aa)
                if ds is None:
                    continue
                # Include stop codon (*) if the next codon is a stop
                orf_display = orf_aa
                if de + 3 <= len(dna):
                    next_aa = str(Seq(dna[de:de + 3].upper()).translate())
                    if next_aa == '*':
                        orf_display = orf_aa + '*'
                for aa_idx, aa_char in enumerate(orf_display):
                    for nt_off, marker in [(0, '▓'), (1, aa_char), (2, '▓')]:
                        dp = ds + aa_idx * 3 + nt_off
                        rp = dna_to_ref.get(dp)
                        if rp is not None and rp not in orf_aa_row:
                            orf_aa_row[rp] = marker
            fwd_orf_aa_row = orf_aa_row
        else:
            # Reverse read: dna is already the reverse-complemented T7-Term read
            # (RC was applied in process_pair before storing in raw_reads).
            # Search ALL 3 frames of the full dna — the ORF start (ATG) may lie
            # before the aligned region, so searching only aln_lo:aln_hi misses it.
            best_orf_aa: Optional[str] = None
            best_rc_start: Optional[int] = None
            for offset in range(3):
                fragment = dna[offset:]
                n = (len(fragment) // 3) * 3
                if n == 0:
                    continue
                aa = str(Seq(fragment[:n]).translate())
                for orf_seq, start_aa in find_all_orfs(aa, min_length=1):
                    if 'X' not in orf_seq:
                        # Include stop codon (*) if present immediately after ORF
                        stop_pos = start_aa + len(orf_seq)
                        if stop_pos < len(aa) and aa[stop_pos] == '*':
                            orf_with_stop = orf_seq + '*'
                        else:
                            orf_with_stop = orf_seq
                        if best_orf_aa is None or len(orf_with_stop) > len(best_orf_aa):
                            best_orf_aa = orf_with_stop
                            best_rc_start = offset + start_aa * 3

            if best_orf_aa and best_rc_start is not None:
                for aa_idx, aa_char in enumerate(best_orf_aa):
                    codon_start = best_rc_start + aa_idx * 3
                    for nt_off, marker in [(0, '▓'), (1, aa_char), (2, '▓')]:
                        dp = codon_start + nt_off
                        rp = dna_to_ref.get(dp)
                        if rp is not None and rp not in orf_aa_row:
                            orf_aa_row[rp] = marker

        read_data.append({
            "label":             rd["label"],
            "read_arr":          read_arr,
            "lcs_ref_positions": lcs_ref_positions,
            "orf_aa_row":        orf_aa_row,
            "ref_start":         ref_offset,
            "ref_end":           ref_end_pos,
        })

    if not read_data:
        return f"<p>Could not align reads for sample: {sample_name}</p>"

    vis_start = max(0,       min(rd["ref_start"] for rd in read_data) - 5)
    vis_end   = min(ref_len, max(rd["ref_end"]   for rd in read_data) + 5)
    lw = max(len(rd["label"]) for rd in read_data) + 6

    html_parts = [
        """<style>
        .gi-match    { color:#155724; font-weight:bold; background:#d4edda; }
        .gi-mismatch { color:#721c24; font-weight:bold; background:#f8d7da; }
        .gi-del      { color:#856404; background:#fff3cd; }
        .gi-absent   { color:#bbb; }
        .gi-lcs      { border-bottom:2px solid #e8a000; }
        .gi-ref      { color:#264653; }
        .gi-orf      { }
        .gi-orf-aa   { color:#2a6496; font-weight:bold; }
        </style>""",
        "<div style='font-family:monospace;font-size:13px;line-height:1.8;white-space:pre;overflow-x:auto;'>",
        f"<div style='font-weight:bold;margin-bottom:8px;white-space:normal;color:#264653;'>{sample_name}</div>",
    ]
    if lcs:
        html_parts.append(
            f"<div style='margin-bottom:10px;white-space:normal;color:#856404;'>"
            f"LCS overlap: {len(lcs)} bp — "
            f"<span style='border-bottom:2px solid #e8a000;padding-bottom:1px;"
            f"font-family:inherit;'>underlined</span></div>"
        )

    for block_start in range(vis_start, vis_end, block_size):
        block_end = min(block_start + block_size, vis_end)

        ref_label = f"Reference {block_start+1}-{block_end}:".ljust(lw)
        ref_html  = "".join(f"<span class='gi-ref'>{ref_dna[ri]}</span>"
                            for ri in range(block_start, block_end))

        def render_read_row(rd):
            q_label  = (rd["label"][:lw - 3] + ":").ljust(lw)
            arr      = rd["read_arr"]
            lcs_pos  = rd["lcs_ref_positions"]
            row_html = ""
            for ri in range(block_start, block_end):
                c      = arr[ri]
                ref_c  = ref_dna[ri]
                lcs_cl = " gi-lcs" if ri in lcs_pos else ""
                if c is None:
                    row_html += f"<span class='gi-absent{lcs_cl}'>·</span>"
                elif c == "-":
                    row_html += f"<span class='gi-del{lcs_cl}'>-</span>"
                elif c.upper() == ref_c.upper():
                    row_html += f"<span class='gi-match{lcs_cl}'>{c}</span>"
                else:
                    row_html += f"<span class='gi-mismatch{lcs_cl}'>{c}</span>"
            html_parts.append(f"<div>{q_label}{row_html}</div>")

        def render_orf_row(rd, orf_label):
            q_label  = (orf_label[:lw - 3] + ":").ljust(lw)
            aa_row   = rd["orf_aa_row"]
            row_html = ""
            for ri in range(block_start, block_end):
                ch = aa_row.get(ri)
                if ch is None:
                    row_html += "<span class='gi-absent'>·</span>"
                elif ch == '▓':
                    row_html += "<span class='gi-orf'> </span>"
                else:
                    row_html += f"<span class='gi-orf-aa'>{ch}</span>"
            html_parts.append(f"<div class='gi-orf-row'>{q_label}{row_html}</div>")

        html_parts.append("<div style='margin-bottom:12px;'>")
        # ORF row above forward, reference in the middle, reverse read(s) below with ORF rows
        if read_data:
            fwd = read_data[0]
            if fwd["orf_aa_row"]:
                render_orf_row(fwd, "ORF (fwd)")
            render_read_row(fwd)
        html_parts.append(f"<div style='font-weight:bold;'>{ref_label}{ref_html}</div>")
        for rd in read_data[1:]:
            render_read_row(rd)
            if rd["orf_aa_row"]:
                render_orf_row(rd, "ORF (rev)")

        html_parts.append("</div>")

    html_parts.append("</div>")
    return "".join(html_parts)


def get_longest_orf_six_frames(seq: str) -> Tuple[str, Optional[str]]:
    """
    Find the longest ORF across all six reading frames.

    Args:
        seq: DNA sequence string

    Returns:
        Tuple of (longest_orf_sequence, frame_name) where frame_name
        is one of +1, +2, +3, -1, -2, -3, or None if no ORF found
    """
    translations = translate_six_frames(seq)
    longest_orf = ""
    best_frame = None
    for frame, prot in translations.items():
        orf = longest_ORF(prot)
        if len(orf) > len(longest_orf):
            longest_orf = orf
            best_frame = frame
    return longest_orf, best_frame

def find_pairs(files: List[str]) -> List[Tuple[str, str, str]]:
    """
    Identify forward/reverse sequence pairs based on filename patterns.
    Forward = '-T7' in filename (before optional ' (n)' suffix and .seq)
    Reverse = '-T7-Term' in filename (before optional ' (n)' suffix and .seq)
    Returns list of tuples: (forward_file, reverse_file, prefix)
    """
    # Strip optional Windows duplicate-copy suffix like " (1)" before .seq
    def normalize(name):
        return re.sub(r'\s*\(\d+\)(?=\.seq$)', '', name, flags=re.IGNORECASE)

    # Map normalized-lowercase → original filename for reverse lookup
    norm_map = {normalize(f).lower(): f for f in files}

    pairs = []
    for f in files:
        norm_f = normalize(f)
        if re.search(r'-t7\.seq$', norm_f, re.IGNORECASE):
            prefix = re.sub(r'-t7\.seq$', '', norm_f, flags=re.IGNORECASE)
            rev_candidate = prefix + "-T7-Term.seq"
            rev_file = norm_map.get(rev_candidate.lower())
            if rev_file:
                pairs.append((f, rev_file, prefix))
    return pairs


def process_pair(forward_file: Union[str, Any], reverse_file: Union[str, Any],
                  output_prefix: str, min_lcs: int = 12,
                  ref_dict: Optional[Dict[str, str]] = None,
                  ref_dna_dict: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
    """
    Process a pair of forward/reverse DNA sequence files.

    Reads both files, merges them using LCS, then either:
    - If ref_dna_dict provided: Aligns merged DNA to DNA reference, then translates
    - Otherwise: Translates in six frames and finds longest ORF, compares to AA reference

    Args:
        forward_file: Path or file-like object for forward sequence
        reverse_file: Path or file-like object for reverse sequence
        output_prefix: Sample name/prefix for output
        min_lcs: Minimum LCS length required for merging (default: 12 bp)
        ref_dict: Optional dictionary of AA reference sequences {sample_name: sequence}
        ref_dna_dict: Optional dictionary of DNA reference sequences {sample_name: sequence}

    Returns:
        Dictionary containing:
            - merged_dna: Merged DNA sequence or None
            - longest_orf: Longest ORF amino acid sequence or None
            - orf_length: Length of longest ORF
            - frame: Reading frame of longest ORF (None if DNA alignment used)
            - merge_info: Details about the merge (LCS, trimmed lengths)
            - alignment_text: Alignment visualization or None
            - identity_percent: Percent identity to reference or 0.0
    """
    try:
        fwd_dna = read_seq_file(forward_file)
        rev_dna_raw = read_seq_file(reverse_file)
    except Exception as e:
        raise ValueError(f"Error reading sequence files: {e}")

    try:
        rev_dna = str(Seq(rev_dna_raw).reverse_complement())
        merge_info = merge_by_lcs_dna(fwd_dna, rev_dna, min_lcs=min_lcs)
        merged_dna = merge_info["merged"]
    except Exception as e:
        raise ValueError(f"Error processing sequences: {e}")

    sample_name = os.path.basename(output_prefix)
    sample_base = re.split(r'[_-]', sample_name, 1)[0].lower()

    # Build per-read inspections regardless of merge success
    _ref_dna_seq = ref_dna_dict.get(sample_base) if ref_dna_dict else None
    _ref_aa_seq = ref_dict.get(sample_base) if ref_dict else None
    _use_dna = bool(ref_dna_dict)
    raw_reads = [
        inspect_single_read(fwd_dna, f"{output_prefix} — T7 (forward)",
                            ref_dna=_ref_dna_seq, ref_aa=_ref_aa_seq,
                            use_dna_alignment=_use_dna),
        inspect_single_read(rev_dna, f"{output_prefix} — T7-Term (rev-comp)",
                            ref_dna=_ref_dna_seq, ref_aa=_ref_aa_seq,
                            use_dna_alignment=_use_dna),
    ]
    _lcs = merge_info.get("lcs") or None
    investigate_html = (
        build_investigate_html(raw_reads, _ref_dna_seq, lcs=_lcs, sample_name=sample_name)
        if _ref_dna_seq else None
    )

    if not merged_dna:
        return {
            "merged_dna": None,
            "longest_orf": None,
            "orf_dna": None,
            "orf_length": 0,
            "frame": None,
            "merge_info": merge_info,
            "alignment_text": None,
            "identity_percent": 0.0,
            "raw_reads": raw_reads,
            "investigate_html": investigate_html,
        }

    alignment_text = None
    identity_percent = 0.0
    longest_orf = None
    frame = None

    # DNA-to-DNA alignment mode: align DNA first, then translate
    final_orf_data = {"orf_count": 0, "alignment_with_orfs": None}
    if ref_dna_dict:
        ref_dna_seq = ref_dna_dict.get(sample_base)
        if ref_dna_seq:
            alignment_text, identity_percent, orf_raw, _, final_orf_data = align_dna_then_translate(
                merged_dna,
                ref_dna_seq,
                sample_name
            )
            longest_orf = longest_ORF(orf_raw)
            frame = "DNA-aligned"
        else:
            # No DNA reference found, fall back to six-frame translation
            longest_orf, frame = get_longest_orf_six_frames(merged_dna)
    else:
        # Original mode: translate first, then align AA to AA reference
        longest_orf, frame = get_longest_orf_six_frames(merged_dna)

        if ref_dict and longest_orf:
            ref_seq = ref_dict.get(sample_base)
            if ref_seq:
                alignment_text, identity_percent = visualize_orf_vs_reference(
                    longest_orf,
                    ref_seq,
                    sample_name
                )

    orf_dna = None
    if longest_orf and merged_dna:
        ds, de = _find_orf_dna_span(merged_dna, longest_orf)
        if ds is not None and de is not None:
            orf_dna = merged_dna[ds:de]

    return {
        "merged_dna": merged_dna,
        "longest_orf": longest_orf,
        "orf_dna": orf_dna,
        "orf_length": len(longest_orf) if longest_orf else 0,
        "frame": frame,
        "merge_info": merge_info,
        "alignment_text": alignment_text,
        "identity_percent": identity_percent,
        "orf_count": final_orf_data.get("orf_count", 0),
        "alignment_text_with_orfs": final_orf_data.get("alignment_with_orfs"),
        "raw_reads": raw_reads,
        "investigate_html": investigate_html,
    }



# ───────────────────────────────────────────────────────────────
#  Reference handling and alignment
# ───────────────────────────────────────────────────────────────

def load_reference_xlsx(filepath: Union[str, Any]) -> Dict[str, str]:
    """
    Load reference amino acid sequences from an Excel or CSV file.

    Reads a file with 'Clone number' and 'AA sequence' columns and
    creates a case-insensitive mapping.

    Args:
        filepath: Path to file or file-like object (.xlsx or .csv)

    Returns:
        Dictionary mapping lowercase reference names to amino acid sequences

    Raises:
        ValueError: If required columns are missing from the file
    """
    df = _read_ref_table(filepath)

    if df.empty:
        raise ValueError("Excel file is empty")

    df.columns = [c.strip().lower() for c in df.columns]
    if "clone number" not in df.columns or "aa sequence" not in df.columns:
        raise ValueError("Project file must have headers 'Clone number' and 'AA sequence' (case-insensitive).")

    def _clone_key(val):
        s = str(val).strip()
        if s.endswith(".0"):
            s = s[:-2]
        return s.lower()

    ref_dict = {
        _clone_key(row["clone number"]): str(row["aa sequence"]).strip()
        for _, row in df.iterrows()
        if pd.notna(row["clone number"]) and pd.notna(row["aa sequence"])
    }

    if not ref_dict:
        raise ValueError("No valid reference sequences found in Excel file")

    return ref_dict


def load_reference_dna_xlsx(filepath: Union[str, Any]) -> Dict[str, str]:
    """
    Load reference DNA sequences from an Excel or CSV file.

    Reads a file with 'Clone number' and 'DNA sequence' columns and
    creates a case-insensitive mapping.

    Args:
        filepath: Path to file or file-like object (.xlsx or .csv)

    Returns:
        Dictionary mapping lowercase reference names to DNA sequences

    Raises:
        ValueError: If required columns are missing from the file
    """
    df = _read_ref_table(filepath)

    if df.empty:
        raise ValueError("Excel file is empty")

    df.columns = [c.strip().lower() for c in df.columns]
    if "clone number" not in df.columns:
        raise ValueError("Project file must have 'Clone number' column (case-insensitive).")

    seq_col = None
    if "dna sequence" in df.columns:
        seq_col = "dna sequence"
    else:
        raise ValueError("Project file must have 'DNA sequence' column (case-insensitive).")

    def _clone_key(val):
        s = str(val).strip()
        if s.endswith(".0"):
            s = s[:-2]
        return s.lower()

    ref_dict = {
        _clone_key(row["clone number"]): str(row[seq_col]).strip().upper()
        for _, row in df.iterrows()
        if pd.notna(row["clone number"]) and pd.notna(row[seq_col])
    }

    if not ref_dict:
        raise ValueError("No valid DNA reference sequences found in Excel file")

    return ref_dict



def visualize_orf_vs_reference(orf_seq: str, ref_seq: str, sample_name: str,
                                block_size: int = 50,
                                query_circular: bool = False) -> Tuple[str, float]:
    """
    Create an HTML-formatted alignment between query ORF and reference sequence.

    Performs global or local pairwise alignment and generates a color-coded alignment
    visualization with position numbers. Matches are highlighted in green,
    mismatches in red, and gaps in gray.

    Args:
        orf_seq: Query ORF amino acid sequence
        ref_seq: Reference amino acid sequence
        sample_name: Name of the sample for display
        block_size: Number of amino acids per line in output (default: 50)
        query_circular: If True, use local alignment (better for circular plasmids)

    Returns:
        Tuple of (alignment_html, identity_percent) where alignment_html
        is an HTML formatted string and identity_percent is 0-100
    """
    original_query_len = len(orf_seq)

    # For circular queries, double the sequence to allow wraparound alignments
    query_for_align = orf_seq + orf_seq if query_circular else orf_seq

    # Step 1: Alignment using pairwise2
    # Use local alignment for circular queries (better for aligning reference to whole plasmid)
    # Use global alignment for linear sequences
    if query_circular:
        alignments = pairwise2.align.localms(query_for_align, ref_seq, 2, -1, -2, -0.5)
    else:
        alignments = pairwise2.align.globalms(query_for_align, ref_seq, 2, -1, -2, -0.5)
    if not alignments:
        return f"<p>No alignment found for sample: {sample_name}</p>", 0.0
    alignment = alignments[0]  # take the first/best alignment

    seqA_full = alignment[0]  # Query sequence with gaps
    seqB_full = alignment[1]  # Reference sequence with gaps

    # For local alignment, extract only the aligned region (begin to end)
    if query_circular:
        begin_pos = int(alignment[3])
        end_pos = int(alignment[4])
        seqA_str = seqA_full[begin_pos:end_pos]
        seqB_str = seqB_full[begin_pos:end_pos]
        # Track where in the query the alignment starts (for position numbering)
        query_start_offset = len(seqA_full[:begin_pos].replace('-', ''))
    else:
        seqA_str = seqA_full
        seqB_str = seqB_full
        query_start_offset = 0

    # Step 2: Compute identity
    matches = sum(a == b and a != '-' for a, b in zip(seqA_str, seqB_str))
    aligned_len = max(len(seqA_str), len(seqB_str))
    identity = (matches / aligned_len) * 100 if aligned_len else 0.0

    # Step 3: Build HTML output
    circular_note = " (circular query, local alignment)" if query_circular else ""

    html_parts = [
        """<style>
        .align-match { color: #155724; font-weight: bold; background-color: #d4edda; }
        .align-mismatch { color: #721c24; font-weight: bold; background-color: #f8d7da; }
        .align-gap { color: #856404; background-color: #fff3cd; }
        .align-match-pipe { color: #2a9d8f; }
        .ref-match { color: #155724; font-weight: bold; background-color: #d4edda; }
        .ref-mismatch { color: #721c24; background-color: #f8d7da; }
        </style>""",
        "<div style='font-family: monospace; font-size: 14px; line-height: 1.6; white-space: pre;'>",
        f"<div style='font-weight: bold; margin-bottom: 10px; white-space: normal; color: #264653;'>Alignment for sample: {sample_name}{circular_note}</div>",
        f"<div style='margin-bottom: 15px; white-space: normal; color: #264653;'>Alignment length: {aligned_len} | Identity: {identity:.1f}%</div>",
    ]

    # Show original length for circular query
    if query_circular:
        html_parts.append(f"<div style='margin-bottom: 10px; white-space: normal; color: #264653;'>Original query length: {original_query_len} aa | Alignment starts at position {(query_start_offset % original_query_len) + 1}</div>")
    q_counter = query_start_offset + 1  # Start from alignment position in query
    r_counter = 1
    label_width = 20

    for i in range(0, len(seqA_str), block_size):
        blockA = seqA_str[i:i+block_size]
        blockB = seqB_str[i:i+block_size]

        # Count non-gap characters to get end positions
        q_end = q_counter + len([c for c in blockA if c != '-']) - 1
        r_end = r_counter + len([c for c in blockB if c != '-']) - 1

        q_label = f"Query {q_counter}-{q_end}:".ljust(label_width)
        r_label = f"Ref   {r_counter}-{r_end}:".ljust(label_width)

        # Build color-coded sequences
        query_html = ""
        ref_html = ""
        match_html = ""

        for a, b in zip(blockA, blockB):
            # Query sequence coloring
            if a == '-':
                query_html += f"<span class='align-gap'>{a}</span>"
            elif a == b:
                query_html += f"<span class='align-match'>{a}</span>"  # Green for match
            else:
                query_html += f"<span class='align-mismatch'>{a}</span>"  # Red for mismatch

            # Reference sequence coloring
            if b == '-':
                ref_html += f"<span class='align-gap'>{b}</span>"
            elif a == b:
                ref_html += f"<span class='ref-match'>{b}</span>"  # Green for match
            else:
                ref_html += f"<span class='ref-mismatch'>{b}</span>"  # Red for mismatch

            # Match line
            if a == b and a != '-':
                match_html += "<span class='align-match-pipe'>|</span>"
            else:
                match_html += " "

        # Add block to HTML
        html_parts.append("<div style='margin-bottom: 12px;'>")
        html_parts.append(f"<div>{q_label}{query_html}</div>")
        html_parts.append(f"<div>{' ' * label_width}{match_html}</div>")
        html_parts.append(f"<div>{r_label}{ref_html}</div>")
        html_parts.append("</div>")

        # Update counters for next block
        q_counter += len([c for c in blockA if c != '-'])
        r_counter += len([c for c in blockB if c != '-'])

    html_parts.append("</div>")
    return "".join(html_parts), identity


def align_dna_then_translate(query_dna: str, ref_dna: str, sample_name: str,
                              block_size: int = 50,
                              query_circular: bool = False,
                              features: Optional[List[DetectedFeature]] = None) -> Tuple[str, float, str, str, dict]:
    """
    Align query DNA to reference DNA, then translate the aligned sequences.

    Performs DNA-to-DNA local alignment first to find the best matching region,
    then trims to only the region where the reference has sequence (no leading/
    trailing reference gaps). Only that region is scored and translated.

    Args:
        query_dna: Query DNA sequence
        ref_dna: Reference DNA sequence
        sample_name: Name of the sample for display
        block_size: Number of nucleotides per line in output (default: 50)
        query_circular: If True, double the query sequence to allow wraparound alignments
        features: Optional list of detected features to highlight in alignment

    Returns:
        Tuple of (alignment_html, identity_percent, aligned_query_aa, aligned_ref_aa, orf_data)
        where alignment_html is an HTML formatted string, identity_percent is 0-100,
        the aa sequences are the translated aligned regions, and orf_data is a dict with
        keys 'orf_count' and 'alignment_with_orfs'.
    """
    original_query_len = len(query_dna)

    # For circular queries, double the sequence to allow wraparound alignments
    query_for_align = query_dna + query_dna if query_circular else query_dna

    # Step 1: DNA alignment using local alignment to find best matching region
    alignments = pairwise2.align.localms(query_for_align, ref_dna, 2, -1, -2, -0.5)

    if not alignments:
        return f"<p>No alignment found for sample: {sample_name}</p>", 0.0, "", "", {"orf_count": 0, "alignment_with_orfs": None}

    alignment = alignments[0]  # take the first/best alignment

    seqA_full = alignment[0]  # Query DNA with gaps
    seqB_full = alignment[1]  # Reference DNA with gaps

    # Extract the aligned region from local alignment
    begin_pos = int(alignment[3])
    end_pos = int(alignment[4])
    seqA_aligned = seqA_full[begin_pos:end_pos]
    seqB_aligned = seqB_full[begin_pos:end_pos]

    # Step 2: Trim to only the region where reference has sequence (no leading/trailing ref gaps)
    # Find first and last position where reference is not a gap
    ref_start = 0
    ref_end = len(seqB_aligned)

    for i, c in enumerate(seqB_aligned):
        if c != '-':
            ref_start = i
            break

    for i in range(len(seqB_aligned) - 1, -1, -1):
        if seqB_aligned[i] != '-':
            ref_end = i + 1
            break

    # Extract only the reference-covered region
    seqA_str = seqA_aligned[ref_start:ref_end]
    seqB_str = seqB_aligned[ref_start:ref_end]

    # Calculate query start offset (for position numbering)
    query_start_offset = len(seqA_full[:begin_pos].replace('-', '')) + len(seqA_aligned[:ref_start].replace('-', ''))

    # Step 3: Compute DNA identity only on reference-covered region
    matches = sum(a == b and a != '-' for a, b in zip(seqA_str, seqB_str))
    # Use reference length (without gaps) as denominator for identity
    ref_len_no_gaps = len(seqB_str.replace('-', ''))
    identity = (matches / ref_len_no_gaps) * 100 if ref_len_no_gaps else 0.0

    # Step 4: Translate only the reference-covered aligned DNA sequences
    # Try all 3 reading frame offsets and pick the one with the longest ORF
    query_no_gaps = seqA_str.replace('-', '')
    ref_no_gaps = seqB_str.replace('-', '')

    def _best_frame(dna):
        best_aa, best_off, best_len = "", 0, 0
        for off in range(3):
            sub = dna[off:]
            tl = len(sub) - len(sub) % 3
            if tl == 0:
                continue
            aa = str(Seq(sub[:tl]).translate())
            orf = longest_ORF(aa)
            if len(orf) > best_len:
                best_len = len(orf)
                best_aa = aa
                best_off = off
        return best_aa, best_off

    aligned_query_aa, query_frame_offset = _best_frame(query_no_gaps)
    aligned_ref_aa, ref_frame_offset = _best_frame(ref_no_gaps)

    # Compute ORFs > 70 aa in each translated sequence
    query_orfs = find_all_orfs(aligned_query_aa)
    ref_orfs   = find_all_orfs(aligned_ref_aa)
    total_orf_count = len(query_orfs) + len(ref_orfs)
    best_query_orf = longest_ORF(aligned_query_aa)

    # Precompute ungapped→gapped position maps for ORF row generation
    u2g_query = build_ungapped_to_gapped_map(seqA_str)
    u2g_ref   = build_ungapped_to_gapped_map(seqB_str)

    # Map features to gapped alignment positions if provided
    feature_map = {}
    if features:
        mapped_features = map_features_to_gapped_alignment(features, seqA_str)
        feature_map = build_feature_position_map(mapped_features, use_gapped=True)

    # Step 5: Build HTML output
    circular_note = " (circular query)" if query_circular else ""
    aligned_region_len = len(seqA_str)

    html_parts = [
        """<style>
        .align-match { color: #155724; font-weight: bold; background-color: #d4edda; }
        .align-mismatch { color: #721c24; font-weight: bold; background-color: #f8d7da; }
        .align-gap { color: #856404; background-color: #fff3cd; }
        .align-match-pipe { color: #2a9d8f; }
        .ref-match { color: #155724; font-weight: bold; background-color: #d4edda; }
        .ref-mismatch { color: #721c24; background-color: #f8d7da; }
        .dna-match { color: #155724; font-weight: bold; background-color: #d4edda; }
        .dna-mismatch { color: #721c24; background-color: #f8d7da; }
        .dna-gap { color: #856404; background-color: #fff3cd; }
        </style>""",
        "<div style='font-family: monospace; font-size: 14px; line-height: 1.6; white-space: pre;'>",
        f"<div style='font-weight: bold; margin-bottom: 10px; white-space: normal; color: #264653;'>DNA Alignment for sample: {sample_name}{circular_note}</div>",
        f"<div style='margin-bottom: 10px; white-space: normal; color: #264653;'>Reference-aligned region: {aligned_region_len} bp | DNA Identity: {identity:.1f}% (vs {ref_len_no_gaps} bp reference)</div>",
    ]

    # Show translation info
    html_parts.append(f"<div style='margin-bottom: 10px; white-space: normal; color: #2a9d8f;'>Translated query: {len(aligned_query_aa)} aa | Translated reference: {len(aligned_ref_aa)} aa</div>")

    if query_circular:
        start_pos = (query_start_offset % original_query_len) + 1 if original_query_len > 0 else query_start_offset + 1
        html_parts.append(f"<div style='margin-bottom: 10px; white-space: normal; color: #264653;'>Original query length: {original_query_len} bp | Alignment starts at query position {start_pos}</div>")
    else:
        html_parts.append(f"<div style='margin-bottom: 10px; white-space: normal; color: #264653;'>Alignment starts at query position {query_start_offset + 1}</div>")

    q_counter = query_start_offset + 1
    r_counter = 1
    label_width = 22

    def orf_row_for_block(orfs, u2g_map, block_start, bsize, lw, label_prefix, frame_offset=0):
        rows = []
        for orf_idx, (orf_aa, orf_start_aa) in enumerate(orfs):
            chars = [' '] * bsize
            for aa_idx, aa_char in enumerate(orf_aa):
                mid_ungapped = frame_offset + orf_start_aa * 3 + aa_idx * 3 + 1  # middle nt of codon
                if mid_ungapped >= len(u2g_map):
                    break
                gapped_pos = u2g_map[mid_ungapped]
                col = gapped_pos - block_start
                if 0 <= col < bsize:
                    chars[col] = aa_char
            label = f"{label_prefix}{orf_idx+1}:".ljust(lw)
            styled = f"<span style='color:#2a6496;font-weight:bold'>{''.join(chars)}</span>"
            rows.append(f"<div>{label}{styled}</div>")
        return rows

    # orf_parts builds a parallel version of the alignment with ORF rows interspersed
    orf_parts = list(html_parts) if total_orf_count > 0 else []

    for i in range(0, len(seqA_str), block_size):
        blockA = seqA_str[i:i+block_size]
        blockB = seqB_str[i:i+block_size]

        # Count non-gap characters to get end positions
        q_end = q_counter + len([c for c in blockA if c != '-']) - 1
        r_end = r_counter + len([c for c in blockB if c != '-']) - 1

        q_label = f"Query  {q_counter}-{q_end}:".ljust(label_width)
        r_label = f"Ref    {r_counter}-{r_end}:".ljust(label_width)

        # Build color-coded DNA sequences
        query_html = ""
        ref_html = ""
        match_html = ""

        for j, (a, b) in enumerate(zip(blockA, blockB)):
            gapped_pos = i + j

            # Check for feature at this position
            feature = feature_map.get(gapped_pos)
            if feature:
                bg_style = f"background-color: {feature.color}40;"  # 40 = 25% opacity
                title_attr = f" title='{feature.name}'"
            else:
                bg_style = ""
                title_attr = ""

            # Query sequence coloring
            if a == '-':
                query_html += f"<span class='dna-gap' style='{bg_style}'{title_attr}>{a}</span>"
            elif a == b:
                query_html += f"<span class='dna-match' style='{bg_style}'{title_attr}>{a}</span>"
            else:
                query_html += f"<span class='dna-mismatch' style='{bg_style}'{title_attr}>{a}</span>"

            # Reference sequence coloring (no feature highlighting on reference)
            if b == '-':
                ref_html += f"<span class='dna-gap'>{b}</span>"
            elif a == b:
                ref_html += f"<span class='dna-match'>{b}</span>"
            else:
                ref_html += f"<span class='dna-mismatch'>{b}</span>"

            # Match line
            if a == b and a != '-':
                match_html += "<span class='align-match-pipe'>|</span>"
            else:
                match_html += " "

        # Add block to HTML
        html_parts.append("<div style='margin-bottom: 12px;'>")
        html_parts.append(f"<div>{q_label}{query_html}</div>")
        html_parts.append(f"<div>{' ' * label_width}{match_html}</div>")
        html_parts.append(f"<div>{r_label}{ref_html}</div>")
        html_parts.append("</div>")

        # Add ORF-annotated block to orf_parts
        if total_orf_count > 0:
            bsize = len(blockA)
            orf_parts.append("<div style='margin-bottom: 12px;'>")
            orf_parts.extend(orf_row_for_block(query_orfs, u2g_query, i, bsize, label_width, "QryORF", query_frame_offset))
            orf_parts.append(f"<div>{q_label}{query_html}</div>")
            orf_parts.append(f"<div>{' ' * label_width}{match_html}</div>")
            orf_parts.append(f"<div>{r_label}{ref_html}</div>")
            orf_parts.extend(orf_row_for_block(ref_orfs, u2g_ref, i, bsize, label_width, "RefORF", ref_frame_offset))
            orf_parts.append("</div>")

        q_counter += len([c for c in blockA if c != '-'])
        r_counter += len([c for c in blockB if c != '-'])

    # Step 5: Add protein translation alignment section (only when ORF > 70 aa)
    if best_query_orf and len(best_query_orf) > 70:
        html_parts.append("<hr style='margin: 20px 0; border: 1px solid #ccc;'>")
        html_parts.append("<div style='font-weight: bold; margin-bottom: 10px; white-space: normal;'>Translated Protein Alignment:</div>")
        if total_orf_count > 0:
            orf_parts.append("<hr style='margin: 20px 0; border: 1px solid #ccc;'>")
            orf_parts.append("<div style='font-weight: bold; margin-bottom: 10px; white-space: normal;'>Translated Protein Alignment:</div>")

    # Do a protein alignment of the translated sequences (gated on ORF > 70 aa)
    if best_query_orf and len(best_query_orf) > 70 and aligned_query_aa and aligned_ref_aa:
        prot_alignments = pairwise2.align.globalms(aligned_query_aa, aligned_ref_aa, 2, -1, -2, -0.5)
        if prot_alignments:
            prot_align = prot_alignments[0]
            prot_A = prot_align[0]
            prot_B = prot_align[1]

            prot_matches = sum(a == b and a != '-' for a, b in zip(prot_A, prot_B))
            prot_len = max(len(prot_A), len(prot_B))
            prot_identity = (prot_matches / prot_len) * 100 if prot_len else 0.0

            html_parts.append(f"<div style='margin-bottom: 10px; white-space: normal;'>Protein Identity: {prot_identity:.1f}%</div>")
            if total_orf_count > 0:
                orf_parts.append(f"<div style='margin-bottom: 10px; white-space: normal;'>Protein Identity: {prot_identity:.1f}%</div>")

            pq_counter = 1
            pr_counter = 1

            for i in range(0, len(prot_A), block_size):
                blockA = prot_A[i:i+block_size]
                blockB = prot_B[i:i+block_size]

                pq_end = pq_counter + len([c for c in blockA if c != '-']) - 1
                pr_end = pr_counter + len([c for c in blockB if c != '-']) - 1

                pq_label = f"Query  {pq_counter}-{pq_end}:".ljust(label_width)
                pr_label = f"Ref    {pr_counter}-{pr_end}:".ljust(label_width)

                pquery_html = ""
                pref_html = ""
                pmatch_html = ""

                for a, b in zip(blockA, blockB):
                    if a == '-':
                        pquery_html += f"<span class='align-gap'>{a}</span>"
                    elif a == b:
                        pquery_html += f"<span class='ref-match'>{a}</span>"
                    else:
                        pquery_html += f"<span class='ref-mismatch'>{a}</span>"

                    if b == '-':
                        pref_html += f"<span class='align-gap'>{b}</span>"
                    elif a == b:
                        pref_html += f"<span class='ref-match'>{b}</span>"
                    else:
                        pref_html += f"<span class='ref-mismatch'>{b}</span>"

                    if a == b and a != '-':
                        pmatch_html += "<span class='align-match-pipe'>|</span>"
                    else:
                        pmatch_html += " "

                html_parts.append("<div style='margin-bottom: 12px;'>")
                html_parts.append(f"<div>{pq_label}{pquery_html}</div>")
                html_parts.append(f"<div>{' ' * label_width}{pmatch_html}</div>")
                html_parts.append(f"<div>{pr_label}{pref_html}</div>")
                html_parts.append("</div>")
                if total_orf_count > 0:
                    orf_parts.append("<div style='margin-bottom: 12px;'>")
                    orf_parts.append(f"<div>{pq_label}{pquery_html}</div>")
                    orf_parts.append(f"<div>{' ' * label_width}{pmatch_html}</div>")
                    orf_parts.append(f"<div>{pr_label}{pref_html}</div>")
                    orf_parts.append("</div>")

                pq_counter += len([c for c in blockA if c != '-'])
                pr_counter += len([c for c in blockB if c != '-'])

    html_parts.append("</div>")
    if total_orf_count > 0:
        orf_parts.append("</div>")
    orf_data = {
        "orf_count": total_orf_count,
        "alignment_with_orfs": "".join(orf_parts) if total_orf_count > 0 else None,
    }
    return "".join(html_parts), identity, aligned_query_aa, aligned_ref_aa, orf_data


# ───────────────────────────────────────────────────────────────
#  Main processing
# ───────────────────────────────────────────────────────────────

def main(seq_files: List[Any], ref_file: Union[str, Any], min_lcs: int = 12,
         use_dna_alignment: bool = True) -> List[Dict[str, Any]]:
    """
    Main processing function for batch sequence merging and analysis.

    Identifies forward/reverse pairs from uploaded files, merges them,
    then either aligns DNA to DNA reference and translates (default),
    or translates to amino acids first and compares to AA reference.

    Args:
        seq_files: List of uploaded sequence files (file-like objects with .name attribute)
        ref_file: Excel file with reference sequences (file-like object or path)
        min_lcs: Minimum LCS length for merging (default: 12 bp)
        use_dna_alignment: If True (default), align DNA to DNA reference then translate.
                          If False, translate first then align AA to AA reference.

    Returns:
        List of result dictionaries, one per processed pair, each containing:
            - sample: Sample name/prefix
            - merged_dna: Merged DNA sequence or None
            - longest_orf: Longest ORF amino acid sequence or None
            - orf_length: Length of ORF
            - frame: Reading frame or "DNA-aligned"
            - merge_info: Merge details (LCS info)
            - reference: Reference sequence or None
            - alignment_text: Alignment visualization or None
            - identity_percent: Percent identity to reference
    """
    # Validate inputs
    if not seq_files:
        raise ValueError("No sequence files provided")

    # Load reference - DNA for DNA alignment mode, AA otherwise
    if use_dna_alignment:
        ref_dna_dict = load_reference_dna_xlsx(ref_file)
        ref_dict = None
    else:
        ref_dict = load_reference_xlsx(ref_file)
        ref_dna_dict = None

    # Find pairs
    files = [f.name for f in seq_files]  # names for pairing
    pairs = find_pairs(files)

    if not pairs:
        raise ValueError("No valid forward/reverse pairs found. Files should follow naming pattern: *-T7.seq and *-T7-Term.seq")

    all_results = []

    for fwd_name, rev_name, prefix in pairs:
        # Map names back to uploaded files
        fwd_file = next(f for f in seq_files if f.name == fwd_name)
        rev_file = next(f for f in seq_files if f.name == rev_name)

        result = process_pair(fwd_file, rev_file, prefix, min_lcs=min_lcs,
                              ref_dict=ref_dict, ref_dna_dict=ref_dna_dict)
        result["sample"] = prefix

        sample_base = re.split(r'[_-]', prefix, 1)[0].lower()
        if ref_dna_dict:
            result["reference"] = ref_dna_dict.get(sample_base)
        elif ref_dict:
            result["reference"] = ref_dict.get(sample_base)
        all_results.append(result)

    return all_results


def process_circular_queries(seq_files: List[Any], ref_file: Union[str, Any],
                              use_dna_alignment: bool = True,
                              selected_categories: Optional[List[str]] = None) -> List[Dict[str, Any]]:
    """
    Process circular query sequences (whole plasmids in FASTA format).

    Skips pair-finding and merging. When use_dna_alignment is True (default),
    aligns query DNA to DNA reference first, then translates the aligned sequence.
    Matches query files to references by filename prefix (text before first underscore).
    If initial alignment has poor identity, tries reverse complement of DNA.

    Args:
        seq_files: List of uploaded FASTA files (file-like objects with .name attribute)
        ref_file: Excel file with reference sequences (file-like object or path)
        use_dna_alignment: If True (default), align DNA to DNA reference then translate.
                          If False, use original AA-to-AA alignment mode.
        selected_categories: Optional list of feature categories to detect

    Returns:
        List of result dictionaries, one per processed file, each containing:
            - sample: Filename
            - merged_dna: The query DNA sequence (in best orientation)
            - longest_orf: Longest ORF amino acid sequence
            - orf_length: Length of ORF
            - frame: Reading frame or "DNA-aligned"
            - merge_info: Empty dict
            - reference: Reference sequence or None
            - alignment_text: Alignment visualization or None
            - identity_percent: Percent identity to reference
            - orientation: "forward" or "reverse complement"
            - detected_features: List of detected DNA features
    """
    if not seq_files:
        raise ValueError("No sequence files provided")

    # Load reference sequences - DNA for DNA alignment mode, AA otherwise
    if use_dna_alignment:
        ref_dict = load_reference_dna_xlsx(ref_file)
    else:
        ref_dict = load_reference_xlsx(ref_file)

    results = []
    for seq_file in seq_files:
        filename = seq_file.name

        try:
            # Read the FASTA sequence
            query_dna = read_seq_file(seq_file)

            # Match by filename prefix - try multiple strategies
            base_name = filename.rsplit(".", 1)[0].lower()  # Remove extension
            prefix_underscore = re.split(r'[_-]', filename, 1)[0].lower()  # Before first _ or -

            # Try multiple matching strategies in order of specificity
            ref_seq = None
            matched_ref_name = None
            candidates = [
                prefix_underscore,                          # e.g., "11421984-1"
                base_name,                                  # e.g., "11421984-1_contig"
                prefix_underscore.replace("-", ""),         # e.g., "114219841" (no hyphen)
                prefix_underscore.split("-")[0],            # e.g., "11421984" (before hyphen)
            ]

            for candidate in candidates:
                if candidate in ref_dict:
                    ref_seq = ref_dict[candidate]
                    matched_ref_name = candidate
                    break

            # If still no match, try partial matching
            if not ref_seq:
                for ref_name, ref_sequence in ref_dict.items():
                    if any(c in ref_name or ref_name in c for c in candidates if c):
                        ref_seq = ref_sequence
                        matched_ref_name = ref_name
                        break

            # Initialize results
            alignment_text = None
            identity_percent = 0.0
            orientation = "forward"
            final_dna = query_dna
            final_orf = None
            final_frame = None
            detected_features = []
            final_orf_data = {"orf_count": 0, "alignment_with_orfs": None}
            raw_reads = []
            investigate_html = None

            # Detect features on both orientations
            features_fwd = detect_features(query_dna, selected_categories=selected_categories)
            query_dna_rc = str(Seq(query_dna).reverse_complement())
            features_rc = detect_features(query_dna_rc, selected_categories=selected_categories)

            if use_dna_alignment and ref_seq:
                # DNA-to-DNA alignment mode: align DNA first, then translate
                # Always try both orientations and use the better alignment
                alignment_text_fwd, identity_fwd, orf_fwd_raw, _, orf_data_fwd = align_dna_then_translate(
                    query_dna, ref_seq, filename,
                    query_circular=True,
                    features=features_fwd
                )
                orf_fwd = longest_ORF(orf_fwd_raw)

                alignment_text_rc, identity_rc, orf_rc_raw, _, orf_data_rc = align_dna_then_translate(
                    query_dna_rc, ref_seq, f"{filename} (reverse complement)",
                    query_circular=True,
                    features=features_rc
                )
                orf_rc = longest_ORF(orf_rc_raw)

                # Use the orientation with better identity
                if identity_rc > identity_fwd:
                    alignment_text = alignment_text_rc
                    identity_percent = identity_rc
                    orientation = "reverse complement"
                    final_dna = query_dna_rc
                    final_orf = orf_rc
                    final_orf_data = orf_data_rc
                    detected_features = features_rc
                else:
                    alignment_text = alignment_text_fwd
                    identity_percent = identity_fwd
                    final_orf = orf_fwd
                    final_orf_data = orf_data_fwd
                    detected_features = features_fwd

                final_frame = "DNA-aligned"

                _o_fwd = find_all_orfs(orf_fwd_raw)
                _o_rc  = find_all_orfs(orf_rc_raw)
                raw_reads = [
                    {
                        "label": f"{filename} (forward)",
                        "dna": query_dna, "length": len(query_dna),
                        "alignment_html": alignment_text_fwd,
                        "alignment_html_with_orfs": orf_data_fwd.get("alignment_with_orfs"),
                        "identity": identity_fwd,
                        "orf_count": orf_data_fwd.get("orf_count", 0),
                        "best_orf_length": max((len(o) for o, _ in _o_fwd), default=0),
                        "all_orfs": [{"sequence": o, "length": len(o), "start": s} for o, s in _o_fwd],
                    },
                    {
                        "label": f"{filename} (reverse complement)",
                        "dna": query_dna_rc, "length": len(query_dna_rc),
                        "alignment_html": alignment_text_rc,
                        "alignment_html_with_orfs": orf_data_rc.get("alignment_with_orfs"),
                        "identity": identity_rc,
                        "orf_count": orf_data_rc.get("orf_count", 0),
                        "best_orf_length": max((len(o) for o, _ in _o_rc), default=0),
                        "all_orfs": [{"sequence": o, "length": len(o), "start": s} for o, s in _o_rc],
                    },
                ]
                investigate_html = build_investigate_html(
                    raw_reads, ref_seq, lcs=None, sample_name=filename
                )
            else:
                # Original mode: translate first, find ORF, then align AA to AA
                # Always try both orientations and use the better alignment
                longest_orf_fwd, frame_fwd = get_longest_orf_six_frames(query_dna)
                longest_orf_rc, frame_rc = get_longest_orf_six_frames(query_dna_rc)

                if ref_seq:
                    # Try forward orientation
                    identity_fwd = 0.0
                    alignment_text_fwd = None
                    if longest_orf_fwd:
                        alignment_text_fwd, identity_fwd = visualize_orf_vs_reference(
                            longest_orf_fwd, ref_seq, filename,
                            query_circular=True
                        )

                    # Try reverse complement orientation
                    identity_rc = 0.0
                    alignment_text_rc = None
                    if longest_orf_rc:
                        alignment_text_rc, identity_rc = visualize_orf_vs_reference(
                            longest_orf_rc, ref_seq, f"{filename} (reverse complement)",
                            query_circular=True
                        )

                    # Use the orientation with better identity
                    if identity_rc > identity_fwd:
                        alignment_text = alignment_text_rc
                        identity_percent = identity_rc
                        orientation = "reverse complement"
                        final_dna = query_dna_rc
                        final_orf = longest_orf_rc
                        final_frame = frame_rc
                        detected_features = features_rc
                    else:
                        alignment_text = alignment_text_fwd
                        identity_percent = identity_fwd
                        final_orf = longest_orf_fwd
                        final_frame = frame_fwd
                        detected_features = features_fwd

                    _o_fwd = find_all_orfs(longest_orf_fwd or "")
                    _o_rc  = find_all_orfs(longest_orf_rc or "")
                    raw_reads = [
                        {
                            "label": f"{filename} (forward)",
                            "dna": query_dna, "length": len(query_dna),
                            "alignment_html": alignment_text_fwd,
                            "alignment_html_with_orfs": None,
                            "identity": identity_fwd,
                            "orf_count": len(_o_fwd),
                            "best_orf_length": max((len(o) for o, _ in _o_fwd), default=0),
                            "all_orfs": [{"sequence": o, "length": len(o), "start": s} for o, s in _o_fwd],
                        },
                        {
                            "label": f"{filename} (reverse complement)",
                            "dna": query_dna_rc, "length": len(query_dna_rc),
                            "alignment_html": alignment_text_rc,
                            "alignment_html_with_orfs": None,
                            "identity": identity_rc,
                            "orf_count": len(_o_rc),
                            "best_orf_length": max((len(o) for o, _ in _o_rc), default=0),
                            "all_orfs": [{"sequence": o, "length": len(o), "start": s} for o, s in _o_rc],
                        },
                    ]
                    investigate_html = build_investigate_html(
                        raw_reads, ref_seq, lcs=None, sample_name=filename
                    )
                else:
                    # No reference - just use forward orientation
                    final_orf = longest_orf_fwd
                    final_frame = frame_fwd
                    detected_features = features_fwd
                    _o_fwd = find_all_orfs(longest_orf_fwd or "")
                    raw_reads = [
                        {
                            "label": f"{filename} (forward)",
                            "dna": query_dna, "length": len(query_dna),
                            "alignment_html": None, "alignment_html_with_orfs": None,
                            "identity": 0.0,
                            "orf_count": len(_o_fwd),
                            "best_orf_length": max((len(o) for o, _ in _o_fwd), default=0),
                            "all_orfs": [{"sequence": o, "length": len(o), "start": s} for o, s in _o_fwd],
                        },
                    ]

            _orf_dna = None
            if final_orf and final_dna:
                _ds, _de = _find_orf_dna_span(final_dna, final_orf)
                if _ds is not None and _de is not None:
                    _orf_dna = final_dna[_ds:_de]

            results.append({
                "sample": filename,
                "merged_dna": final_dna,
                "longest_orf": final_orf,
                "orf_dna": _orf_dna,
                "orf_length": len(final_orf) if final_orf else 0,
                "frame": final_frame,
                "merge_info": {},
                "reference": ref_seq,
                "matched_ref_name": matched_ref_name,
                "orientation": orientation,
                "available_refs": list(ref_dict.keys()) if not ref_seq else None,
                "tried_prefixes": candidates if not ref_seq else None,
                "alignment_text": alignment_text,
                "identity_percent": identity_percent,
                "detected_features": detected_features,
                "orf_count": final_orf_data.get("orf_count", 0),
                "alignment_text_with_orfs": final_orf_data.get("alignment_with_orfs"),
                "raw_reads": raw_reads,
                "investigate_html": investigate_html,
            })

        except Exception as e:
            results.append({
                "sample": filename,
                "error": str(e),
                "merged_dna": None,
                "longest_orf": None,
                "orf_dna": None,
                "orf_length": 0,
                "frame": None,
                "merge_info": {},
                "reference": None,
                "alignment_text": None,
                "identity_percent": 0.0,
                "detected_features": [],
                "orf_count": 0,
                "alignment_text_with_orfs": None,
                "raw_reads": [],
                "investigate_html": None,
            })

    return results


ORF_MIN_AA = 30  # Minimum ORF length (amino acids) to be considered significant


def process_unidirectional(seq_files: List[Any], ref_file: Union[str, Any],
                            use_dna_alignment: bool = True,
                            selected_categories: Optional[List[str]] = None) -> List[Dict[str, Any]]:
    """
    Process unidirectional sequences (single reads like T7-only or M13R-only).

    Uses local alignment since unidirectional sequences are unlikely to span
    the entire reference sequence. Tries both forward and reverse complement
    orientations and uses the better alignment.

    Args:
        seq_files: List of uploaded sequence files (file-like objects with .name attribute)
        ref_file: Excel file with reference sequences (file-like object or path)
        use_dna_alignment: If True (default), align DNA to DNA reference then translate.
                          If False, translate first then align AA to AA reference.
        selected_categories: Optional list of feature categories to detect

    Returns:
        List of result dictionaries, one per processed file, each containing:
            - sample: Filename
            - merged_dna: The query DNA sequence (in best orientation)
            - longest_orf: Translated amino acid sequence from aligned region
            - orf_length: Length of translated sequence
            - frame: "DNA-aligned" or reading frame
            - merge_info: Empty dict
            - reference: Reference sequence or None
            - alignment_text: Alignment visualization or None
            - identity_percent: Percent identity to reference
            - orientation: "forward" or "reverse complement"
            - detected_features: List of detected DNA features
    """
    if not seq_files:
        raise ValueError("No sequence files provided")

    # Load reference sequences - DNA for DNA alignment mode, AA otherwise
    if use_dna_alignment:
        ref_dict = load_reference_dna_xlsx(ref_file)
    else:
        ref_dict = load_reference_xlsx(ref_file)

    results = []
    for seq_file in seq_files:
        filename = seq_file.name

        try:
            # Read the sequence
            query_dna = read_seq_file(seq_file)

            # Check if this is a T7-term read (sequences in reverse direction)
            # These need to be reverse-complemented before alignment
            is_reverse_primer = bool(re.search(r't7[-_]?term', filename, re.IGNORECASE))
            if is_reverse_primer:
                query_dna = str(Seq(query_dna).reverse_complement())

            # Match by filename prefix - try multiple strategies
            # Remove extension first
            base_name = filename.rsplit(".", 1)[0].lower()
            base_name = re.sub(r'\s*\(\d+\)$', '', base_name)  # strip Windows dup suffix e.g. " (1)"
            # Remove common suffixes like -T7, -T7-Term, -M13R, etc.
            clean_name = re.sub(r'[-_]?(t7[-_]?term|t7|m13r|m13f)$', '', base_name, flags=re.IGNORECASE)
            prefix_underscore = re.split(r'[_-]', clean_name, 1)[0].lower()

            # Try multiple matching strategies in order of specificity
            ref_seq = None
            matched_ref_name = None
            candidates = [
                clean_name,                                 # e.g., "sample1"
                prefix_underscore,                          # e.g., "sample1" (before underscore)
                prefix_underscore.replace("-", ""),         # e.g., "sample1" (no hyphen)
                prefix_underscore.split("-")[0],            # e.g., "sample" (before hyphen)
            ]

            for candidate in candidates:
                if candidate in ref_dict:
                    ref_seq = ref_dict[candidate]
                    matched_ref_name = candidate
                    break

            # If still no match, try partial matching
            if not ref_seq:
                for ref_name, ref_sequence in ref_dict.items():
                    if any(c in ref_name or ref_name in c for c in candidates if c):
                        ref_seq = ref_sequence
                        matched_ref_name = ref_name
                        break

            # Initialize results
            alignment_text = None
            identity_percent = 0.0
            # Track orientation - note if auto-flipped due to primer type
            orientation = "T7-term (auto-flipped)" if is_reverse_primer else "forward"
            final_dna = query_dna
            final_orf = None
            final_frame = None
            detected_features = []
            final_orf_data = {"orf_count": 0, "alignment_with_orfs": None}
            raw_reads = []
            investigate_html = None

            # Detect features on both orientations
            features_fwd = detect_features(query_dna, selected_categories=selected_categories)
            query_dna_rc = str(Seq(query_dna).reverse_complement())
            features_rc = detect_features(query_dna_rc, selected_categories=selected_categories)

            if use_dna_alignment and ref_seq:
                # DNA-to-DNA local alignment mode
                # Always try both orientations and use the better alignment
                alignment_text_fwd, identity_fwd, orf_fwd_raw, _, orf_data_fwd = align_dna_then_translate(
                    query_dna, ref_seq, filename,
                    query_circular=False,
                    features=features_fwd
                )
                orf_fwd = longest_ORF(orf_fwd_raw)

                alignment_text_rc, identity_rc, orf_rc_raw, _, orf_data_rc = align_dna_then_translate(
                    query_dna_rc, ref_seq, f"{filename} (reverse complement)",
                    query_circular=False,
                    features=features_rc
                )
                orf_rc = longest_ORF(orf_rc_raw)

                # Use the orientation with better identity
                if identity_rc > identity_fwd:
                    alignment_text = alignment_text_rc
                    identity_percent = identity_rc
                    # Update orientation label based on original primer type
                    orientation = "forward (flipped back)" if is_reverse_primer else "reverse complement"
                    final_dna = query_dna_rc
                    final_orf = orf_rc
                    final_orf_data = orf_data_rc
                    detected_features = features_rc
                else:
                    alignment_text = alignment_text_fwd
                    identity_percent = identity_fwd
                    final_orf = orf_fwd
                    final_orf_data = orf_data_fwd
                    detected_features = features_fwd

                final_frame = "DNA-aligned"

                _o_fwd = find_all_orfs(orf_fwd_raw)
                _o_rc  = find_all_orfs(orf_rc_raw)
                raw_reads = [
                    {
                        "label": f"{filename} (forward)",
                        "dna": query_dna, "length": len(query_dna),
                        "alignment_html": alignment_text_fwd,
                        "alignment_html_with_orfs": orf_data_fwd.get("alignment_with_orfs"),
                        "identity": identity_fwd,
                        "orf_count": orf_data_fwd.get("orf_count", 0),
                        "best_orf_length": max((len(o) for o, _ in _o_fwd), default=0),
                        "all_orfs": [{"sequence": o, "length": len(o), "start": s} for o, s in _o_fwd],
                    },
                    {
                        "label": f"{filename} (reverse complement)",
                        "dna": query_dna_rc, "length": len(query_dna_rc),
                        "alignment_html": alignment_text_rc,
                        "alignment_html_with_orfs": orf_data_rc.get("alignment_with_orfs"),
                        "identity": identity_rc,
                        "orf_count": orf_data_rc.get("orf_count", 0),
                        "best_orf_length": max((len(o) for o, _ in _o_rc), default=0),
                        "all_orfs": [{"sequence": o, "length": len(o), "start": s} for o, s in _o_rc],
                    },
                ]
                investigate_html = build_investigate_html(
                    raw_reads, ref_seq, lcs=None, sample_name=filename
                )

            elif ref_seq:
                # AA alignment mode: translate first, then align
                # Always try both orientations and use the better alignment
                longest_orf_fwd, frame_fwd = get_longest_orf_six_frames(query_dna)
                longest_orf_rc, frame_rc = get_longest_orf_six_frames(query_dna_rc)

                # Try forward orientation
                identity_fwd = 0.0
                alignment_text_fwd = None
                if longest_orf_fwd:
                    alignment_text_fwd, identity_fwd = visualize_orf_vs_reference(
                        longest_orf_fwd, ref_seq, filename,
                        query_circular=True
                    )

                # Try reverse complement orientation
                identity_rc = 0.0
                alignment_text_rc = None
                if longest_orf_rc:
                    alignment_text_rc, identity_rc = visualize_orf_vs_reference(
                        longest_orf_rc, ref_seq, f"{filename} (reverse complement)",
                        query_circular=True
                    )

                # Use the orientation with better identity
                if identity_rc > identity_fwd:
                    alignment_text = alignment_text_rc
                    identity_percent = identity_rc
                    orientation = "forward (flipped back)" if is_reverse_primer else "reverse complement"
                    final_dna = query_dna_rc
                    final_orf = longest_orf_rc
                    final_frame = frame_rc
                    detected_features = features_rc
                else:
                    alignment_text = alignment_text_fwd
                    identity_percent = identity_fwd
                    final_orf = longest_orf_fwd
                    final_frame = frame_fwd
                    detected_features = features_fwd

                _o_fwd = find_all_orfs(longest_orf_fwd or "")
                _o_rc  = find_all_orfs(longest_orf_rc or "")
                raw_reads = [
                    {
                        "label": f"{filename} (forward)",
                        "dna": query_dna, "length": len(query_dna),
                        "alignment_html": alignment_text_fwd,
                        "alignment_html_with_orfs": None,
                        "identity": identity_fwd,
                        "orf_count": len(_o_fwd),
                        "best_orf_length": max((len(o) for o, _ in _o_fwd), default=0),
                        "all_orfs": [{"sequence": o, "length": len(o), "start": s} for o, s in _o_fwd],
                    },
                    {
                        "label": f"{filename} (reverse complement)",
                        "dna": query_dna_rc, "length": len(query_dna_rc),
                        "alignment_html": alignment_text_rc,
                        "alignment_html_with_orfs": None,
                        "identity": identity_rc,
                        "orf_count": len(_o_rc),
                        "best_orf_length": max((len(o) for o, _ in _o_rc), default=0),
                        "all_orfs": [{"sequence": o, "length": len(o), "start": s} for o, s in _o_rc],
                    },
                ]
                investigate_html = build_investigate_html(
                    raw_reads, ref_seq, lcs=None, sample_name=filename
                )

            else:
                # No reference found - just translate
                final_orf, final_frame = get_longest_orf_six_frames(query_dna)
                detected_features = features_fwd
                _o_fwd = find_all_orfs(final_orf or "")
                raw_reads = [
                    {
                        "label": f"{filename} (forward)",
                        "dna": query_dna, "length": len(query_dna),
                        "alignment_html": None, "alignment_html_with_orfs": None,
                        "identity": 0.0,
                        "orf_count": len(_o_fwd),
                        "best_orf_length": max((len(o) for o, _ in _o_fwd), default=0),
                        "all_orfs": [{"sequence": o, "length": len(o), "start": s} for o, s in _o_fwd],
                    },
                ]

            _orf_dna = None
            if final_orf and final_dna:
                _ds, _de = _find_orf_dna_span(final_dna, final_orf)
                if _ds is not None and _de is not None:
                    _orf_dna = final_dna[_ds:_de]

            results.append({
                "sample": filename,
                "merged_dna": final_dna,
                "longest_orf": final_orf,
                "orf_dna": _orf_dna,
                "orf_length": len(final_orf) if final_orf else 0,
                "frame": final_frame,
                "merge_info": {},
                "reference": ref_seq,
                "matched_ref_name": matched_ref_name,
                "orientation": orientation,
                "available_refs": list(ref_dict.keys()) if not ref_seq else None,
                "tried_prefixes": candidates if not ref_seq else None,
                "alignment_text": alignment_text,
                "identity_percent": identity_percent,
                "detected_features": detected_features,
                "orf_count": final_orf_data.get("orf_count", 0),
                "alignment_text_with_orfs": final_orf_data.get("alignment_with_orfs"),
                "raw_reads": raw_reads,
                "investigate_html": investigate_html,
            })

        except Exception as e:
            results.append({
                "sample": filename,
                "error": str(e),
                "merged_dna": None,
                "longest_orf": None,
                "orf_dna": None,
                "orf_length": 0,
                "frame": None,
                "merge_info": {},
                "reference": None,
                "alignment_text": None,
                "identity_percent": 0.0,
                "detected_features": [],
                "orf_count": 0,
                "alignment_text_with_orfs": None,
                "raw_reads": [],
                "investigate_html": None,
            })

    return results


# ───────────────────────────────────────────────────────────────
#  Primer type auto-detection
# ───────────────────────────────────────────────────────────────

def detect_primer_type(filenames: List[str]) -> str:
    """
    Auto-detect primer type from a list of filenames.

    Returns "protein" if T7/T7-Term patterns found,
    "intron" if pUC19/M13R patterns found, "unknown" otherwise.
    """
    has_t7 = any(re.search(r't7', f, re.IGNORECASE) for f in filenames)
    has_puc19 = any(re.search(r'puc19', f, re.IGNORECASE) for f in filenames)
    has_m13r = any(re.search(r'm13r', f, re.IGNORECASE) for f in filenames)

    if has_t7:
        return "protein"
    if has_puc19 or has_m13r:
        return "intron"
    return "unknown"


def detect_single_primer_type(filename: str) -> str:
    """
    Auto-detect primer type from a single filename.

    Returns "protein" if T7/T7-Term patterns found,
    "intron" if pUC19/M13R patterns found, "unknown" otherwise.
    """
    if re.search(r't7', filename, re.IGNORECASE):
        return "protein"
    if re.search(r'(puc19|m13r)', filename, re.IGNORECASE):
        return "intron"
    return "unknown"


# ───────────────────────────────────────────────────────────────
#  Unified processing pipeline
# ───────────────────────────────────────────────────────────────

def _add_orf_info(result: Dict[str, Any], dna: Optional[str]) -> Dict[str, Any]:
    """Run ORF detection on DNA and add results to the dict."""
    if not dna:
        result.setdefault("longest_orf", None)
        result.setdefault("orf_dna", None)
        result.setdefault("orf_length", 0)
        result.setdefault("frame", None)
        result["has_orf"] = False
        return result

    # If result already has an ORF from protein pipeline, check significance
    existing_orf = result.get("longest_orf")
    if existing_orf:
        result["orf_length"] = len(existing_orf)
        result["has_orf"] = len(existing_orf) >= ORF_MIN_AA
        if "orf_dna" not in result:
            ds, de = _find_orf_dna_span(dna, existing_orf)
            result["orf_dna"] = dna[ds:de] if ds is not None and de is not None else None
        return result

    # Run ORF detection
    orf, frame = get_longest_orf_six_frames(dna)
    result["longest_orf"] = orf if orf else None
    result["orf_length"] = len(orf) if orf else 0
    result["frame"] = frame
    result["has_orf"] = len(orf) >= ORF_MIN_AA if orf else False
    if orf:
        ds, de = _find_orf_dna_span(dna, orf)
        result["orf_dna"] = dna[ds:de] if ds is not None and de is not None else None
    else:
        result["orf_dna"] = None
    return result


def unified_main(seq_files: List[Any], ref_file: Union[str, Any],
                 min_lcs: int = 12, primer_type: str = "auto",
                 use_aa_refs: bool = False,
                 selected_categories: Optional[List[str]] = None) -> List[Dict[str, Any]]:
    """
    Unified processing for paired sequences.

    Auto-detects primer type (T7/T7-Term vs pUC19/M13R) and routes to
    the appropriate merge/alignment pipeline. Always runs ORF detection.
    """
    from intron_modules import (
        find_puc19_m13r_pairs,
        load_reference_xlsx as intron_load_ref,
        align_dna_local,
    )

    if not seq_files:
        raise ValueError("No sequence files provided")

    filenames = [f.name for f in seq_files]

    # Auto-detect primer type
    if primer_type == "auto":
        primer_type = detect_primer_type(filenames)
        if primer_type == "unknown":
            primer_type = "protein"  # default fallback

    if primer_type == "intron":
        # ─── Intron pipeline (pUC19/M13R pairs) ───
        pairs = find_puc19_m13r_pairs(filenames)
        if not pairs:
            raise ValueError(
                "No valid pUC19/M13R pairs found. "
                "Files should contain 'pUC19' and 'M13R' in their names."
            )

        ref_dict = intron_load_ref(ref_file)
        all_results = []

        for fwd_name, rev_name, prefix in pairs:
            fwd_file = next(f for f in seq_files if f.name == fwd_name)
            rev_file = next(f for f in seq_files if f.name == rev_name)

            try:
                fwd_dna = read_seq_file(fwd_file)
                rev_dna_raw = read_seq_file(rev_file)
                rev_dna = str(Seq(rev_dna_raw).reverse_complement())

                merge_info = merge_by_lcs_dna(fwd_dna, rev_dna, min_lcs=min_lcs)
                merged_dna = merge_info["merged"]

                if not merged_dna:
                    result = {
                        "sample": prefix, "merged_dna": None, "final_dna": None,
                        "merge_info": merge_info, "crop_info": {},
                        "reference": None, "alignment_text": None,
                        "identity_percent": 0.0, "orientation": "forward",
                        "pipeline": "intron", "detected_features": [],
                    }
                    _add_orf_info(result, None)
                    all_results.append(result)
                    continue

                # Match reference
                sample_base = re.split(r'[_-]', prefix, 1)[0].lower()
                candidates = [
                    sample_base, sample_base.replace("-", ""),
                    sample_base.split("-")[0], prefix.lower(),
                ]
                ref_seq, matched_ref_name = None, None
                for c in candidates:
                    if c in ref_dict:
                        ref_seq = ref_dict[c]
                        matched_ref_name = c
                        break
                if not ref_seq:
                    for rn, rs in ref_dict.items():
                        if any(c in rn or rn in c for c in candidates if c):
                            ref_seq = rs
                            matched_ref_name = rn
                            break

                # Detect features & align both orientations
                features_fwd = detect_features(merged_dna, selected_categories=selected_categories)
                merged_rc = str(Seq(merged_dna).reverse_complement())
                features_rc = detect_features(merged_rc, selected_categories=selected_categories)

                alignment_text = None
                identity_percent = 0.0
                orientation = "forward"
                final_dna = merged_dna
                detected_features = features_fwd

                if ref_seq:
                    at_fwd, id_fwd = align_dna_local(merged_dna, ref_seq, prefix, features=features_fwd)
                    at_rc, id_rc = align_dna_local(merged_rc, ref_seq, f"{prefix} (RC)", features=features_rc)
                    if id_rc > id_fwd:
                        alignment_text, identity_percent = at_rc, id_rc
                        orientation = "reverse complement"
                        final_dna = merged_rc
                        detected_features = features_rc
                    else:
                        alignment_text, identity_percent = at_fwd, id_fwd

                result = {
                    "sample": prefix, "merged_dna": merged_dna,
                    "final_dna": final_dna, "merge_info": merge_info,
                    "crop_info": {}, "reference": ref_seq,
                    "matched_ref_name": matched_ref_name,
                    "alignment_text": alignment_text,
                    "identity_percent": identity_percent,
                    "orientation": orientation, "pipeline": "intron",
                    "detected_features": detected_features,
                }
                _add_orf_info(result, final_dna)
                all_results.append(result)

            except Exception as e:
                result = {
                    "sample": prefix, "error": str(e),
                    "merged_dna": None, "final_dna": None,
                    "merge_info": {}, "crop_info": {},
                    "reference": None, "alignment_text": None,
                    "identity_percent": 0.0, "orientation": "forward",
                    "pipeline": "intron", "detected_features": [],
                }
                _add_orf_info(result, None)
                all_results.append(result)

        return all_results

    else:
        # ─── Protein pipeline (T7/T7-Term pairs) ───
        pairs = find_pairs(filenames)
        if not pairs:
            raise ValueError(
                "No valid T7/T7-Term pairs found. "
                "Files should follow naming pattern: *-T7.seq and *-T7-Term.seq"
            )

        # Load references
        if use_aa_refs:
            ref_dict = load_reference_xlsx(ref_file)
            ref_dna_dict = None
        else:
            ref_dna_dict = load_reference_dna_xlsx(ref_file)
            ref_dict = None

        all_results = []

        for fwd_name, rev_name, prefix in pairs:
            fwd_file = next(f for f in seq_files if f.name == fwd_name)
            rev_file = next(f for f in seq_files if f.name == rev_name)

            try:
                result = process_pair(fwd_file, rev_file, prefix, min_lcs=min_lcs,
                                      ref_dict=ref_dict, ref_dna_dict=ref_dna_dict)
                result["sample"] = prefix
                result["pipeline"] = "protein"
                result["orientation"] = "forward"
                result["detected_features"] = detect_features(
                    result.get("merged_dna") or "", selected_categories=selected_categories
                ) if result.get("merged_dna") else []
                result["crop_info"] = {}
                result["final_dna"] = result.get("merged_dna")

                sample_base = re.split(r'[_-]', prefix, 1)[0].lower()
                if ref_dna_dict:
                    result["reference"] = ref_dna_dict.get(sample_base)
                elif ref_dict:
                    result["reference"] = ref_dict.get(sample_base)

                _add_orf_info(result, result.get("merged_dna"))
                all_results.append(result)

            except Exception as e:
                result = {
                    "sample": prefix, "error": str(e),
                    "merged_dna": None, "final_dna": None,
                    "longest_orf": None, "orf_dna": None, "orf_length": 0,
                    "frame": None, "has_orf": False,
                    "merge_info": {}, "crop_info": {},
                    "reference": None, "alignment_text": None,
                    "identity_percent": 0.0, "orientation": "forward",
                    "pipeline": "protein", "detected_features": [],
                }
                all_results.append(result)

        return all_results


def unified_circular(seq_files: List[Any], ref_file: Union[str, Any],
                     use_aa_refs: bool = False,
                     selected_categories: Optional[List[str]] = None) -> List[Dict[str, Any]]:
    """
    Unified processing for circular plasmid sequences.

    Processes using DNA alignment with both orientations, always runs ORF detection.
    """
    results = process_circular_queries(seq_files, ref_file,
                                       use_dna_alignment=not use_aa_refs,
                                       selected_categories=selected_categories)

    for r in results:
        r["pipeline"] = "protein"
        r["crop_info"] = {}
        r["final_dna"] = r.get("merged_dna")
        _add_orf_info(r, r.get("merged_dna"))

    return results


def unified_unidirectional(seq_files: List[Any], ref_file: Union[str, Any],
                           primer_type: str = "auto",
                           use_aa_refs: bool = False,
                           selected_categories: Optional[List[str]] = None) -> List[Dict[str, Any]]:
    """
    Unified processing for unidirectional sequences.

    Auto-detects primer type per file and processes accordingly.
    Always runs ORF detection.
    """
    from intron_modules import (
        load_reference_xlsx as intron_load_ref,
        align_dna_local,
    )

    if not seq_files:
        raise ValueError("No sequence files provided")

    # Determine which pipeline per file
    protein_files = []
    intron_files = []

    for f in seq_files:
        if primer_type == "auto":
            pt = detect_single_primer_type(f.name)
        else:
            pt = primer_type

        if pt == "intron":
            intron_files.append(f)
        else:
            protein_files.append(f)

    all_results = []

    # Process protein-type unidirectional files
    if protein_files:
        prot_results = process_unidirectional(
            protein_files, ref_file,
            use_dna_alignment=not use_aa_refs,
            selected_categories=selected_categories
        )
        for r in prot_results:
            r["pipeline"] = "protein"
            r["crop_info"] = {}
            r["final_dna"] = r.get("merged_dna")
            _add_orf_info(r, r.get("merged_dna"))
        all_results.extend(prot_results)

    # Process intron-type unidirectional files
    if intron_files:
        ref_dict = intron_load_ref(ref_file)

        for seq_file in intron_files:
            filename = seq_file.name

            try:
                query_dna = read_seq_file(seq_file)

                # Check if reverse primer
                is_reverse = bool(re.search(r'm13r', filename, re.IGNORECASE))
                if is_reverse:
                    query_dna = str(Seq(query_dna).reverse_complement())

                # Match reference
                base_name = filename.rsplit(".", 1)[0].lower()
                base_name = re.sub(r'\s*\(\d+\)$', '', base_name)  # strip Windows dup suffix e.g. " (1)"
                clean_name = re.sub(r'[-_]?(puc19(seqback)?|m13r|m13f)$', '', base_name, re.IGNORECASE)
                prefix = re.split(r'[_-]', clean_name, 1)[0].lower()
                candidates = [clean_name, prefix, prefix.replace("-", ""), prefix.split("-")[0]]

                ref_seq, matched_ref_name = None, None
                for c in candidates:
                    if c in ref_dict:
                        ref_seq = ref_dict[c]
                        matched_ref_name = c
                        break
                if not ref_seq:
                    for rn, rs in ref_dict.items():
                        if any(c in rn or rn in c for c in candidates if c):
                            ref_seq = rs
                            matched_ref_name = rn
                            break

                # Detect features both orientations
                features_fwd = detect_features(query_dna, selected_categories=selected_categories)
                query_rc = str(Seq(query_dna).reverse_complement())
                features_rc = detect_features(query_rc, selected_categories=selected_categories)

                alignment_text = None
                identity_percent = 0.0
                orientation = "M13R (auto-flipped)" if is_reverse else "forward"
                final_dna = query_dna
                detected_features = features_fwd

                if ref_seq:
                    at_fwd, id_fwd = align_dna_local(query_dna, ref_seq, filename, features=features_fwd)
                    at_rc, id_rc = align_dna_local(query_rc, ref_seq, f"{filename} (RC)", features=features_rc)
                    if id_rc > id_fwd:
                        alignment_text, identity_percent = at_rc, id_rc
                        orientation = "forward (flipped back)" if is_reverse else "reverse complement"
                        final_dna = query_rc
                        detected_features = features_rc
                    else:
                        alignment_text, identity_percent = at_fwd, id_fwd

                result = {
                    "sample": filename, "merged_dna": final_dna,
                    "final_dna": final_dna, "merge_info": {},
                    "crop_info": {}, "reference": ref_seq,
                    "matched_ref_name": matched_ref_name,
                    "alignment_text": alignment_text,
                    "identity_percent": identity_percent,
                    "orientation": orientation, "pipeline": "intron",
                    "detected_features": detected_features,
                }
                _add_orf_info(result, final_dna)
                all_results.append(result)

            except Exception as e:
                result = {
                    "sample": filename, "error": str(e),
                    "merged_dna": None, "final_dna": None,
                    "merge_info": {}, "crop_info": {},
                    "reference": None, "alignment_text": None,
                    "identity_percent": 0.0, "orientation": "forward",
                    "pipeline": "intron", "detected_features": [],
                }
                _add_orf_info(result, None)
                all_results.append(result)

    return all_results
