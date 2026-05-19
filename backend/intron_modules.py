#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 24 13:28:27 2025

@author: JessicaGeneWeaver

Intron/Plasmid DNA sequence processing pipeline:
- Merge forward (pUC19) and reverse (M13R) DNA sequences using LCS strategy
- Crop merged sequences using start and end motifs
- Compare cropped sequences with reference sequences
- Generate alignment visualizations and identity metrics

Reference Excel file requirements:
    refseq.xlsx with columns:
        "ID"  → sample prefix (before underscore, case-insensitive)
        "Sequence"       → reference DNA sequence
"""

import re
from typing import Union, Tuple, Dict, List, Any, Optional
from io import StringIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
import pandas as pd

from feature_detection import (
    detect_features,
    map_features_to_gapped_alignment,
    build_feature_position_map,
    DetectedFeature,
    FeatureDefinition,
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

    Finds the LCS between forward and reverse sequences and merges them by
    trimming the forward sequence at the LCS end and the reverse sequence
    at the LCS start.

    Args:
        fwd_dna: Forward DNA sequence
        rev_dna: Reverse DNA sequence (already reverse complemented)
        min_lcs: Minimum LCS length required for merging (default: 12 bp)

    Returns:
        Dictionary with keys:
            - merged: Merged DNA sequence or None if LCS too short
            - lcs: The longest common substring found
            - lcs_length: Length of the LCS
            - fwd_trimmed_len: Length of forward sequence used
            - rev_trimmed_len: Length of reverse sequence used

    Raises:
        ValueError: If input sequences are empty
    """
    if not fwd_dna or not rev_dna:
        raise ValueError("Forward and reverse sequences cannot be empty")

    lcs, fwd_start, fwd_end = longest_common_substring(fwd_dna, rev_dna)
    if not lcs or len(lcs) < min_lcs:
        return {
            "merged": None,
            "lcs": "",
            "lcs_length": 0,
            "fwd_trimmed_len": 0,
            "rev_trimmed_len": 0
        }

    rev_pos = rev_dna.find(lcs)
    trimmed_fwd = fwd_dna[:fwd_end]
    trimmed_rev = rev_dna[rev_pos + len(lcs):]
    merged = trimmed_fwd + trimmed_rev

    return {
        "merged": merged,
        "lcs": lcs,
        "lcs_length": len(lcs),
        "fwd_trimmed_len": len(trimmed_fwd),
        "rev_trimmed_len": len(trimmed_rev)
    }


def crop_sequence(seq: str, start_motif: str, end_motif: str) -> Tuple[Optional[str], Optional[int], Optional[int]]:
    """
    Crop a DNA sequence between start and end motifs.

    Searches for the first occurrence of start_motif and the last occurrence
    of end_motif, then extracts the sequence between them (inclusive).

    Args:
        seq: DNA sequence to crop
        start_motif: Starting motif sequence
        end_motif: Ending motif sequence

    Returns:
        Tuple of (cropped_sequence, start_position, end_position).
        Returns (None, None, None) if either motif is not found or if
        end position is before start position.
        Positions are 1-indexed.
    """
    if not seq:
        raise ValueError("Sequence cannot be empty")
    if not start_motif or not end_motif:
        raise ValueError("Start and end motifs cannot be empty")

    start_index = seq.find(start_motif)
    end_index = seq.rfind(end_motif)
    if start_index == -1 or end_index == -1:
        return None, None, None
    end_index += len(end_motif)
    if end_index <= start_index:
        return None, None, None
    cropped = seq[start_index:end_index]
    return cropped, start_index + 1, end_index

def find_puc19_m13r_pairs(files: List[str]) -> List[Tuple[str, str, str]]:
    """
    Identify forward/reverse sequence pairs based on filename patterns.

    Forward files contain 'pUC19' (case-insensitive) in the filename.
    Reverse files contain 'M13R' (case-insensitive) in the filename.
    Files are matched based on their common prefix.

    Args:
        files: List of filenames to search

    Returns:
        List of tuples: (forward_file, reverse_file, prefix)
    """
    pairs = []
    # Case-insensitive search for pUC19 and M13R patterns
    fwd_files = [f for f in files if re.search(r'puc19', f, re.IGNORECASE)]
    rev_files = [f for f in files if re.search(r'm13r', f, re.IGNORECASE)]

    for fwd in fwd_files:
        # Extract prefix before pUC19 (case-insensitive)
        match_obj = re.search(r'(.+?)puc19', fwd, re.IGNORECASE)
        if not match_obj:
            continue
        prefix = match_obj.group(1)

        # Try to find matching reverse file with same prefix
        # Check both exact prefix match and with common separators
        for rev in rev_files:
            if rev.lower().startswith(prefix.lower()):
                pairs.append((fwd, rev, prefix.rstrip('_-')))
                break

    return pairs

# ───────────────────────────────────────────────────────────────
#  Reference handling and alignment
# ───────────────────────────────────────────────────────────────

def load_reference_xlsx(filepath: Union[str, Any]) -> Dict[str, str]:
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
    if "clone number" not in df.columns or "dna sequence" not in df.columns:
        raise ValueError("Project file must have headers 'Clone number' and 'DNA sequence' (case-insensitive).")

    def _clone_key(val):
        s = str(val).strip()
        if s.endswith(".0"):
            s = s[:-2]
        return s.lower()

    ref_dict = {
        _clone_key(row["clone number"]): str(row["dna sequence"]).strip()
        for _, row in df.iterrows()
        if pd.notna(row["clone number"]) and pd.notna(row["dna sequence"])
    }

    if not ref_dict:
        raise ValueError("No valid reference sequences found in Excel file")

    return ref_dict


def visualize_merge_vs_reference(merge_seq: str, ref_seq: str, sample_name: str,
                                block_size: int = 50,
                                query_circular: bool = False,
                                features: Optional[List[DetectedFeature]] = None) -> Tuple[str, float]:
    """
    Create an HTML-formatted alignment between query sequence and reference sequence.

    Performs global pairwise alignment and generates a color-coded alignment
    visualization with position numbers. Matches are highlighted in green,
    mismatches in red, and gaps in gray. Detected features are highlighted
    with colored backgrounds.

    For circular query sequences, the sequence is doubled before alignment to allow
    alignments that span the arbitrary start/end boundary.

    Args:
        merge_seq: Query DNA sequence
        ref_seq: Reference DNA sequence
        sample_name: Name of the sample for display
        block_size: Number of nucleotides per line in output (default: 50)
        query_circular: If True, treat query as circular (doubles query for alignment)
        features: Optional list of detected features to highlight

    Returns:
        Tuple of (alignment_html, identity_percent) where alignment_html
        is an HTML formatted string and identity_percent is 0-100
    """
    original_query_len = len(merge_seq)

    # For circular queries, double the sequence to allow wraparound alignments
    query_for_align = merge_seq + merge_seq if query_circular else merge_seq

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

    # Map features to gapped alignment positions if provided
    feature_map = {}
    if features:
        mapped_features = map_features_to_gapped_alignment(features, seqA_str)
        feature_map = build_feature_position_map(mapped_features, use_gapped=True)

    # Step 2: Compute identity
    matches = sum(a == b and a != '-' for a, b in zip(seqA_str, seqB_str))
    aligned_len = max(len(seqA_str), len(seqB_str))
    identity = (matches / aligned_len) * 100 if aligned_len else 0.0

    # Step 3: Build HTML output
    circular_note = " (circular query, local alignment)" if query_circular else ""

    html_parts = [
        """<style>
        .dna-match { color: #2a9d8f; font-weight: bold; }
        .dna-mismatch { color: #e76f51; font-weight: bold; }
        .dna-gap { color: #f4a261; }
        .dna-match-bg { color: #264653; background-color: #b7e4dd; font-weight: bold; }
        .dna-mismatch-bg { color: #264653; background-color: #f4c4b5; font-weight: bold; }
        .align-match-pipe { color: #2a9d8f; }
        </style>""",
        "<div style='font-family: monospace; font-size: 14px; line-height: 1.6; white-space: pre;'>",
        f"<div style='font-weight: bold; margin-bottom: 10px; white-space: normal; color: #264653;'>Alignment for sample: {sample_name}{circular_note}</div>",
        f"<div style='margin-bottom: 15px; white-space: normal; color: #264653;'>Alignment length: {aligned_len} | Identity: {identity:.1f}%</div>",
    ]

    # Show original length for circular query
    if query_circular:
        start_pos = (query_start_offset % original_query_len) + 1 if original_query_len > 0 else query_start_offset + 1
        html_parts.append(f"<div style='margin-bottom: 10px; white-space: normal; color: #264653;'>Original query length: {original_query_len} bp | Alignment starts at position {start_pos}</div>")

    q_counter = query_start_offset + 1  # Start from alignment position in query
    r_counter = 1
    label_width = 20

    for i in range(0, len(seqA_str), block_size):
        blockA = seqA_str[i:i+block_size]
        blockB = seqB_str[i:i+block_size]

        # Count non-gap characters to get end positions
        q_chars_in_block = len([c for c in blockA if c != '-'])
        r_chars_in_block = len([c for c in blockB if c != '-'])

        # For circular sequences, wrap position numbers
        if query_circular and original_query_len > 0:
            q_start_display = ((q_counter - 1) % original_query_len) + 1
            q_end_display = ((q_counter + q_chars_in_block - 2) % original_query_len) + 1
        else:
            q_start_display = q_counter
            q_end_display = q_counter + q_chars_in_block - 1

        r_start_display = r_counter
        r_end_display = r_counter + r_chars_in_block - 1

        q_label = f"Query {q_start_display}-{q_end_display}:".ljust(label_width)
        r_label = f"Ref   {r_start_display}-{r_end_display}:".ljust(label_width)

        # Build color-coded sequences
        query_html = ""
        ref_html = ""
        match_html = ""

        for j, (a, b) in enumerate(zip(blockA, blockB)):
            gapped_pos = i + j
            feature = feature_map.get(gapped_pos)

            # Build background style for feature highlighting
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
                query_html += f"<span class='dna-match-bg' style='{bg_style}'{title_attr}>{a}</span>"
            else:
                query_html += f"<span class='dna-mismatch-bg' style='{bg_style}'{title_attr}>{a}</span>"

            # Reference sequence coloring (no feature highlighting)
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

        # Update counters for next block
        q_counter += q_chars_in_block
        r_counter += r_chars_in_block

    html_parts.append("</div>")
    return "".join(html_parts), identity


def align_dna_local(query_dna: str, ref_dna: str, sample_name: str,
                    block_size: int = 50,
                    features: Optional[List[DetectedFeature]] = None) -> Tuple[str, float]:
    """
    Perform local DNA-to-DNA alignment for partial/unidirectional sequences.

    Uses local alignment to find the best matching region between query and reference,
    then trims to only the region where the reference has sequence. This is ideal for
    unidirectional sequences that won't cover the full reference.

    Args:
        query_dna: Query DNA sequence
        ref_dna: Reference DNA sequence
        sample_name: Name of the sample for display
        block_size: Number of nucleotides per line in output (default: 50)
        features: Optional list of detected features to highlight

    Returns:
        Tuple of (alignment_html, identity_percent) where alignment_html
        is an HTML formatted string and identity_percent is 0-100
    """
    original_query_len = len(query_dna)

    # Use local alignment to find best matching region
    alignments = pairwise2.align.localms(query_dna, ref_dna, 2, -1, -2, -0.5)

    if not alignments:
        return f"<p>No alignment found for sample: {sample_name}</p>", 0.0

    alignment = alignments[0]

    seqA_full = alignment[0]  # Query DNA with gaps
    seqB_full = alignment[1]  # Reference DNA with gaps

    # Extract the aligned region from local alignment
    begin_pos = int(alignment[3])
    end_pos = int(alignment[4])
    seqA_aligned = seqA_full[begin_pos:end_pos]
    seqB_aligned = seqB_full[begin_pos:end_pos]

    # Trim to only the region where reference has sequence (no leading/trailing ref gaps)
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

    # Map features to gapped alignment positions if provided
    feature_map = {}
    if features:
        # Map features to the gapped sequence
        mapped_features = map_features_to_gapped_alignment(features, seqA_str)
        feature_map = build_feature_position_map(mapped_features, use_gapped=True)

    # Compute DNA identity only on reference-covered region
    matches = sum(a == b and a != '-' for a, b in zip(seqA_str, seqB_str))
    ref_len_no_gaps = len(seqB_str.replace('-', ''))
    identity = (matches / ref_len_no_gaps) * 100 if ref_len_no_gaps else 0.0

    # Build HTML output
    aligned_region_len = len(seqA_str)

    html_parts = [
        """<style>
        .dna-match { color: #2a9d8f; font-weight: bold; }
        .dna-mismatch { color: #e76f51; font-weight: bold; }
        .dna-gap { color: #f4a261; }
        .dna-match-bg { color: #264653; background-color: #b7e4dd; font-weight: bold; }
        .dna-mismatch-bg { color: #264653; background-color: #f4c4b5; font-weight: bold; }
        .align-match-pipe { color: #2a9d8f; }
        </style>""",
        "<div style='font-family: monospace; font-size: 14px; line-height: 1.6; white-space: pre;'>",
        f"<div style='font-weight: bold; margin-bottom: 10px; white-space: normal; color: #264653;'>DNA Alignment for sample: {sample_name} (local alignment)</div>",
        f"<div style='margin-bottom: 10px; white-space: normal; color: #264653;'>Reference-aligned region: {aligned_region_len} bp | DNA Identity: {identity:.1f}% (vs {ref_len_no_gaps} bp reference)</div>",
        f"<div style='margin-bottom: 10px; white-space: normal; color: #264653;'>Original query length: {original_query_len} bp | Alignment starts at query position {query_start_offset + 1}</div>",
    ]

    q_counter = query_start_offset + 1
    r_counter = 1
    label_width = 22

    for i in range(0, len(seqA_str), block_size):
        blockA = seqA_str[i:i+block_size]
        blockB = seqB_str[i:i+block_size]

        q_end = q_counter + len([c for c in blockA if c != '-']) - 1
        r_end = r_counter + len([c for c in blockB if c != '-']) - 1

        q_label = f"Query  {q_counter}-{q_end}:".ljust(label_width)
        r_label = f"Ref    {r_counter}-{r_end}:".ljust(label_width)

        query_html = ""
        ref_html = ""
        match_html = ""

        for j, (a, b) in enumerate(zip(blockA, blockB)):
            gapped_pos = i + j
            feature = feature_map.get(gapped_pos)

            # Build background style for feature highlighting
            if feature:
                bg_style = f"background-color: {feature.color}40;"  # 40 = 25% opacity
                title_attr = f" title='{feature.name}'"
            else:
                bg_style = ""
                title_attr = ""

            if a == '-':
                query_html += f"<span class='dna-gap' style='{bg_style}'{title_attr}>{a}</span>"
            elif a == b:
                query_html += f"<span class='dna-match-bg' style='{bg_style}'{title_attr}>{a}</span>"
            else:
                query_html += f"<span class='dna-mismatch-bg' style='{bg_style}'{title_attr}>{a}</span>"

            if b == '-':
                ref_html += f"<span class='dna-gap'>{b}</span>"
            elif a == b:
                ref_html += f"<span class='dna-match'>{b}</span>"
            else:
                ref_html += f"<span class='dna-mismatch'>{b}</span>"

            if a == b and a != '-':
                match_html += "<span class='align-match-pipe'>|</span>"
            else:
                match_html += " "

        html_parts.append("<div style='margin-bottom: 12px;'>")
        html_parts.append(f"<div>{q_label}{query_html}</div>")
        html_parts.append(f"<div>{' ' * label_width}{match_html}</div>")
        html_parts.append(f"<div>{r_label}{ref_html}</div>")
        html_parts.append("</div>")

        q_counter += len([c for c in blockA if c != '-'])
        r_counter += len([c for c in blockB if c != '-'])

    html_parts.append("</div>")
    return "".join(html_parts), identity


# ───────────────────────────────────────────────────────────────
#  Main processing
# ───────────────────────────────────────────────────────────────

def process_circular_queries(seq_files: List[Any], ref_file: Union[str, Any],
                             selected_categories: Optional[List[str]] = None) -> List[Dict[str, Any]]:
    """
    Process circular query sequences (whole plasmids in FASTA format).

    Skips pair-finding, merging, and cropping - aligns directly to references.
    Matches query files to references by filename prefix (text before first underscore).

    Args:
        seq_files: List of uploaded FASTA files (file-like objects with .name attribute)
        ref_file: Excel file with reference sequences (file-like object or path)
        selected_categories: List of feature categories to detect (None = all)

    Returns:
        List of result dictionaries, one per processed file, each containing:
            - sample: Filename
            - merged_dna: None (no merging for circular)
            - cropped_dna: The query sequence
            - merge_info: Empty dict
            - crop_info: Empty dict
            - reference: Reference sequence or None
            - alignment_text: Alignment visualization or None
            - identity_percent: Percent identity to reference
            - detected_features: List of detected DNA features
    """
    if not seq_files:
        raise ValueError("No sequence files provided")

    # Load reference sequences
    ref_dict = load_reference_xlsx(ref_file)

    results = []
    for seq_file in seq_files:
        filename = seq_file.name

        try:
            # Read the FASTA sequence
            query_seq = read_seq_file(seq_file)

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

            # If still no match, try partial matching (reference name contains candidate or vice versa)
            if not ref_seq:
                for ref_name, ref_sequence in ref_dict.items():
                    if any(c in ref_name or ref_name in c for c in candidates if c):
                        ref_seq = ref_sequence
                        matched_ref_name = ref_name
                        break

            # Detect features on both orientations
            features_fwd = detect_features(query_seq, selected_categories=selected_categories)
            query_seq_rc = str(Seq(query_seq).reverse_complement())
            features_rc = detect_features(query_seq_rc, selected_categories=selected_categories)

            # Generate alignment if reference found
            # Always try both orientations and use the better alignment
            alignment_text = None
            identity_percent = 0.0
            orientation = "forward"
            final_query_seq = query_seq
            detected_features = features_fwd

            if ref_seq:
                # Try forward orientation
                alignment_text_fwd, identity_fwd = visualize_merge_vs_reference(
                    query_seq, ref_seq, filename,
                    query_circular=True, features=features_fwd
                )

                # Try reverse complement
                alignment_text_rc, identity_rc = visualize_merge_vs_reference(
                    query_seq_rc, ref_seq, f"{filename} (reverse complement)",
                    query_circular=True, features=features_rc
                )

                # Use the orientation with better identity
                if identity_rc > identity_fwd:
                    alignment_text = alignment_text_rc
                    identity_percent = identity_rc
                    orientation = "reverse complement"
                    final_query_seq = query_seq_rc
                    detected_features = features_rc
                else:
                    alignment_text = alignment_text_fwd
                    identity_percent = identity_fwd

            results.append({
                "sample": filename,
                "merged_dna": None,
                "cropped_dna": final_query_seq,
                "merge_info": {},
                "crop_info": {},
                "reference": ref_seq,
                "matched_ref_name": matched_ref_name,
                "orientation": orientation,
                "available_refs": list(ref_dict.keys()) if not ref_seq else None,
                "tried_prefixes": candidates if not ref_seq else None,
                "alignment_text": alignment_text,
                "identity_percent": identity_percent,
                "detected_features": detected_features
            })

        except Exception as e:
            results.append({
                "sample": filename,
                "error": str(e),
                "merged_dna": None,
                "cropped_dna": None,
                "merge_info": {},
                "crop_info": {},
                "reference": None,
                "alignment_text": None,
                "identity_percent": 0.0,
                "detected_features": []
            })

    return results


def main(seq_files: List[Any], ref_file: Union[str, Any], min_lcs: int = 12,
         selected_categories: Optional[List[str]] = None) -> List[Dict[str, Any]]:
    """
    Main processing function for intron/plasmid sequence analysis.

    Identifies pUC19/M13R pairs from uploaded files, merges them,
    and compares to reference sequences using local alignment.

    Args:
        seq_files: List of uploaded sequence files (file-like objects with .name attribute)
        ref_file: Excel file with reference sequences (file-like object or path)
        min_lcs: Minimum LCS length for merging (default: 12 bp)
        selected_categories: List of feature categories to detect (None = all)

    Returns:
        List of result dictionaries, one per processed pair, each containing:
            - sample: Sample name/prefix
            - merged_dna: Merged DNA sequence or None
            - cropped_dna: Same as merged_dna (no cropping)
            - merge_info: Details about the merge (LCS info)
            - crop_info: Empty dict (no cropping)
            - reference: Reference sequence or None
            - alignment_text: Alignment visualization or None
            - identity_percent: Percent identity to reference
            - detected_features: List of detected DNA features
    """
    # Validate inputs
    if not seq_files:
        raise ValueError("No sequence files provided")

    # Load reference sequences
    ref_dict = load_reference_xlsx(ref_file)

    # Find pairs
    files = [f.name for f in seq_files]  # names for pairing
    pairs = find_puc19_m13r_pairs(files)

    if not pairs:
        # Provide detailed error message with file names
        puc19_files = [f for f in files if "pUC19" in f or "puc19" in f.lower()]
        m13r_files = [f for f in files if "M13R" in f or "m13r" in f.lower()]
        error_msg = "No valid pUC19/M13R pairs found.\n"
        error_msg += f"Files uploaded: {files}\n"
        error_msg += f"Files with pUC19: {puc19_files}\n"
        error_msg += f"Files with M13R: {m13r_files}\n"
        error_msg += "Files should contain 'pUC19' and 'M13R' in their names."
        raise ValueError(error_msg)

    # Process each pair
    all_results = []
    for fwd_name, rev_name, prefix in pairs:
        # Map names back to uploaded files
        fwd_file = next(f for f in seq_files if f.name == fwd_name)
        rev_file = next(f for f in seq_files if f.name == rev_name)

        try:
            # Read sequences
            fwd_dna = read_seq_file(fwd_file)
            rev_dna_raw = read_seq_file(rev_file)

            # Merge sequences
            rev_dna = str(Seq(rev_dna_raw).reverse_complement())
            merge_info = merge_by_lcs_dna(fwd_dna, rev_dna, min_lcs=min_lcs)
            merged_dna = merge_info["merged"]

            if not merged_dna:
                result = {
                    "sample": prefix,
                    "merged_dna": None,
                    "cropped_dna": None,
                    "merge_info": merge_info,
                    "crop_info": {"start": None, "end": None},
                    "alignment_text": None,
                    "identity_percent": 0.0
                }
                all_results.append(result)
                continue

            # Match by sample prefix - try multiple strategies
            sample_base = re.split(r'[_-]', prefix, 1)[0].lower()
            candidates = [
                sample_base,                                # e.g., "10868469-16"
                sample_base.replace("-", ""),               # e.g., "1086846916" (no hyphen)
                sample_base.split("-")[0],                  # e.g., "10868469" (before hyphen)
                prefix.lower(),                             # Full prefix
            ]

            ref_seq = None
            matched_ref_name = None
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

            # Detect features on merged DNA sequence
            features_fwd = detect_features(merged_dna, selected_categories=selected_categories)
            merged_rc = str(Seq(merged_dna).reverse_complement())
            features_rc = detect_features(merged_rc, selected_categories=selected_categories)

            # Compare with reference using local alignment - always try both orientations
            alignment_text = None
            identity_percent = 0.0
            orientation = "forward"
            final_merged = merged_dna
            detected_features = features_fwd

            if ref_seq:
                # Try forward orientation (local alignment)
                alignment_text_fwd, identity_fwd = align_dna_local(
                    merged_dna, ref_seq, prefix, features=features_fwd
                )

                # Try reverse complement
                alignment_text_rc, identity_rc = align_dna_local(
                    merged_rc, ref_seq, f"{prefix} (reverse complement)", features=features_rc
                )

                # Use the orientation with better identity
                if identity_rc > identity_fwd:
                    alignment_text = alignment_text_rc
                    identity_percent = identity_rc
                    orientation = "reverse complement"
                    final_merged = merged_rc
                    detected_features = features_rc
                else:
                    alignment_text = alignment_text_fwd
                    identity_percent = identity_fwd

            result = {
                "sample": prefix,
                "merged_dna": final_merged,
                "cropped_dna": final_merged,  # No cropping - use merged sequence
                "merge_info": merge_info,
                "crop_info": {},  # No cropping info
                "reference": ref_seq,
                "matched_ref_name": matched_ref_name,
                "orientation": orientation,
                "available_refs": list(ref_dict.keys()) if not ref_seq else None,
                "tried_prefixes": candidates if not ref_seq else None,
                "alignment_text": alignment_text,
                "identity_percent": identity_percent,
                "detected_features": detected_features
            }
            all_results.append(result)

        except Exception as e:
            # Log error but continue processing other pairs
            result = {
                "sample": prefix,
                "error": str(e),
                "merged_dna": None,
                "cropped_dna": None,
                "merge_info": {},
                "crop_info": {},
                "alignment_text": None,
                "identity_percent": 0.0,
                "detected_features": []
            }
            all_results.append(result)

    return all_results


def process_unidirectional(seq_files: List[Any], ref_file: Union[str, Any],
                           selected_categories: Optional[List[str]] = None) -> List[Dict[str, Any]]:
    """
    Process unidirectional sequences (single reads like M13R-only or pUC19-only).

    Uses local alignment since unidirectional sequences are unlikely to span
    the entire reference sequence. Tries both forward and reverse complement
    orientations and uses the better alignment.

    Args:
        seq_files: List of uploaded sequence files (file-like objects with .name attribute)
        ref_file: Excel file with reference sequences (file-like object or path)

    Returns:
        List of result dictionaries, one per processed file, each containing:
            - sample: Filename
            - merged_dna: None (no merging for unidirectional)
            - cropped_dna: The query DNA sequence (in best orientation)
            - merge_info: Empty dict
            - crop_info: Empty dict
            - reference: Reference sequence or None
            - alignment_text: Alignment visualization or None
            - identity_percent: Percent identity to reference
            - orientation: "forward" or "reverse complement"
    """
    if not seq_files:
        raise ValueError("No sequence files provided")

    # Load reference sequences
    ref_dict = load_reference_xlsx(ref_file)

    results = []
    for seq_file in seq_files:
        filename = seq_file.name

        try:
            # Read the sequence
            query_dna = read_seq_file(seq_file)

            # Check if this is an M13R read (sequences in reverse direction)
            # These need to be reverse-complemented before alignment
            is_reverse_primer = bool(re.search(r'm13r', filename, re.IGNORECASE))
            if is_reverse_primer:
                query_dna = str(Seq(query_dna).reverse_complement())

            # Match by filename prefix - try multiple strategies
            # Remove extension first
            base_name = filename.rsplit(".", 1)[0].lower()
            base_name = re.sub(r'\s*\(\d+\)$', '', base_name)  # strip Windows dup suffix e.g. " (1)"
            # Remove common suffixes like -pUC19, -M13R, etc.
            clean_name = re.sub(r'[-_]?(puc19(seqback)?|m13r|m13f)$', '', base_name, flags=re.IGNORECASE)
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

            # Detect features on both orientations
            features_fwd = detect_features(query_dna, selected_categories=selected_categories)
            query_dna_rc = str(Seq(query_dna).reverse_complement())
            features_rc = detect_features(query_dna_rc, selected_categories=selected_categories)

            # Initialize results
            # Always try both orientations and use the better alignment
            alignment_text = None
            identity_percent = 0.0
            # Track orientation - note if auto-flipped due to primer type
            orientation = "M13R (auto-flipped)" if is_reverse_primer else "forward"
            final_dna = query_dna
            detected_features = features_fwd

            if ref_seq:
                # Try forward orientation (local alignment)
                alignment_text_fwd, identity_fwd = align_dna_local(
                    query_dna, ref_seq, filename, features=features_fwd
                )

                # Try reverse complement
                alignment_text_rc, identity_rc = align_dna_local(
                    query_dna_rc, ref_seq, f"{filename} (reverse complement)", features=features_rc
                )

                # Use the orientation with better identity
                if identity_rc > identity_fwd:
                    alignment_text = alignment_text_rc
                    identity_percent = identity_rc
                    # Update orientation label based on original primer type
                    orientation = "forward (flipped back)" if is_reverse_primer else "reverse complement"
                    final_dna = query_dna_rc
                    detected_features = features_rc
                else:
                    alignment_text = alignment_text_fwd
                    identity_percent = identity_fwd

            results.append({
                "sample": filename,
                "merged_dna": None,
                "cropped_dna": final_dna,
                "merge_info": {},
                "crop_info": {},
                "reference": ref_seq,
                "matched_ref_name": matched_ref_name,
                "orientation": orientation,
                "available_refs": list(ref_dict.keys()) if not ref_seq else None,
                "tried_prefixes": candidates if not ref_seq else None,
                "alignment_text": alignment_text,
                "identity_percent": identity_percent,
                "detected_features": detected_features
            })

        except Exception as e:
            results.append({
                "sample": filename,
                "error": str(e),
                "merged_dna": None,
                "cropped_dna": None,
                "merge_info": {},
                "crop_info": {},
                "reference": None,
                "alignment_text": None,
                "identity_percent": 0.0,
                "detected_features": []
            })

    return results