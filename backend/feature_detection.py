#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Feature Detection Module for PA-app2

Detects common DNA features (promoters, tags, terminators, etc.) in sequences
and provides position mapping for alignment visualization.
"""

import re
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from Bio.Seq import Seq


# ───────────────────────────────────────────────────────────────
#  Data Structures
# ───────────────────────────────────────────────────────────────

@dataclass
class FeatureDefinition:
    """Definition of a DNA feature to search for."""
    name: str           # Display name, e.g., "T7 Promoter"
    sequence: str       # DNA sequence (5' to 3')
    category: str       # Category, e.g., "Promoter", "Tag"
    color: str          # Hex color for highlighting
    description: str = ""  # Optional description


@dataclass
class DetectedFeature:
    """A feature found in a sequence."""
    name: str           # Feature name
    category: str       # Feature category
    start: int          # 1-indexed start position in raw sequence
    end: int            # 1-indexed end position in raw sequence
    strand: str         # "+" (forward) or "-" (reverse complement)
    sequence: str       # Actual sequence matched
    color: str          # Color for highlighting
    gapped_start: Optional[int] = None  # Position in gapped alignment
    gapped_end: Optional[int] = None    # Position in gapped alignment


# ───────────────────────────────────────────────────────────────
#  Color Scheme
# ───────────────────────────────────────────────────────────────

CATEGORY_COLORS = {
    "Promoter": "#3b82f6",       # Blue
    "Terminator": "#ef4444",     # Red
    "Tag": "#8b5cf6",            # Purple
    "RBS": "#22c55e",            # Green
    "Restriction Site": "#f59e0b",  # Amber
    "Linker": "#06b6d4",         # Cyan
    "Cleavage Site": "#ec4899",  # Pink
    "Custom": "#6b7280",         # Gray
}


# ───────────────────────────────────────────────────────────────
#  Built-in Feature Database
# ───────────────────────────────────────────────────────────────

BUILTIN_FEATURES: Dict[str, FeatureDefinition] = {
    # --- Promoters ---
    "t7_promoter": FeatureDefinition(
        name="T7 Promoter",
        sequence="TAATACGACTCACTATAGG",
        category="Promoter",
        color=CATEGORY_COLORS["Promoter"],
        description="T7 RNA polymerase promoter"
    ),
    "t7_promoter_alt": FeatureDefinition(
        name="T7 Promoter (17bp)",
        sequence="TAATACGACTCACTATAG",
        category="Promoter",
        color=CATEGORY_COLORS["Promoter"],
        description="T7 RNA polymerase promoter (shorter variant)"
    ),
    "sp6_promoter": FeatureDefinition(
        name="SP6 Promoter",
        sequence="ATTTAGGTGACACTATAG",
        category="Promoter",
        color=CATEGORY_COLORS["Promoter"],
        description="SP6 RNA polymerase promoter"
    ),
    "t3_promoter": FeatureDefinition(
        name="T3 Promoter",
        sequence="AATTAACCCTCACTAAAG",
        category="Promoter",
        color=CATEGORY_COLORS["Promoter"],
        description="T3 RNA polymerase promoter"
    ),

    # --- Terminators ---
    "t7_terminator": FeatureDefinition(
        name="T7 Terminator",
        sequence="GCTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG",
        category="Terminator",
        color=CATEGORY_COLORS["Terminator"],
        description="T7 transcription terminator"
    ),
    "t7_terminator_short": FeatureDefinition(
        name="T7 Term (short)",
        sequence="CTAGCATAACCCCTTGGGGCCTCT",
        category="Terminator",
        color=CATEGORY_COLORS["Terminator"],
        description="T7 terminator core sequence"
    ),

    # --- Ribosome Binding Sites ---
    "kozak": FeatureDefinition(
        name="Kozak",
        sequence="GCCGCCACCATGG",
        category="RBS",
        color=CATEGORY_COLORS["RBS"],
        description="Kozak consensus sequence for translation initiation"
    ),
    "kozak_minimal": FeatureDefinition(
        name="Kozak (minimal)",
        sequence="ACCATGG",
        category="RBS",
        color=CATEGORY_COLORS["RBS"],
        description="Minimal Kozak sequence"
    ),
    "shine_dalgarno": FeatureDefinition(
        name="Shine-Dalgarno",
        sequence="AGGAGG",
        category="RBS",
        color=CATEGORY_COLORS["RBS"],
        description="Prokaryotic ribosome binding site"
    ),

    # --- Affinity Tags ---
    "his_tag_cat": FeatureDefinition(
        name="6x His-tag (CAT)",
        sequence="CATCATCATCATCATCAT",
        category="Tag",
        color=CATEGORY_COLORS["Tag"],
        description="Hexahistidine tag (CAT codons)"
    ),
    "his_tag_cac": FeatureDefinition(
        name="6x His-tag (CAC)",
        sequence="CACCACCACCACCACCAC",
        category="Tag",
        color=CATEGORY_COLORS["Tag"],
        description="Hexahistidine tag (CAC codons)"
    ),
    "flag_tag": FeatureDefinition(
        name="FLAG tag",
        sequence="GACTACAAGGACGACGATGACAAG",
        category="Tag",
        color=CATEGORY_COLORS["Tag"],
        description="FLAG epitope tag (DYKDDDDK)"
    ),
    "ha_tag": FeatureDefinition(
        name="HA tag",
        sequence="TACCCATACGATGTTCCAGATTACGCT",
        category="Tag",
        color=CATEGORY_COLORS["Tag"],
        description="Hemagglutinin epitope tag"
    ),
    "myc_tag": FeatureDefinition(
        name="Myc tag",
        sequence="GAACAAAAACTCATCTCAGAAGAGGATCTG",
        category="Tag",
        color=CATEGORY_COLORS["Tag"],
        description="c-Myc epitope tag"
    ),
    "strep_tag": FeatureDefinition(
        name="Strep-tag II",
        sequence="TGGTCCCACCCGCAGTTCGAAAAA",
        category="Tag",
        color=CATEGORY_COLORS["Tag"],
        description="Strep-tag II for streptavidin purification"
    ),

    # --- Restriction Sites ---
    "ndei": FeatureDefinition(
        name="NdeI",
        sequence="CATATG",
        category="Restriction Site",
        color=CATEGORY_COLORS["Restriction Site"],
        description="NdeI restriction site"
    ),
    "ncoi": FeatureDefinition(
        name="NcoI",
        sequence="CCATGG",
        category="Restriction Site",
        color=CATEGORY_COLORS["Restriction Site"],
        description="NcoI restriction site"
    ),
    "bamhi": FeatureDefinition(
        name="BamHI",
        sequence="GGATCC",
        category="Restriction Site",
        color=CATEGORY_COLORS["Restriction Site"],
        description="BamHI restriction site"
    ),
    "xhoi": FeatureDefinition(
        name="XhoI",
        sequence="CTCGAG",
        category="Restriction Site",
        color=CATEGORY_COLORS["Restriction Site"],
        description="XhoI restriction site"
    ),
    "hindiii": FeatureDefinition(
        name="HindIII",
        sequence="AAGCTT",
        category="Restriction Site",
        color=CATEGORY_COLORS["Restriction Site"],
        description="HindIII restriction site"
    ),
    "ecori": FeatureDefinition(
        name="EcoRI",
        sequence="GAATTC",
        category="Restriction Site",
        color=CATEGORY_COLORS["Restriction Site"],
        description="EcoRI restriction site"
    ),
    "noti": FeatureDefinition(
        name="NotI",
        sequence="GCGGCCGC",
        category="Restriction Site",
        color=CATEGORY_COLORS["Restriction Site"],
        description="NotI restriction site"
    ),
    "sali": FeatureDefinition(
        name="SalI",
        sequence="GTCGAC",
        category="Restriction Site",
        color=CATEGORY_COLORS["Restriction Site"],
        description="SalI restriction site"
    ),

    # --- Linkers ---
    "gs_linker": FeatureDefinition(
        name="GS Linker",
        sequence="GGTTCTGGTGGTTCTGGTGGTTCT",
        category="Linker",
        color=CATEGORY_COLORS["Linker"],
        description="Glycine-Serine flexible linker"
    ),

    # --- Cleavage Sites ---
    "tev_site": FeatureDefinition(
        name="TEV site",
        sequence="GAAAACCTGTATTTTCAGGGC",
        category="Cleavage Site",
        color=CATEGORY_COLORS["Cleavage Site"],
        description="TEV protease cleavage site (ENLYFQG)"
    ),
    "enterokinase_site": FeatureDefinition(
        name="Enterokinase site",
        sequence="GATGATGATGATAAG",
        category="Cleavage Site",
        color=CATEGORY_COLORS["Cleavage Site"],
        description="Enterokinase cleavage site (DDDDK)"
    ),
}


# ───────────────────────────────────────────────────────────────
#  Feature Detection Functions
# ───────────────────────────────────────────────────────────────

def detect_features(
    sequence: str,
    features_to_search: Optional[Dict[str, FeatureDefinition]] = None,
    search_reverse_complement: bool = True,
    min_feature_length: int = 6,
    selected_categories: Optional[List[str]] = None
) -> List[DetectedFeature]:
    """
    Detect DNA features in a sequence.

    Args:
        sequence: DNA sequence to search (will be uppercased)
        features_to_search: Feature definitions to search for.
                           If None, uses BUILTIN_FEATURES.
        search_reverse_complement: Also search for reverse complement matches
        min_feature_length: Minimum feature length to consider
        selected_categories: Only search these categories (None = all)

    Returns:
        List of DetectedFeature objects, sorted by start position
    """
    if features_to_search is None:
        features_to_search = BUILTIN_FEATURES

    detected = []
    sequence = sequence.upper()

    for feature_id, feature_def in features_to_search.items():
        # Skip if category not selected
        if selected_categories and feature_def.category not in selected_categories:
            continue

        # Skip very short features
        if len(feature_def.sequence) < min_feature_length:
            continue

        pattern = feature_def.sequence.upper()

        # Search forward strand
        for match in re.finditer(re.escape(pattern), sequence):
            detected.append(DetectedFeature(
                name=feature_def.name,
                category=feature_def.category,
                start=match.start() + 1,  # 1-indexed
                end=match.end(),
                strand="+",
                sequence=match.group(),
                color=feature_def.color
            ))

        # Search reverse complement
        if search_reverse_complement:
            rc_pattern = str(Seq(pattern).reverse_complement())
            for match in re.finditer(re.escape(rc_pattern), sequence):
                detected.append(DetectedFeature(
                    name=feature_def.name,
                    category=feature_def.category,
                    start=match.start() + 1,
                    end=match.end(),
                    strand="-",
                    sequence=match.group(),
                    color=feature_def.color
                ))

    # Sort by start position
    detected.sort(key=lambda x: x.start)
    return detected


def map_features_to_gapped_alignment(
    features: List[DetectedFeature],
    gapped_sequence: str
) -> List[DetectedFeature]:
    """
    Map feature positions from ungapped sequence to gapped alignment.

    Args:
        features: List of detected features with ungapped positions
        gapped_sequence: Aligned sequence with gaps ('-')

    Returns:
        Features with updated gapped_start and gapped_end positions
    """
    # Build position mapping: ungapped_pos (0-indexed) -> gapped_pos
    ungapped_to_gapped = {}
    ungapped_pos = 0

    for gapped_pos, char in enumerate(gapped_sequence):
        if char != '-':
            ungapped_to_gapped[ungapped_pos] = gapped_pos
            ungapped_pos += 1

    # Map each feature
    for feature in features:
        # Convert from 1-indexed to 0-indexed for mapping
        start_0 = feature.start - 1
        end_0 = feature.end - 1

        if start_0 in ungapped_to_gapped:
            feature.gapped_start = ungapped_to_gapped[start_0]
        if end_0 in ungapped_to_gapped:
            feature.gapped_end = ungapped_to_gapped[end_0]

    return features


def parse_custom_features(
    text_input: str,
    default_category: str = "Custom"
) -> Dict[str, FeatureDefinition]:
    """
    Parse custom features from user text input.

    Expected format (one per line):
        FeatureName: SEQUENCE
        FeatureName, Category: SEQUENCE

    Lines starting with # are comments.

    Args:
        text_input: Multi-line text with feature definitions
        default_category: Category to use if not specified

    Returns:
        Dictionary of FeatureDefinition objects
    """
    custom_features = {}

    for line in text_input.strip().split('\n'):
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        if ':' not in line:
            continue

        name_part, sequence = line.split(':', 1)
        sequence = sequence.strip().upper()

        # Validate sequence (DNA only)
        if not re.match(r'^[ACGTN]+$', sequence):
            continue

        # Skip very short sequences
        if len(sequence) < 4:
            continue

        # Parse name and optional category
        if ',' in name_part:
            name, category = name_part.split(',', 1)
            name = name.strip()
            category = category.strip()
        else:
            name = name_part.strip()
            category = default_category

        # Get color for category
        color = CATEGORY_COLORS.get(category, CATEGORY_COLORS["Custom"])

        feature_id = name.lower().replace(' ', '_')
        custom_features[feature_id] = FeatureDefinition(
            name=name,
            sequence=sequence,
            category=category,
            color=color
        )

    return custom_features


def get_active_features(
    use_builtin: bool = True,
    custom_features: Optional[Dict[str, FeatureDefinition]] = None,
    selected_categories: Optional[List[str]] = None
) -> Dict[str, FeatureDefinition]:
    """
    Combine built-in and custom features based on user preferences.

    Args:
        use_builtin: Include built-in feature database
        custom_features: User-defined custom features
        selected_categories: Only include these categories (None = all)

    Returns:
        Combined dictionary of features to search
    """
    active = {}

    if use_builtin:
        active.update(BUILTIN_FEATURES)

    if custom_features:
        active.update(custom_features)

    # Filter by category if specified
    if selected_categories:
        active = {
            k: v for k, v in active.items()
            if v.category in selected_categories
        }

    return active


def get_all_categories() -> List[str]:
    """Return list of all available feature categories."""
    return list(CATEGORY_COLORS.keys())


def build_feature_position_map(
    features: List[DetectedFeature],
    use_gapped: bool = False
) -> Dict[int, DetectedFeature]:
    """
    Build a position-to-feature map for quick lookup during HTML generation.

    Args:
        features: List of detected features
        use_gapped: Use gapped positions instead of raw positions

    Returns:
        Dictionary mapping position -> feature (first feature at each position)
    """
    feature_map = {}

    for feature in features:
        if use_gapped:
            if feature.gapped_start is None or feature.gapped_end is None:
                continue
            start = feature.gapped_start
            end = feature.gapped_end
        else:
            start = feature.start - 1  # Convert to 0-indexed
            end = feature.end - 1

        for pos in range(start, end + 1):
            if pos not in feature_map:
                feature_map[pos] = feature

    return feature_map
