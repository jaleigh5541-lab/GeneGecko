import re
from typing import List, Tuple

from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
from pydantic import BaseModel

from Bio import pairwise2
from Bio.Seq import Seq

router = APIRouter(prefix="/api")


# ── helpers ──────────────────────────────────────────────────────

def _longest_orf(aa_seq: str) -> str:
    """Longest M→* (or M→end) ORF in an amino acid sequence."""
    orfs = re.findall(r"M[^*]*", aa_seq)
    if orfs:
        return max(orfs, key=len)
    m = aa_seq.find("M")
    return aa_seq[m:] if m != -1 else ""


def _find_all_orfs(aa_seq: str, min_length: int = 1) -> List[Tuple[str, int]]:
    """All non-overlapping ORFs (M→* or M→end) >= min_length."""
    orfs, pos = [], 0
    while pos < len(aa_seq):
        m = aa_seq.find("M", pos)
        if m == -1:
            break
        stop = aa_seq.find("*", m)
        orf = aa_seq[m:stop] if stop != -1 else aa_seq[m:]
        if len(orf) >= min_length:
            orfs.append((orf, m))
        pos = stop + 1 if stop != -1 else len(aa_seq)
    return orfs


def _ungapped_to_gapped(gapped_seq: str) -> List[int]:
    """Map ungapped index → gapped position."""
    return [gi for gi, c in enumerate(gapped_seq) if c != "-"]


def _clean_seq(text: str, alphabet: str) -> str:
    """Strip whitespace / FASTA headers and keep only valid characters."""
    lines = text.strip().splitlines()
    out = []
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            continue
        out.append(line)
    raw = "".join(out).upper()
    return re.sub(f"[^{alphabet}]", "", raw)


# ── alignment + HTML builders ────────────────────────────────────

def _align_dna(seq1: str, seq2: str, block_size: int = 60):
    """Align two DNA sequences, translate, find ORFs, produce HTML."""
    alns = pairwise2.align.localms(seq1, seq2, 2, -1, -2, -0.5)
    if not alns:
        empty = "<p>No alignment found.</p>"
        return {
            "alignment_html": empty,
            "alignment_html_with_orfs": empty,
            "identity_percent": 0.0,
            "orf_count": 0,
            "longest_orf_length": 0,
            "longest_orf_seq": "",
        }

    aln = alns[0]
    seqA_full, seqB_full = aln[0], aln[1]
    begin, end = int(aln[3]), int(aln[4])
    seqA_reg = seqA_full[begin:end]
    seqB_reg = seqB_full[begin:end]

    # trim to reference-covered region (no leading/trailing gaps in seq2)
    r_start = next((i for i, c in enumerate(seqB_reg) if c != "-"), 0)
    r_end = len(seqB_reg) - next(
        (i for i, c in enumerate(reversed(seqB_reg)) if c != "-"), 0
    )
    seqA_str = seqA_reg[r_start:r_end]
    seqB_str = seqB_reg[r_start:r_end]

    q_offset = len(seqA_full[:begin].replace("-", "")) + len(
        seqA_reg[:r_start].replace("-", "")
    )

    # identity
    matches = sum(a == b and a != "-" for a, b in zip(seqA_str, seqB_str))
    ref_len = len(seqB_str.replace("-", ""))
    identity = (matches / ref_len) * 100 if ref_len else 0.0

    # translate
    q_nogap = seqA_str.replace("-", "")
    r_nogap = seqB_str.replace("-", "")
    q_trim = len(q_nogap) - len(q_nogap) % 3
    r_trim = len(r_nogap) - len(r_nogap) % 3
    q_aa = str(Seq(q_nogap[:q_trim]).translate()) if q_trim else ""
    r_aa = str(Seq(r_nogap[:r_trim]).translate()) if r_trim else ""

    # ORFs (no min length)
    q_orfs = _find_all_orfs(q_aa)
    r_orfs = _find_all_orfs(r_aa)
    total_orfs = len(q_orfs) + len(r_orfs)

    longest = _longest_orf(q_aa)
    longest_r = _longest_orf(r_aa)
    if len(longest_r) > len(longest):
        longest = longest_r

    u2g_q = _ungapped_to_gapped(seqA_str)
    u2g_r = _ungapped_to_gapped(seqB_str)

    # ── build HTML ───
    aln_len = len(seqA_str)
    css = (
        "<style>"
        ".al-m{color:#155724;font-weight:bold;background:#d4edda}"
        ".al-x{color:#721c24;font-weight:bold;background:#f8d7da}"
        ".al-g{color:#856404;background:#fff3cd}"
        ".al-p{color:#2a9d8f}"
        ".al-o{color:#2a6496;font-weight:bold}"
        "</style>"
    )
    header = (
        f"<div style='font-family:monospace;font-size:14px;line-height:1.6;white-space:pre'>"
        f"<div style='font-weight:bold;margin-bottom:10px;white-space:normal;color:#264653'>"
        f"DNA Alignment — {aln_len} bp aligned region | Identity: {identity:.1f}%</div>"
    )
    html_parts = [css, header]
    orf_parts = [css, header] if total_orfs > 0 else []

    lw = 22
    q_ctr = q_offset + 1
    r_ctr = 1

    def _orf_rows(orfs, u2g, blk_start, bsize, prefix):
        rows = []
        for idx, (orf_aa, orf_start_aa) in enumerate(orfs):
            chars = [" "] * bsize
            for ai, ac in enumerate(orf_aa):
                mid = orf_start_aa * 3 + ai * 3 + 1
                if mid >= len(u2g):
                    break
                col = u2g[mid] - blk_start
                if 0 <= col < bsize:
                    chars[col] = ac
            label = f"{prefix}{idx+1}:".ljust(lw)
            rows.append(
                f"<div class='gi-orf-row'>{label}<span class='al-o'>{''.join(chars)}</span></div>"
            )
        return rows

    for i in range(0, len(seqA_str), block_size):
        bA = seqA_str[i : i + block_size]
        bB = seqB_str[i : i + block_size]
        qe = q_ctr + sum(1 for c in bA if c != "-") - 1
        re_ = r_ctr + sum(1 for c in bB if c != "-") - 1
        ql = f"Seq1   {q_ctr}-{qe}:".ljust(lw)
        rl = f"Seq2   {r_ctr}-{re_}:".ljust(lw)

        qh, rh, mh = "", "", ""
        for a, b in zip(bA, bB):
            if a == "-":
                qh += f"<span class='al-g'>{a}</span>"
            elif a == b:
                qh += f"<span class='al-m'>{a}</span>"
            else:
                qh += f"<span class='al-x'>{a}</span>"
            if b == "-":
                rh += f"<span class='al-g'>{b}</span>"
            elif a == b:
                rh += f"<span class='al-m'>{b}</span>"
            else:
                rh += f"<span class='al-x'>{b}</span>"
            mh += "<span class='al-p'>|</span>" if a == b and a != "-" else " "

        blk = (
            f"<div style='margin-bottom:12px'>"
            f"<div>{ql}{qh}</div>"
            f"<div>{' ' * lw}{mh}</div>"
            f"<div>{rl}{rh}</div>"
            f"</div>"
        )
        html_parts.append(blk)

        if total_orfs > 0:
            bsize = len(bA)
            orf_parts.append("<div style='margin-bottom:12px'>")
            orf_parts.extend(_orf_rows(q_orfs, u2g_q, i, bsize, "S1-ORF"))
            orf_parts.append(f"<div>{ql}{qh}</div>")
            orf_parts.append(f"<div>{' ' * lw}{mh}</div>")
            orf_parts.append(f"<div>{rl}{rh}</div>")
            orf_parts.extend(_orf_rows(r_orfs, u2g_r, i, bsize, "S2-ORF"))
            orf_parts.append("</div>")

        q_ctr += sum(1 for c in bA if c != "-")
        r_ctr += sum(1 for c in bB if c != "-")

    # protein alignment section (always shown if longest ORF exists)
    if longest and q_aa and r_aa:
        prot_alns = pairwise2.align.globalms(q_aa, r_aa, 2, -1, -2, -0.5)
        if prot_alns:
            pa = prot_alns[0]
            pA, pB = pa[0], pa[1]
            pm = sum(a == b and a != "-" for a, b in zip(pA, pB))
            pl = max(len(pA), len(pB))
            pi = (pm / pl) * 100 if pl else 0.0

            sep = "<hr style='margin:20px 0;border:1px solid #ccc'>"
            prot_hdr = (
                f"<div style='font-weight:bold;margin-bottom:10px;white-space:normal'>"
                f"Translated Protein Alignment — Identity: {pi:.1f}%</div>"
            )
            html_parts.append(sep + prot_hdr)
            if total_orfs > 0:
                orf_parts.append(sep + prot_hdr)

            pq, pr = 1, 1
            for j in range(0, len(pA), block_size):
                bpA = pA[j : j + block_size]
                bpB = pB[j : j + block_size]
                pqe = pq + sum(1 for c in bpA if c != "-") - 1
                pre = pr + sum(1 for c in bpB if c != "-") - 1
                pql = f"Seq1   {pq}-{pqe}:".ljust(lw)
                prl = f"Seq2   {pr}-{pre}:".ljust(lw)
                pqh, prh, pmh = "", "", ""
                for a, b in zip(bpA, bpB):
                    if a == "-":
                        pqh += f"<span class='al-g'>{a}</span>"
                    elif a == b:
                        pqh += f"<span class='al-m'>{a}</span>"
                    else:
                        pqh += f"<span class='al-x'>{a}</span>"
                    if b == "-":
                        prh += f"<span class='al-g'>{b}</span>"
                    elif a == b:
                        prh += f"<span class='al-m'>{b}</span>"
                    else:
                        prh += f"<span class='al-x'>{b}</span>"
                    pmh += "<span class='al-p'>|</span>" if a == b and a != "-" else " "
                pblk = (
                    f"<div style='margin-bottom:12px'>"
                    f"<div>{pql}{pqh}</div>"
                    f"<div>{' ' * lw}{pmh}</div>"
                    f"<div>{prl}{prh}</div></div>"
                )
                html_parts.append(pblk)
                if total_orfs > 0:
                    orf_parts.append(pblk)
                pq += sum(1 for c in bpA if c != "-")
                pr += sum(1 for c in bpB if c != "-")

    html_parts.append("</div>")
    if total_orfs > 0:
        orf_parts.append("</div>")

    return {
        "alignment_html": "".join(html_parts),
        "alignment_html_with_orfs": "".join(orf_parts) if total_orfs > 0 else None,
        "identity_percent": round(identity, 2),
        "orf_count": total_orfs,
        "longest_orf_length": len(longest),
        "longest_orf_seq": longest,
    }


def _align_protein(seq1: str, seq2: str, block_size: int = 60):
    """Align two amino acid sequences, find ORFs, produce HTML."""
    alns = pairwise2.align.globalms(seq1, seq2, 2, -1, -2, -0.5)
    if not alns:
        empty = "<p>No alignment found.</p>"
        return {
            "alignment_html": empty,
            "alignment_html_with_orfs": empty,
            "identity_percent": 0.0,
            "orf_count": 0,
            "longest_orf_length": 0,
            "longest_orf_seq": "",
        }

    aln = alns[0]
    seqA, seqB = aln[0], aln[1]

    matches = sum(a == b and a != "-" for a, b in zip(seqA, seqB))
    aln_len = max(len(seqA), len(seqB))
    identity = (matches / aln_len) * 100 if aln_len else 0.0

    # ORFs from ungapped sequences
    s1_nogap = seqA.replace("-", "")
    s2_nogap = seqB.replace("-", "")
    s1_orfs = _find_all_orfs(s1_nogap)
    s2_orfs = _find_all_orfs(s2_nogap)
    total_orfs = len(s1_orfs) + len(s2_orfs)

    longest = _longest_orf(s1_nogap)
    longest2 = _longest_orf(s2_nogap)
    if len(longest2) > len(longest):
        longest = longest2

    u2g_1 = _ungapped_to_gapped(seqA)
    u2g_2 = _ungapped_to_gapped(seqB)

    css = (
        "<style>"
        ".al-m{color:#155724;font-weight:bold;background:#d4edda}"
        ".al-x{color:#721c24;font-weight:bold;background:#f8d7da}"
        ".al-g{color:#856404;background:#fff3cd}"
        ".al-p{color:#2a9d8f}"
        ".al-o{color:#2a6496;font-weight:bold}"
        "</style>"
    )
    header = (
        f"<div style='font-family:monospace;font-size:14px;line-height:1.6;white-space:pre'>"
        f"<div style='font-weight:bold;margin-bottom:10px;white-space:normal;color:#264653'>"
        f"Protein Alignment — {aln_len} positions | Identity: {identity:.1f}%</div>"
        f"<div style='margin-bottom:10px;white-space:normal;color:#2a9d8f'>"
        f"Seq 1: {len(s1_nogap)} aa | Seq 2: {len(s2_nogap)} aa"
        f"{' | Longest ORF: ' + str(len(longest)) + ' aa' if longest else ''}</div>"
    )
    html_parts = [css, header]
    orf_parts = [css, header] if total_orfs > 0 else []

    lw = 22
    q_ctr, r_ctr = 1, 1

    def _aa_orf_rows(orfs, u2g, blk_start, bsize, prefix):
        rows = []
        for idx, (orf_aa, orf_start) in enumerate(orfs):
            chars = [" "] * bsize
            for ai, ac in enumerate(orf_aa):
                ungapped = orf_start + ai
                if ungapped >= len(u2g):
                    break
                col = u2g[ungapped] - blk_start
                if 0 <= col < bsize:
                    chars[col] = ac
            label = f"{prefix}{idx+1}:".ljust(lw)
            rows.append(
                f"<div class='gi-orf-row'>{label}<span class='al-o'>{''.join(chars)}</span></div>"
            )
        return rows

    for i in range(0, len(seqA), block_size):
        bA = seqA[i : i + block_size]
        bB = seqB[i : i + block_size]
        qe = q_ctr + sum(1 for c in bA if c != "-") - 1
        re_ = r_ctr + sum(1 for c in bB if c != "-") - 1
        ql = f"Seq1   {q_ctr}-{qe}:".ljust(lw)
        rl = f"Seq2   {r_ctr}-{re_}:".ljust(lw)

        qh, rh, mh = "", "", ""
        for a, b in zip(bA, bB):
            if a == "-":
                qh += f"<span class='al-g'>{a}</span>"
            elif a == b:
                qh += f"<span class='al-m'>{a}</span>"
            else:
                qh += f"<span class='al-x'>{a}</span>"
            if b == "-":
                rh += f"<span class='al-g'>{b}</span>"
            elif a == b:
                rh += f"<span class='al-m'>{b}</span>"
            else:
                rh += f"<span class='al-x'>{b}</span>"
            mh += "<span class='al-p'>|</span>" if a == b and a != "-" else " "

        blk = (
            f"<div style='margin-bottom:12px'>"
            f"<div>{ql}{qh}</div>"
            f"<div>{' ' * lw}{mh}</div>"
            f"<div>{rl}{rh}</div></div>"
        )
        html_parts.append(blk)

        if total_orfs > 0:
            bsize = len(bA)
            orf_parts.append("<div style='margin-bottom:12px'>")
            orf_parts.extend(_aa_orf_rows(s1_orfs, u2g_1, i, bsize, "S1-ORF"))
            orf_parts.append(f"<div>{ql}{qh}</div>")
            orf_parts.append(f"<div>{' ' * lw}{mh}</div>")
            orf_parts.append(f"<div>{rl}{rh}</div>")
            orf_parts.extend(_aa_orf_rows(s2_orfs, u2g_2, i, bsize, "S2-ORF"))
            orf_parts.append("</div>")

        q_ctr += sum(1 for c in bA if c != "-")
        r_ctr += sum(1 for c in bB if c != "-")

    html_parts.append("</div>")
    if total_orfs > 0:
        orf_parts.append("</div>")

    return {
        "alignment_html": "".join(html_parts),
        "alignment_html_with_orfs": "".join(orf_parts) if total_orfs > 0 else None,
        "identity_percent": round(identity, 2),
        "orf_count": total_orfs,
        "longest_orf_length": len(longest),
        "longest_orf_seq": longest,
    }


# ── endpoint ─────────────────────────────────────────────────────

class AlignRequest(BaseModel):
    seq1: str
    seq2: str
    seq_type: str = "dna"  # "dna" or "protein"


@router.post("/align")
async def align(req: AlignRequest) -> JSONResponse:
    if req.seq_type == "dna":
        s1 = _clean_seq(req.seq1, "ACGTN")
        s2 = _clean_seq(req.seq2, "ACGTN")
    else:
        s1 = _clean_seq(req.seq1, "ACDEFGHIKLMNPQRSTVWY*")
        s2 = _clean_seq(req.seq2, "ACDEFGHIKLMNPQRSTVWY*")

    if len(s1) < 3 or len(s2) < 3:
        raise HTTPException(status_code=400, detail="Both sequences must be at least 3 characters.")

    try:
        if req.seq_type == "dna":
            result = _align_dna(s1, s2)
        else:
            result = _align_protein(s1, s2)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e

    return JSONResponse({"success": True, **result})
