import React, { useState } from "react";
import { submitAlign, type AlignResponse } from "../api/client";

export default function AlignTab() {
  const [seq1, setSeq1] = useState("");
  const [seq2, setSeq2] = useState("");
  const [seqType, setSeqType] = useState<"dna" | "protein">("dna");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<AlignResponse | null>(null);
  const [showOrfs, setShowOrfs] = useState(false);

  async function handleAlign() {
    setError(null);
    if (!seq1.trim() || !seq2.trim()) {
      setError("Paste both sequences before aligning.");
      return;
    }
    setLoading(true);
    try {
      const res = await submitAlign(seq1, seq2, seqType);
      setResult(res);
      setShowOrfs(false);
    } catch (e) {
      setError(e instanceof Error ? e.message : "Unknown error");
    } finally {
      setLoading(false);
    }
  }

  const html =
    showOrfs && result?.alignment_html_with_orfs
      ? result.alignment_html_with_orfs
      : result?.alignment_html;

  return (
    <div className="align-tab">
      <div className="form-section">
        <fieldset className="form-group">
          <legend>Sequence type</legend>
          <div className="radio-row">
            {(["dna", "protein"] as const).map((t) => (
              <label key={t} className="radio-label">
                <input
                  type="radio"
                  name="align-seq-type"
                  value={t}
                  checked={seqType === t}
                  onChange={() => {
                    setSeqType(t);
                    setResult(null);
                  }}
                />
                {t === "dna" ? "DNA" : "Amino acid"}
              </label>
            ))}
          </div>
        </fieldset>

        <div className="align-inputs">
          <div className="align-input-group">
            <label className="file-label">Sequence 1</label>
            <textarea
              className="align-textarea"
              rows={6}
              placeholder={
                seqType === "dna"
                  ? "Paste DNA sequence or FASTA..."
                  : "Paste amino acid sequence or FASTA..."
              }
              value={seq1}
              onChange={(e) => setSeq1(e.target.value)}
            />
          </div>
          <div className="align-input-group">
            <label className="file-label">Sequence 2</label>
            <textarea
              className="align-textarea"
              rows={6}
              placeholder={
                seqType === "dna"
                  ? "Paste DNA sequence or FASTA..."
                  : "Paste amino acid sequence or FASTA..."
              }
              value={seq2}
              onChange={(e) => setSeq2(e.target.value)}
            />
          </div>
        </div>

        {error && <div className="error-banner">{error}</div>}

        <button
          className="process-btn"
          disabled={loading}
          onClick={handleAlign}
        >
          {loading ? "Aligning..." : "Align"}
        </button>
      </div>

      {result && (
        <div className="align-results">
          {/* summary cards */}
          <div className="panel">
            <div className="panel-body">
              <div className="metrics-row">
                <div className="metric-card">
                  <span className="metric-value">
                    {result.identity_percent.toFixed(1)}%
                  </span>
                  <span className="metric-label">Identity</span>
                </div>
                {result.orf_count > 0 && (
                  <div className="metric-card">
                    <span className="metric-value">{result.orf_count}</span>
                    <span className="metric-label">ORFs detected</span>
                  </div>
                )}
              </div>
            </div>
          </div>

          {/* alignment panel */}
          <div className="panel">
            <div
              className="panel-header"
              style={{
                display: "flex",
                alignItems: "center",
                justifyContent: "space-between",
              }}
            >
              <span>Alignment</span>
              {result.orf_count > 0 && (
                <label
                  style={{
                    fontWeight: "normal",
                    fontSize: "0.85rem",
                    display: "flex",
                    alignItems: "center",
                    gap: "0.3rem",
                  }}
                >
                  <input
                    type="checkbox"
                    checked={showOrfs}
                    onChange={(e) => setShowOrfs(e.target.checked)}
                  />
                  Show ORFs ({result.orf_count} detected)
                </label>
              )}
            </div>
            <div className="panel-body alignment-body">
              {html ? (
                <div dangerouslySetInnerHTML={{ __html: html }} />
              ) : (
                <p className="muted">No alignment available.</p>
              )}
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
