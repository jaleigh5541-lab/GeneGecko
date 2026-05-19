import React, { useState } from "react";

interface Result {
  sample?: string;
  merged_dna?: string;
  cropped_dna?: string;
  identity_percent?: number;
  longest_orf?: string;
  orf_dna?: string;
  orientation?: string;
  frame?: string;
  error?: string;
}

function escapeCsvField(value: string): string {
  if (value.includes(",") || value.includes('"') || value.includes("\n")) {
    return `"${value.replace(/"/g, '""')}"`;
  }
  return value;
}

function downloadCsv(results: Result[], sequenceType: string) {
  const header = ["Sample", "DNA Length", "% Identity", "ORF Length (aa)", "DNA Contig Sequence", "ORF DNA Sequence", "ORF Amino Acid Sequence"];
  const rows = results.map((r) => {
    const dna = r.merged_dna ?? r.cropped_dna ?? "";
    const dnaLen = dna.length || "";
    const identity = r.error ? "" : (r.identity_percent ?? 0).toFixed(1);
    const orf = r.longest_orf ?? "";
    const orfLen = orf.length > 70 ? orf.length : "";
    const orfDna = orf.length > 70 ? (r.orf_dna ?? "") : "";
    return [
      r.sample ?? "",
      String(dnaLen),
      identity,
      String(orfLen),
      dna,
      orfDna,
      orf.length > 70 ? orf : "",
    ].map(escapeCsvField).join(",");
  });
  const csv = [header.join(","), ...rows].join("\n");
  const blob = new Blob([csv], { type: "text/csv" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = "summary.csv";
  a.click();
  URL.revokeObjectURL(url);
}

interface Props {
  results: Result[];
  sequenceType: string;
  primerType: string;
  onSelectSample: (index: number) => void;
}

export default function SummaryTable({ results, sequenceType, primerType, onSelectSample }: Props) {
  const [open, setOpen] = useState(true);

  const isProtein = primerType === "protein";
  const isNonPaired = sequenceType !== "paired";

  function getDnaLength(r: Result): number | string {
    if (!isProtein && sequenceType === "paired") {
      return (r.merged_dna ?? "").length || "—";
    }
    return (r.merged_dna ?? r.cropped_dna ?? "").length || "—";
  }

  function getOrfLength(r: Result): number | null {
    const len = (r.longest_orf ?? "").length;
    return len > 70 ? len : null;
  }

  function getIdentityClass(r: Result): string {
    if (r.error) return "st-error";
    if ((r.identity_percent ?? 0) >= 95) return "st-high";
    if ((r.identity_percent ?? 0) >= 80) return "st-mid";
    return "st-low";
  }

  return (
    <div className="panel summary-table-panel">
      <div className="panel-header" style={{ display: "flex", alignItems: "center", justifyContent: "space-between" }}>
        <button className="panel-header-toggle" onClick={() => setOpen(!open)} style={{ background: "none", border: "none", color: "inherit", font: "inherit", cursor: "pointer", padding: 0 }}>
          {open ? "▼" : "▶"} All Samples Summary ({results.length})
        </button>
        <button
          className="dl-btn"
          style={{ fontSize: "0.8rem", padding: "4px 10px" }}
          onClick={(e) => { e.stopPropagation(); downloadCsv(results, sequenceType); }}
          title="Download summary as CSV"
        >
          Download CSV
        </button>
      </div>
      {open && (
        <div className="panel-body" style={{ padding: 0 }}>
          <div className="summary-table-wrap">
            <table className="summary-table">
              <thead>
                <tr>
                  <th>Sample</th>
                  <th>DNA Length</th>
                  {!isProtein && sequenceType === "paired" && <th>Cropped Length</th>}
                  <th>% Identity</th>
                  {isProtein && <th>ORF Length (aa)</th>}
                  {isNonPaired && <th>Orientation</th>}
                  {isProtein && isNonPaired && <th>Frame</th>}
                </tr>
              </thead>
              <tbody>
                {results.map((r, i) => {
                  const orfLen = getOrfLength(r);
                  return (
                    <tr
                      key={i}
                      className="summary-table-row"
                      onClick={() => onSelectSample(i)}
                      title="Click to view this sample"
                    >
                      <td className="st-sample">{r.sample}</td>
                      <td>{getDnaLength(r)}</td>
                      {!isProtein && sequenceType === "paired" && (
                        <td>{(r.cropped_dna ?? "").length || "—"}</td>
                      )}
                      <td>
                        <span className={`st-identity ${getIdentityClass(r)}`}>
                          {r.error ? "—" : `${(r.identity_percent ?? 0).toFixed(1)}%`}
                        </span>
                      </td>
                      {isProtein && (
                        <td>
                          {orfLen != null ? (
                            orfLen
                          ) : (
                            <span className="st-badge-dim">{"< 70"}</span>
                          )}
                        </td>
                      )}
                      {isNonPaired && (
                        <td>
                          <span className="st-orientation">
                            {r.orientation ?? "forward"}
                          </span>
                        </td>
                      )}
                      {isProtein && isNonPaired && <td>{r.frame ?? "—"}</td>}
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>
        </div>
      )}
    </div>
  );
}
