import React, { useState } from "react";

function downloadBlob(content: string, filename: string) {
  const blob = new Blob([content], { type: "text/plain" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

interface Result {
  sample?: string;
  merged_dna?: string;
  cropped_dna?: string;
  longest_orf?: string;
}

interface Props {
  result: Result;
  sequenceType: string;
  primerType: string;
}

export default function DownloadPanel({ result, sequenceType, primerType }: Props) {
  const [open, setOpen] = useState(false);
  const r = result;
  const isProtein = primerType === "protein";

  return (
    <div className="panel">
      <button className="panel-header" onClick={() => setOpen(!open)}>
        {open ? "▼" : "▶"} Download
      </button>
      {open && (
        <div className="panel-body download-body">
          {isProtein && r.merged_dna && (
            <button
              className="dl-btn"
              onClick={() => downloadBlob(r.merged_dna!, `${r.sample}_contig.fasta`)}
            >
              Download DNA contig
            </button>
          )}
          {isProtein && r.longest_orf && (
            <button
              className="dl-btn"
              onClick={() => downloadBlob(r.longest_orf!, `${r.sample}_ORF.fasta`)}
            >
              Download Protein ORF
            </button>
          )}
          {!isProtein && r.merged_dna && (
            <button
              className="dl-btn"
              onClick={() => downloadBlob(r.merged_dna!, `${r.sample}_merged.fasta`)}
            >
              Download merged DNA contig
            </button>
          )}
          {!isProtein && r.cropped_dna && (
            <button
              className="dl-btn"
              onClick={() => {
                const label =
                  sequenceType === "paired"
                    ? `${r.sample}_cropped.fasta`
                    : `${r.sample}_sequence.fasta`;
                downloadBlob(r.cropped_dna!, label);
              }}
            >
              {sequenceType === "paired" ? "Download cropped DNA sequence" : "Download sequence"}
            </button>
          )}
          {!r.merged_dna && !r.cropped_dna && !r.longest_orf && (
            <p className="muted">No downloadable sequences available.</p>
          )}
        </div>
      )}
    </div>
  );
}
