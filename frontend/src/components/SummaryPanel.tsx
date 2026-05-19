import React, { useState } from "react";
import { useCountUp } from "../hooks/useCountUp";

interface MergeInfo {
  lcs_length?: number;
}

interface CropInfo {
  start?: number;
  end?: number;
}

interface Result {
  merged_dna?: string;
  cropped_dna?: string;
  identity_percent?: number;
  longest_orf?: string;
  orientation?: string;
  frame?: string;
  merge_info?: MergeInfo;
  crop_info?: CropInfo;
}

interface Props {
  result: Result;
  sequenceType: string;
  primerType: string;
}

export default function SummaryPanel({ result, sequenceType, primerType }: Props) {
  const [open, setOpen] = useState(true);

  const isProtein = primerType === "protein";
  const r = result;

  const dnaLength    = (r.merged_dna ?? r.cropped_dna ?? "").length;
  const mergedLength = (r.merged_dna  ?? "").length;
  const croppedLength= (r.cropped_dna ?? "").length;
  const orfLength    = (r.longest_orf ?? "").length;
  const identityRaw  = r.identity_percent ?? 0;
  const orientation  = r.orientation ?? "forward";
  const frame        = r.frame ?? "N/A";

  // Count-up animations — always called unconditionally (rules of hooks)
  const animDna     = useCountUp(dnaLength,     600, 0, open);
  const animMerged  = useCountUp(mergedLength,  600, 0, open);
  const animCropped = useCountUp(croppedLength, 600, 0, open);
  const animOrf     = useCountUp(orfLength,     600, 0, open);
  const animIdentity= useCountUp(identityRaw,   600, 1, open);

  return (
    <div className="panel">
      <button className="panel-header" onClick={() => setOpen(!open)}>
        {open ? "▼" : "▶"} Summary
      </button>
      {open && (
        <div className="panel-body">
          <div className="metrics-row">
            {isProtein && sequenceType === "paired" && (
              <>
                <div className="metric-card">
                  <span className="metric-value">{animDna}</span>
                  <span className="metric-label">Assembled DNA length</span>
                </div>
                {orfLength > 70 && (
                  <div className="metric-card">
                    <span className="metric-value">{animOrf}</span>
                    <span className="metric-label">Protein ORF length</span>
                  </div>
                )}
              </>
            )}
            {isProtein && sequenceType !== "paired" && (
              <>
                <div className="metric-card">
                  <span className="metric-value">{animDna}</span>
                  <span className="metric-label">DNA length</span>
                </div>
                {orfLength > 70 && (
                  <div className="metric-card">
                    <span className="metric-value">{animOrf}</span>
                    <span className="metric-label">Protein ORF length</span>
                  </div>
                )}
              </>
            )}
            {!isProtein && sequenceType === "paired" && (
              <>
                <div className="metric-card">
                  <span className="metric-value">{animMerged}</span>
                  <span className="metric-label">Merged DNA length</span>
                </div>
                <div className="metric-card">
                  <span className="metric-value">{animCropped}</span>
                  <span className="metric-label">Cropped DNA length</span>
                </div>
              </>
            )}
            {!isProtein && sequenceType !== "paired" && (
              <div className="metric-card">
                <span className="metric-value">{animCropped}</span>
                <span className="metric-label">Sequence length</span>
              </div>
            )}
            <div className="metric-card">
              <span className="metric-value">{animIdentity}%</span>
              <span className="metric-label">Percent Identity</span>
            </div>
          </div>

          {sequenceType === "circular" && (
            <div className="info-badge">
              Circular plasmid sequence
              {orientation === "reverse complement" && " (reverse complement)"}
              {isProtein && ` | Frame: ${frame}`}
            </div>
          )}
          {sequenceType === "unidirectional" && (
            <div className="info-badge">
              Unidirectional sequence
              {orientation === "reverse complement"
                ? " (reverse complement)"
                : " (local alignment)"}
              {isProtein && ` | Frame: ${frame}`}
            </div>
          )}
          {!isProtein && sequenceType === "paired" && (
            <>
              {r.merge_info?.lcs_length && (
                <div className="info-badge">LCS overlap: {r.merge_info.lcs_length} bp</div>
              )}
              {r.crop_info?.start && r.crop_info?.end && (
                <div className="info-badge">
                  Cropped region: positions {r.crop_info.start}-{r.crop_info.end}
                </div>
              )}
            </>
          )}
        </div>
      )}
    </div>
  );
}
