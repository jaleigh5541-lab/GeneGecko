import React, { useState } from "react";

interface Result {
  orf_count?: number;
  alignment_text?: string;
  alignment_text_with_orfs?: string;
}

interface Props {
  result: Result;
}

export default function AlignmentPanel({ result }: Props) {
  const [open, setOpen] = useState(false);
  const [showOrfs, setShowOrfs] = useState(false);

  const hasOrfs = (result.orf_count ?? 0) > 0;
  const html =
    showOrfs && result.alignment_text_with_orfs
      ? result.alignment_text_with_orfs
      : result.alignment_text;

  return (
    <div className="panel">
      <button className="panel-header" onClick={() => setOpen(!open)}>
        {open ? "▼" : "▶"} Alignment
        {hasOrfs && (
          <label
            onClick={(e) => e.stopPropagation()}
            style={{ marginLeft: "1rem", fontWeight: "normal", fontSize: "0.85rem" }}
          >
            <input
              type="checkbox"
              checked={showOrfs}
              onChange={(e) => setShowOrfs(e.target.checked)}
              style={{ marginRight: "0.3rem" }}
            />
            Show ORFs ({result.orf_count} detected)
          </label>
        )}
      </button>
      {open && (
        <div className="panel-body alignment-body">
          {html ? (
            <div dangerouslySetInnerHTML={{ __html: html }} />
          ) : (
            <p className="muted">No alignment available.</p>
          )}
        </div>
      )}
    </div>
  );
}
