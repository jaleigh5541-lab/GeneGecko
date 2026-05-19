import React, { useState, useEffect } from "react";
import geckoIcon from "../assets/gecko-investigate.svg";

interface Props {
  investigateHtml?: string;
}

export default function GeckoInvestigate({ investigateHtml }: Props) {
  const [open, setOpen] = useState(false);
  const [showOrfs, setShowOrfs] = useState(false);

  if (!investigateHtml) return null;

  // Close on Escape
  // eslint-disable-next-line react-hooks/rules-of-hooks
  useEffect(() => {
    if (!open) return;
    function onKey(e: KeyboardEvent) {
      if (e.key === "Escape") setOpen(false);
    }
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [open]);

  // Prevent background scroll while modal is open
  // eslint-disable-next-line react-hooks/rules-of-hooks
  useEffect(() => {
    document.body.style.overflow = open ? "hidden" : "";
    return () => { document.body.style.overflow = ""; };
  }, [open]);

  return (
    <>
      <button
        className="gi-fab"
        onClick={() => setOpen(true)}
        title="Gecko Investigate"
        aria-label="Open Gecko Investigate"
      >
        <img src={geckoIcon.src ?? geckoIcon} alt="Gecko Investigate" className="gi-fab-icon" />
      </button>

      {open && (
        <div className="gi-overlay" onClick={() => setOpen(false)}>
          <div className="gi-modal" onClick={(e) => e.stopPropagation()}>
            <div className="gi-modal-header">
              <img src={geckoIcon.src ?? geckoIcon} alt="" className="gi-modal-icon" />
              <span className="gi-modal-title">Gecko Investigate</span>
              <label className="gi-orf-toggle">
                <input
                  type="checkbox"
                  checked={showOrfs}
                  onChange={(e) => setShowOrfs(e.target.checked)}
                />
                Show ORFs
              </label>
              <button
                className="gi-modal-close"
                onClick={() => setOpen(false)}
                aria-label="Close"
              >
                ✕
              </button>
            </div>
            <div className="gi-modal-body">
              <div
                className={showOrfs ? undefined : "gi-hide-orfs"}
                dangerouslySetInnerHTML={{ __html: investigateHtml }}
              />
            </div>
          </div>
        </div>
      )}
    </>
  );
}
