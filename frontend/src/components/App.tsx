import React, { useState, useRef, useEffect } from "react";
import Lottie from "lottie-react";
import { useLoadingText } from "../hooks/useLoadingText";
import { submitProcess } from "../api/client";
import Header from "./Header";
import SequenceTypeSelector from "./SequenceTypeSelector";
import PrimerTypeSelector from "./PrimerTypeSelector";
import ReferenceTypeCheckbox from "./ReferenceTypeCheckbox";
import FileUploader from "./FileUploader";
import MinLcsSlider from "./MinLcsSlider";
import FeatureCategoryPicker from "./FeatureCategoryPicker";
import ProcessButton from "./ProcessButton";
import ResultsTabs from "./ResultsTabs";
import AlignTab from "./AlignTab";

// Lottie JSON assets
import geckoAnimation from "../assets/lilthang.json";

// SVG for gecko footprint (gecko facing/moving left, toes pointing forward-left)
const FOOTPRINT_SVG = `<svg width="30" height="18" viewBox="0 0 30 18" xmlns="http://www.w3.org/2000/svg">
  <ellipse cx="23" cy="9"  rx="5"   ry="4.5" fill="rgba(62,114,29,0.55)"/>
  <ellipse cx="14" cy="9"  rx="3.5" ry="3"   fill="rgba(62,114,29,0.40)"/>
  <ellipse cx="4"  cy="2"  rx="2.2" ry="4"   fill="rgba(62,114,29,0.55)" transform="rotate(28 4 2)"/>
  <ellipse cx="3"  cy="8"  rx="2.2" ry="4"   fill="rgba(62,114,29,0.55)" transform="rotate(8 3 8)"/>
  <ellipse cx="5"  cy="14" rx="2.2" ry="4"   fill="rgba(62,114,29,0.55)" transform="rotate(-12 5 14)"/>
  <ellipse cx="10" cy="17" rx="2.2" ry="3.5" fill="rgba(62,114,29,0.55)" transform="rotate(-32 10 17)"/>
</svg>`;

export default function App() {
  const [activeView, setActiveView] = useState<"pipeline" | "align">("pipeline");
  const [sequenceType, setSequenceType] = useState("paired");
  const [primerType, setPrimerType] = useState("protein");
  const [useAaRefs, setUseAaRefs] = useState(false);
  const [seqFiles, setSeqFiles] = useState<File[]>([]);
  const [refFile, setRefFile] = useState<File | null>(null);
  const [minLcs, setMinLcs] = useState(12);
  const [selectedCategories, setSelectedCategories] = useState<string[] | null>(null);

  const [isProcessing, setIsProcessing] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const loadingText = useLoadingText(isProcessing);

  // Refs for gecko footprint animation
  const overlayRef      = useRef<HTMLDivElement>(null);
  const geckoWalkerRef  = useRef<HTMLDivElement>(null);
  const footCountRef    = useRef(0);

  useEffect(() => {
    if (!isProcessing) {
      footCountRef.current = 0;
      return;
    }
    const id = setInterval(() => {
      if (!geckoWalkerRef.current || !overlayRef.current) return;
      const geckoRect   = geckoWalkerRef.current.getBoundingClientRect();
      const overlayRect = overlayRef.current.getBoundingClientRect();
      if (geckoRect.width === 0) return; // not yet painted

      const count      = footCountRef.current++;
      const isLeftFoot = count % 2 === 0;

      // Tune these two constants if footprints drift from the gecko:
      //   FOOT_X – horizontal fraction across the Lottie frame (0 = left edge, 1 = right edge)
      //   FOOT_Y – vertical fraction down the Lottie frame  (0 = top,  1 = bottom)
      const FOOT_X = 0.50;
      const FOOT_Y = 0.50;

      const x = geckoRect.left - overlayRect.left + geckoRect.width  * FOOT_X;
      const y = geckoRect.top  - overlayRect.top  + geckoRect.height * FOOT_Y;
      // Alternate slightly above/below centre line for left/right feet
      const yOffset = isLeftFoot ? -18 : 18;

const fp = document.createElement("div");
      fp.className = "gecko-footprint";
      fp.style.left = `${x - 15}px`;
      fp.style.top  = `${y + yOffset - 9}px`;
      // Mirror vertically for the right foot
      if (!isLeftFoot) fp.style.transform = "scaleY(-1)";
      fp.innerHTML = FOOTPRINT_SVG;

      overlayRef.current.appendChild(fp);
      setTimeout(() => fp.remove(), 2200);
    }, 350);

    return () => clearInterval(id);
  }, [isProcessing]);

  const [results, setResults] = useState<unknown[] | null>(null);
  const [loadingProgress, setLoadingProgress] = useState<string | null>(null);

  async function handleProcess() {
    setError(null);
    if (seqFiles.length === 0) {
      setError("Gecko has nothing to sniff! Drop at least one sequence file first.");
      return;
    }
    if (!refFile) {
      setError("Gecko needs the reference map too — drop your .xlsx file!");
      return;
    }

    const formData = new FormData();
    seqFiles.forEach((f) => formData.append("seq_files", f));
    formData.append("ref_file", refFile);
    formData.append("sequence_type", sequenceType);
    formData.append("primer_type", primerType);
    formData.append("use_aa_refs", useAaRefs ? "true" : "false");
    formData.append("min_lcs", minLcs.toString());
    if (selectedCategories) {
      formData.append("selected_categories", JSON.stringify(selectedCategories));
    }

    setIsProcessing(true);
    try {
      const res = await submitProcess(formData, (processed, total) => {
        setLoadingProgress(`Processing ${processed} / ${total} files...`);
      });
      setResults(res);
    } catch (e) {
      setError(e instanceof Error ? e.message : "Unknown error");
    } finally {
      setIsProcessing(false);
      setLoadingProgress(null);
    }
  }

  return (
    <div className={`app-container${results ? " has-results" : ""}`}>
      <Header />

      <div className="main-layout">
        <nav className="view-tabs">
          <button
            className={`view-tab ${activeView === "pipeline" ? "active" : ""}`}
            onClick={() => setActiveView("pipeline")}
          >
            Multiple Sequence Checker
          </button>
          <button
            className={`view-tab ${activeView === "align" ? "active" : ""}`}
            onClick={() => setActiveView("align")}
          >
            Sequence Alignment
          </button>
        </nav>

        <div className="main-content">
      {activeView === "align" ? (
        <AlignTab />
      ) : (
        <>
          {isProcessing && (
            <div ref={overlayRef} className="gecko-loading-overlay">
              <p className="gecko-loading-text">{loadingText}</p>
              {loadingProgress && <p className="gecko-loading-progress">{loadingProgress}</p>}
              <div ref={geckoWalkerRef} className="gecko-walker">
                <Lottie animationData={geckoAnimation} loop={true} style={{ width: 180, height: 180 }} />
              </div>
            </div>
          )}

          {!results ? (
            <div className="form-section">
              <SequenceTypeSelector value={sequenceType} onChange={setSequenceType} />
              {sequenceType !== "circular" && (
                <PrimerTypeSelector value={primerType} onChange={setPrimerType} />
              )}
              {primerType === "protein" && sequenceType !== "circular" && (
                <ReferenceTypeCheckbox value={useAaRefs} onChange={setUseAaRefs} />
              )}

              <FileUploader
                sequenceType={sequenceType}
                primerType={primerType}
                useAaRefs={useAaRefs}
                onSeqFilesChange={setSeqFiles}
                onRefFileChange={setRefFile}
                seqFiles={seqFiles}
                refFile={refFile}
              />

              {sequenceType === "paired" && (
                <MinLcsSlider value={minLcs} onChange={setMinLcs} />
              )}

              <FeatureCategoryPicker
                selected={selectedCategories}
                onChange={setSelectedCategories}
              />

              {error && <div className="error-banner">{error}</div>}

              <ProcessButton onClick={handleProcess} isProcessing={isProcessing} loadingText={loadingText} />
            </div>
          ) : (
            <div className="results-section">
              <button
                className="back-btn"
                onClick={() => window.location.reload()}
                title="Change Submission"
              >
                ← New submission
              </button>
              {results.length === 0 ? (
                <div className="warning-banner">Gecko came up empty-handed — no matching sequences found in those files.</div>
              ) : (
                <ResultsTabs
                  results={results}
                  sequenceType={sequenceType}
                  primerType={primerType}
                />
              )}
            </div>
          )}
        </>
      )}
        </div>
      </div>
    </div>
  );
}
