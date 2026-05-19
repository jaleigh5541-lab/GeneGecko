import React from "react";
import SummaryPanel from "./SummaryPanel";
import AlignmentPanel from "./AlignmentPanel";
import GeckoInvestigate from "./GeckoInvestigate";
import FeaturesTable from "./FeaturesTable";
import DownloadPanel from "./DownloadPanel";

interface Feature {
  name: string;
  category: string;
  color: string;
  start: number;
  end: number;
  strand: string;
}

interface Result {
  sample?: string;
  error?: string;
  investigate_html?: string;
  detected_features?: Feature[];
  [key: string]: unknown;
}

interface Props {
  result: Result;
  sequenceType: string;
  primerType: string;
}

export default function SampleTab({ result, sequenceType, primerType }: Props) {
  return (
    <div className="sample-tab">
      <h2 className="sample-title">{result.sample}</h2>
      {result.error && (
        <div className="error-banner">Processing error: {result.error}</div>
      )}
      <SummaryPanel result={result} sequenceType={sequenceType} primerType={primerType} />
      <AlignmentPanel result={result} />
      {primerType === "protein" && sequenceType === "paired" && (
        <GeckoInvestigate investigateHtml={result.investigate_html} />
      )}
      <FeaturesTable features={result.detected_features ?? []} />
      <DownloadPanel result={result} sequenceType={sequenceType} primerType={primerType} />
    </div>
  );
}
