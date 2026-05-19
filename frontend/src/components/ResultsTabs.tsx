import React, { useState } from "react";
import SampleTab from "./SampleTab";
import SummaryTable from "./SummaryTable";

interface Result {
  sample?: string;
  [key: string]: unknown;
}

interface Props {
  results: unknown[];
  sequenceType: string;
  primerType: string;
}

export default function ResultsTabs({ results, sequenceType, primerType }: Props) {
  const [activeIdx, setActiveIdx] = useState(0);
  const typedResults = results as Result[];

  return (
    <div className="results-tabs">
      <SummaryTable
        results={typedResults}
        sequenceType={sequenceType}
        primerType={primerType}
        onSelectSample={setActiveIdx}
      />
      <div className="tab-bar">
        {typedResults.map((r, i) => (
          <button
            key={i}
            className={`tab-btn ${i === activeIdx ? "active" : ""}`}
            onClick={() => setActiveIdx(i)}
          >
            {r.sample}
          </button>
        ))}
      </div>
      <div className="tab-content">
        <SampleTab
          result={typedResults[activeIdx]!}
          sequenceType={sequenceType}
          primerType={primerType}
        />
      </div>
    </div>
  );
}
