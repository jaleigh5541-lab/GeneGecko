import React, { useState, useCallback } from "react";
import Lottie from "lottie-react";
import lizardAnimation from "../assets/toungy-lizard.json";
import GECKO_FACTS from "../assets/geckoFacts";

export default function Header() {
  const [fact, setFact] = useState<string | null>(null);

  const showFact = useCallback(() => {
    setFact(GECKO_FACTS[Math.floor(Math.random() * GECKO_FACTS.length)]);
  }, []);

  const hideFact = useCallback(() => {
    setFact(null);
  }, []);

  return (
    <header className="header">
      <div className="header-title-group">
        <div
          className="gecko-fact-anchor"
          onMouseEnter={showFact}
          onMouseLeave={hideFact}
        >
          <Lottie animationData={lizardAnimation} loop={true} style={{ width: 64, height: 64 }} />
          {fact && (
            <div className="speech-bubble">
              {fact}
            </div>
          )}
        </div>
        <div>
          <h1>Gene-Gecko</h1>
          <p className="subtitle">DNA &amp; Protein Sequence Analysis · Powered by gecko instincts</p>
        </div>
      </div>
    </header>
  );
}
