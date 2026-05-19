import { useState, useEffect, useRef } from "react";

const PHRASES = [
  "Sticking to things...",
  "Eating flies...",
  "Basking...",
  "Licking eyeballs...",
  "Defying gravity...",
  "Scurrying...",
  "Using the force...",
  "Hiding under rocks...",
  "Distracting predators...",
  "Wiggling...",
  "Tree-hopping...",
  "Studying geometry...",
];

const CHAR_INTERVAL_MS = 38;  // ms between each typed character
const HOLD_MS          = 1600; // ms to hold the completed message before clearing

function shuffle(arr: string[]): string[] {
  const a = [...arr];
  for (let i = a.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [a[i], a[j]] = [a[j]!, a[i]!];
  }
  return a;
}

export function useLoadingText(isActive: boolean): string {
  const [displayed, setDisplayed] = useState("");
  const queueRef   = useRef<string[]>([]);
  const currentRef = useRef("");
  const charIdxRef = useRef(0);
  // "typing" → "holding" → move to next phrase
  const phaseRef   = useRef<"typing" | "holding">("typing");

  useEffect(() => {
    if (!isActive) {
      setDisplayed("");
      return;
    }

    function pickNext() {
      if (queueRef.current.length === 0) queueRef.current = shuffle(PHRASES);
      currentRef.current = queueRef.current.pop()!;
      charIdxRef.current = 0;
      phaseRef.current   = "typing";
    }

    pickNext();

    let timeout: ReturnType<typeof setTimeout>;

    function tick() {
      if (phaseRef.current === "typing") {
        charIdxRef.current++;
        const slice = currentRef.current.slice(0, charIdxRef.current);
        setDisplayed(slice);

        if (charIdxRef.current >= currentRef.current.length) {
          // Finished typing — hold, then clear and start next phrase
          phaseRef.current = "holding";
          timeout = setTimeout(() => {
            setDisplayed("");
            timeout = setTimeout(() => {
              pickNext();
              tick();
            }, 180); // brief blank gap between phrases
          }, HOLD_MS);
        } else {
          timeout = setTimeout(tick, CHAR_INTERVAL_MS);
        }
      }
    }

    timeout = setTimeout(tick, CHAR_INTERVAL_MS);
    return () => clearTimeout(timeout);
  }, [isActive]);

  return displayed;
}
