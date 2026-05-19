import { useState, useEffect } from "react";

/**
 * Animates a number from 0 to `target` over `duration` ms with an ease-out
 * cubic curve. Re-triggers whenever `target` or `trigger` changes.
 * When `trigger` is false the returned value resets to 0 immediately.
 */
export function useCountUp(
  target: number,
  duration = 600,
  decimals = 0,
  trigger = true,
): string {
  const [value, setValue] = useState(0);

  useEffect(() => {
    if (!trigger) {
      setValue(0);
      return;
    }

    let startTime: number | null = null;
    let raf: number;

    const animate = (timestamp: number) => {
      if (startTime === null) startTime = timestamp;
      const progress = Math.min((timestamp - startTime) / duration, 1);
      // ease-out cubic
      const eased = 1 - Math.pow(1 - progress, 3);
      setValue(eased * target);
      if (progress < 1) raf = requestAnimationFrame(animate);
    };

    raf = requestAnimationFrame(animate);
    return () => cancelAnimationFrame(raf);
  }, [target, duration, trigger]);

  return value.toFixed(decimals);
}
