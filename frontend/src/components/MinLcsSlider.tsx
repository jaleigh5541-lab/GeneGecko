import React from "react";

interface Props {
  value: number;
  onChange: (value: number) => void;
}

export default function MinLcsSlider({ value, onChange }: Props) {
  return (
    <div className="form-group">
      <label>
        Minimum overlap: <strong>{value} bp</strong>
      </label>
      <input
        type="range"
        min={6}
        max={50}
        value={value}
        onChange={(e) => onChange(Number(e.target.value))}
        className="slider"
      />
      <small className="form-hint">Gecko default is 12 bp.</small>
    </div>
  );
}
