import React from "react";

const OPTIONS = [
  { value: "protein", label: "T7/T7-term" },
  { value: "intron", label: "pUC19-seqback/M13R" },
];

interface Props {
  value: string;
  onChange: (value: string) => void;
}

export default function PrimerTypeSelector({ value, onChange }: Props) {
  return (
    <fieldset className="form-group">
      <legend>Primer(s)</legend>
      <div className="radio-row">
        {OPTIONS.map((opt) => (
          <label key={opt.value} className="radio-label">
            <input
              type="radio"
              name="primer_type"
              value={opt.value}
              checked={value === opt.value}
              onChange={() => onChange(opt.value)}
            />
            {opt.label}
          </label>
        ))}
      </div>
    </fieldset>
  );
}
