import React from "react";

const OPTIONS = [
  { value: "paired", label: "Paired sequences" },
  { value: "unidirectional", label: "Unidirectional sequences" },
  { value: "circular", label: "Circular plasmid sequences" },
];

interface Props {
  value: string;
  onChange: (value: string) => void;
}

export default function SequenceTypeSelector({ value, onChange }: Props) {
  return (
    <fieldset className="form-group">
      <legend>Select sequence type</legend>
      <div className="radio-row">
        {OPTIONS.map((opt) => (
          <label key={opt.value} className="radio-label">
            <input
              type="radio"
              name="sequence_type"
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
