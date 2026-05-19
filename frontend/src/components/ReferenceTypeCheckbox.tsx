import React from "react";

interface Props {
  value: boolean;
  onChange: (value: boolean) => void;
}

export default function ReferenceTypeCheckbox({ value, onChange }: Props) {
  return (
    <div className="form-group">
      <label className="checkbox-label">
        <input
          type="checkbox"
          checked={value}
          onChange={(e) => onChange(e.target.checked)}
        />
        References are amino acid sequences
      </label>
      <small className="form-hint">
        If unchecked (default), aligns query DNA to DNA reference first, then
        translates. If checked, translates to amino acids first, then aligns to
        AA reference.
      </small>
    </div>
  );
}
