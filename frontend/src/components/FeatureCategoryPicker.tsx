import React, { useState, useEffect } from "react";
import { fetchCategories } from "../api/client";

interface Props {
  selected: string[] | null;
  onChange: (value: string[] | null) => void;
}

export default function FeatureCategoryPicker({ selected, onChange }: Props) {
  const [categories, setCategories] = useState<string[]>([]);
  const [colors, setColors] = useState<Record<string, string>>({});
  const [open, setOpen] = useState(false);

  useEffect(() => {
    fetchCategories()
      .then((data) => {
        setCategories(data.categories);
        setColors(data.colors);
      })
      .catch(() => {});
  }, []);

  const allSelected = selected === null;

  function toggle(cat: string) {
    let next: string[] | null;
    if (allSelected) {
      next = categories.filter((c) => c !== cat);
    } else if (selected!.includes(cat)) {
      next = selected!.filter((c) => c !== cat);
    } else {
      next = [...selected!, cat];
    }
    if (next.length === categories.length || next.length === 0) next = null;
    onChange(next);
  }

  function isChecked(cat: string): boolean {
    return allSelected || (selected !== null && selected.includes(cat));
  }

  if (categories.length === 0) return null;

  return (
    <div className="form-group">
      <button
        type="button"
        className="expander-toggle"
        onClick={() => setOpen(!open)}
      >
        {open ? "▼" : "▶"} Feature Detection Settings
      </button>
      {open && (
        <div className="category-grid">
          {categories.map((cat) => (
            <label
              key={cat}
              className="checkbox-label category-chip"
              style={{ borderColor: colors[cat] ?? "#6b7280" }}
            >
              <input
                type="checkbox"
                checked={isChecked(cat)}
                onChange={() => toggle(cat)}
              />
              <span
                className="cat-dot"
                style={{ backgroundColor: colors[cat] ?? "#6b7280" }}
              />
              {cat}
            </label>
          ))}
        </div>
      )}
    </div>
  );
}
