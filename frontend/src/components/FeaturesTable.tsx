import React, { useState } from "react";

interface Feature {
  name: string;
  category: string;
  color: string;
  start: number;
  end: number;
  strand: string;
}

interface Props {
  features: Feature[];
}

export default function FeaturesTable({ features }: Props) {
  const [open, setOpen] = useState(false);

  return (
    <div className="panel">
      <button className="panel-header" onClick={() => setOpen(!open)}>
        {open ? "▼" : "▶"} Detected features ({features.length})
      </button>
      {open && (
        <div className="panel-body">
          {features.length === 0 ? (
            <p className="muted">No features detected in this sequence.</p>
          ) : (
            <>
              <table className="features-table">
                <thead>
                  <tr>
                    <th>Feature</th>
                    <th>Category</th>
                    <th>Position</th>
                    <th>Strand</th>
                    <th>Length</th>
                  </tr>
                </thead>
                <tbody>
                  {features.map((f, i) => (
                    <tr key={i}>
                      <td>{f.name}</td>
                      <td>
                        <span
                          className="cat-badge"
                          style={{ backgroundColor: `${f.color}40`, color: f.color }}
                        >
                          {f.category}
                        </span>
                      </td>
                      <td>{f.start}-{f.end}</td>
                      <td>{f.strand}</td>
                      <td>{f.end - f.start + 1}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
              <div className="category-summary">
                {Object.entries(
                  features.reduce<Record<string, { count: number; color: string }>>(
                    (acc, f) => {
                      if (!acc[f.category]) acc[f.category] = { count: 0, color: f.color };
                      acc[f.category]!.count++;
                      return acc;
                    },
                    {}
                  )
                ).map(([cat, { count, color }]) => (
                  <span
                    key={cat}
                    className="cat-badge"
                    style={{ backgroundColor: `${color}40`, color }}
                  >
                    {cat}: {count}
                  </span>
                ))}
              </div>
            </>
          )}
        </div>
      )}
    </div>
  );
}
