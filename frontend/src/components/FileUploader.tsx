import React, { useRef, useState } from "react";

function acceptFor(sequenceType: string): string {
  if (sequenceType === "circular") return ".fasta,.fa";
  if (sequenceType === "unidirectional") return ".seq,.fasta,.fa";
  return ".seq";
}

/** Recursively read files via File System Access API (async iterator). */
async function readViaHandle(dir: FileSystemDirectoryHandle): Promise<File[]> {
  const out: File[] = [];
  for await (const handle of dir.values()) {
    if (handle.kind === "file") {
      out.push(await (handle as FileSystemFileHandle).getFile());
    } else {
      out.push(...await readViaHandle(handle as FileSystemDirectoryHandle));
    }
  }
  return out;
}

/** Fallback: read files via legacy readEntries API. */
async function readViaEntries(entry: FileSystemDirectoryEntry): Promise<File[]> {
  const reader = entry.createReader();
  const out: File[] = [];
  let batch: FileSystemEntry[];
  do {
    batch = await new Promise<FileSystemEntry[]>((resolve, reject) =>
      reader.readEntries(resolve, reject),
    );
    for (const item of batch) {
      if (item.isFile) {
        const file = await new Promise<File>((resolve, reject) =>
          (item as FileSystemFileEntry).file(resolve, reject),
        );
        out.push(file);
      } else if (item.isDirectory) {
        out.push(...await readViaEntries(item as FileSystemDirectoryEntry));
      }
    }
  } while (batch.length > 0);
  return out;
}

/** Collect File objects from a drop event, handling both files and folders. */
async function filesFromDrop(dt: DataTransfer): Promise<File[]> {
  const items = dt.items;
  if (!items?.length) return Array.from(dt.files);

  // Synchronously initiate handle requests AND grab legacy entries before any await
  type ItemWithHandle = DataTransferItem & {
    getAsFileSystemHandle?: () => Promise<FileSystemHandle>;
  };
  const handlePs: Promise<FileSystemHandle | null>[] = [];
  const legacyEntries: (FileSystemEntry | null)[] = [];
  for (let i = 0; i < items.length; i++) {
    legacyEntries.push(items[i].webkitGetAsEntry?.() ?? null);
    const item = items[i] as ItemWithHandle;
    if (typeof item.getAsFileSystemHandle === "function") {
      handlePs.push(item.getAsFileSystemHandle().catch(() => null));
    }
  }

  // ── Primary: File System Access API (most reliable folder enumeration) ──
  if (handlePs.length) {
    try {
      const handles = await Promise.all(handlePs);
      const collected: File[] = [];
      for (const h of handles) {
        if (!h) continue;
        if (h.kind === "directory") {
          collected.push(...await readViaHandle(h as FileSystemDirectoryHandle));
        } else if (h.kind === "file") {
          collected.push(await (h as FileSystemFileHandle).getFile());
        }
      }
      if (collected.length) return collected;
    } catch {
      /* fall through */
    }
  }

  // ── Fallback: webkitGetAsEntry + readEntries ──
  const collected: File[] = [];
  let hasDir = false;
  for (const entry of legacyEntries) {
    if (entry?.isDirectory) {
      hasDir = true;
      collected.push(...await readViaEntries(entry as FileSystemDirectoryEntry));
    } else if (entry?.isFile) {
      const file = await new Promise<File>((resolve, reject) =>
        (entry as FileSystemFileEntry).file(resolve, reject),
      );
      collected.push(file);
    }
  }
  return hasDir || collected.length ? collected : Array.from(dt.files);
}

function labelFor(sequenceType: string): string {
  if (sequenceType === "circular") return "Sequence files (.fasta)";
  if (sequenceType === "unidirectional") return "Sequence files (.seq or .fasta)";
  return "Paired reads (.seq)";
}

interface Props {
  sequenceType: string;
  primerType: string;
  useAaRefs: boolean;
  onSeqFilesChange: (files: File[]) => void;
  onRefFileChange: (file: File | null) => void;
  seqFiles: File[];
  refFile: File | null;
}

export default function FileUploader({
  sequenceType,
  onSeqFilesChange,
  onRefFileChange,
  seqFiles,
  refFile,
}: Props) {
  const seqRef = useRef<HTMLInputElement>(null);
  const seqFolderRef = useRef<HTMLInputElement>(null);
  const refRef = useRef<HTMLInputElement>(null);
  const [showHelp, setShowHelp] = useState(false);

  const [seqDrag, setSeqDrag] = useState(false);
  const [refDrag, setRefDrag] = useState(false);
  const [seqTick, setSeqTick] = useState(false);
  const [refTick, setRefTick] = useState(false);

  function flashTick(which: "seq" | "ref") {
    if (which === "seq") {
      setSeqTick(true);
      setTimeout(() => setSeqTick(false), 900);
    } else {
      setRefTick(true);
      setTimeout(() => setRefTick(false), 900);
    }
  }

  function handleSeqChange(e: React.ChangeEvent<HTMLInputElement>) {
    const files = Array.from(e.target.files ?? []);
    onSeqFilesChange(files);
    if (files.length) flashTick("seq");
  }

  function handleFolderChange(e: React.ChangeEvent<HTMLInputElement>) {
    const accepted = acceptFor(sequenceType).split(",");
    const files = Array.from(e.target.files ?? []).filter((f) =>
      accepted.some((ext) => f.name.toLowerCase().endsWith(ext)),
    );
    if (files.length) {
      onSeqFilesChange(files);
      flashTick("seq");
    }
  }

  function handleRefChange(e: React.ChangeEvent<HTMLInputElement>) {
    const file = e.target.files?.[0] ?? null;
    onRefFileChange(file);
    if (file) flashTick("ref");
  }

  async function handleSeqDrop(e: React.DragEvent<HTMLDivElement>) {
    e.preventDefault();
    setSeqDrag(false);
    const accepted = acceptFor(sequenceType).split(",");
    const allFiles = await filesFromDrop(e.dataTransfer);
    const files = allFiles.filter((f) =>
      accepted.some((ext) => f.name.toLowerCase().endsWith(ext)),
    );
    if (files.length) {
      onSeqFilesChange(files);
      flashTick("seq");
    }
  }

  async function handleRefDrop(e: React.DragEvent<HTMLDivElement>) {
    e.preventDefault();
    setRefDrag(false);
    const allFiles = await filesFromDrop(e.dataTransfer);
    const file =
      allFiles.find((f) => {
        const name = f.name.toLowerCase();
        return name.endsWith(".xlsx") || name.endsWith(".csv");
      }) ?? null;
    if (file) {
      onRefFileChange(file);
      flashTick("ref");
    }
  }

  const seqZoneClass = [
    "drop-zone",
    seqDrag ? "drop-zone--over" : "",
    seqTick ? "drop-zone--tick" : "",
  ]
    .filter(Boolean)
    .join(" ");

  const refZoneClass = [
    "drop-zone",
    refDrag ? "drop-zone--over" : "",
    refTick ? "drop-zone--tick" : "",
  ]
    .filter(Boolean)
    .join(" ");

  return (
    <div className="form-group">
      {/* ── Sequence files ── */}
      <label className="file-label">{labelFor(sequenceType)}</label>
      <div
        className={seqZoneClass}
        onDragOver={(e) => { e.preventDefault(); setSeqDrag(true); }}
        onDragLeave={() => setSeqDrag(false)}
        onDrop={handleSeqDrop}
        onClick={() => seqRef.current?.click()}
        role="button"
        tabIndex={0}
        onKeyDown={(e) => e.key === "Enter" && seqRef.current?.click()}
        aria-label="Drop sequence files or click to browse"
      >
        <input
          ref={seqRef}
          type="file"
          accept={acceptFor(sequenceType)}
          multiple
          onChange={handleSeqChange}
          style={{ display: "none" }}
        />
        {/* Native folder picker — completely separate browser code path */}
        <input
          ref={seqFolderRef}
          type="file"
          // @ts-expect-error webkitdirectory is non-standard but widely supported
          webkitdirectory=""
          onChange={handleFolderChange}
          style={{ display: "none" }}
        />
        <span className="drop-zone__icon" aria-hidden="true">
          {seqTick ? "✓" : "↓"}
        </span>
        <span className="drop-zone__text">
          {seqTick
            ? "Files received!"
            : seqDrag
            ? "Release to drop"
            : "Drag files here or click to browse"}
        </span>
        <span className="drop-zone__hint">
          {acceptFor(sequenceType)}
          {" | "}
          <button
            type="button"
            className="folder-pick-link"
            onClick={(e) => { e.stopPropagation(); seqFolderRef.current?.click(); }}
          >
            select folder
          </button>
        </span>
      </div>
      {seqFiles.length > 0 && (
        <div className="file-list">
          {seqFiles.map((f, i) => (
            <span key={i} className="file-chip">
              {f.name}
            </span>
          ))}
        </div>
      )}

      {/* ── Reference map ── */}
      <div style={{ marginTop: "0.75rem" }}>
        <span className="file-label" style={{ display: "inline" }}>
          Reference map (.xlsx or .csv)
        </span>
        <button
          type="button"
          className="ref-help-toggle"
          onClick={() => setShowHelp((v) => !v)}
          aria-expanded={showHelp}
          title="Show column format help"
        >
          ⓘ
        </button>
      </div>
      <div
        className={refZoneClass}
        onDragOver={(e) => { e.preventDefault(); setRefDrag(true); }}
        onDragLeave={() => setRefDrag(false)}
        onDrop={handleRefDrop}
        onClick={() => refRef.current?.click()}
        role="button"
        tabIndex={0}
        onKeyDown={(e) => e.key === "Enter" && refRef.current?.click()}
        aria-label="Drop reference file or click to browse"
      >
        <input
          ref={refRef}
          type="file"
          accept=".xlsx,.csv"
          onChange={handleRefChange}
          style={{ display: "none" }}
        />
        <span className="drop-zone__icon" aria-hidden="true">
          {refTick ? "✓" : "↓"}
        </span>
        <span className="drop-zone__text">
          {refTick
            ? "File received!"
            : refDrag
            ? "Release to drop"
            : "Drag file here or click to browse"}
        </span>
        <span className="drop-zone__hint">.xlsx, .csv</span>
      </div>
      {refFile && (
        <div className="file-list">
          <span className="file-chip">{refFile.name}</span>
        </div>
      )}

      {showHelp && (
        <div className="ref-help">
          <p style={{ margin: "0 0 0.6rem" }}>Two column formats are supported:</p>

          <div style={{ marginBottom: "0.75rem" }}>
            <strong>Schema A</strong> — Protein (T7/T7-term) or pUC19/M13R
            <table className="ref-help-table">
              <thead>
                <tr>
                  <th>Clone number</th>
                  <th>DNA sequence</th>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <td>sample1</td>
                  <td>ATGCATCACC…</td>
                </tr>
                <tr>
                  <td>sample2</td>
                  <td>ATGAAAGT…</td>
                </tr>
              </tbody>
            </table>
            <a href="/template_dna.xlsx" download className="ref-help-dl">
              ↓ Download template
            </a>
          </div>

          <div style={{ marginBottom: "0.5rem" }}>
            <strong>Schema B</strong> — Protein (T7/T7-term) with "AA references" checked
            <table className="ref-help-table">
              <thead>
                <tr>
                  <th>Clone number</th>
                  <th>AA sequence</th>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <td>sample1</td>
                  <td>MHHHHHHHMQ…</td>
                </tr>
                <tr>
                  <td>sample2</td>
                  <td>MKVLIGFTKY…</td>
                </tr>
              </tbody>
            </table>
            <a href="/template_aa.xlsx" download className="ref-help-dl">
              ↓ Download template
            </a>
          </div>

          <p className="ref-help-note">
            <strong>Clone number</strong> must match your filename prefix before the
            first <code>_</code> (case-insensitive).<br />
            Example: <code>sample1_T7.seq</code> → Clone number <code>sample1</code>
          </p>
        </div>
      )}
    </div>
  );
}
