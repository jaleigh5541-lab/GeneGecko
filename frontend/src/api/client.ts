const API_BASE = import.meta.env.PUBLIC_API_URL ?? "";

// Wake the backend on page load so it's ready when the user clicks Process
fetch(`${API_BASE}/api/health`).catch(() => {});

export interface CategoriesResponse {
  categories: string[];
  colors: Record<string, string>;
  orf_min_aa: number;
}

export async function fetchCategories(): Promise<CategoriesResponse> {
  const res = await fetch(`${API_BASE}/api/categories`);
  if (!res.ok) throw new Error("Failed to fetch categories");
  return res.json() as Promise<CategoriesResponse>;
}

export interface AlignResponse {
  success: boolean;
  alignment_html: string;
  alignment_html_with_orfs: string | null;
  identity_percent: number;
  orf_count: number;
  longest_orf_length: number;
  longest_orf_seq: string;
}

export async function submitAlign(
  seq1: string,
  seq2: string,
  seqType: "dna" | "protein",
): Promise<AlignResponse> {
  const res = await fetch(`${API_BASE}/api/align`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ seq1, seq2, seq_type: seqType }),
  });
  const data = (await res.json()) as AlignResponse & { detail?: string };
  if (!res.ok || !data.success) {
    throw new Error(data.detail ?? "Alignment failed");
  }
  return data;
}

export type ProgressCallback = (processed: number, total: number) => void;

async function pollJob(
  jobId: string,
  onProgress?: ProgressCallback,
): Promise<unknown[]> {
  while (true) {
    await new Promise((r) => setTimeout(r, 2000));
    const res = await fetch(`${API_BASE}/api/jobs/${jobId}`);
    if (!res.ok) throw new Error(`Job poll failed (HTTP ${res.status})`);
    const data = await res.json();
    if (data.status === "done") return data.results ?? [];
    if (data.status === "error") throw new Error(data.error ?? "Processing failed");
    if (data.progress && onProgress) {
      onProgress(data.progress.processed, data.progress.total);
    }
  }
}

export async function submitProcess(
  formData: FormData,
  onProgress?: ProgressCallback,
): Promise<unknown[]> {
  let res: Response;
  try {
    res = await fetch(`${API_BASE}/api/process`, {
      method: "POST",
      body: formData,
    });
  } catch (err) {
    throw new Error(`Network error: ${err instanceof Error ? err.message : "failed to reach server"}. The server may be starting up — wait 30s and retry.`);
  }
  const text = await res.text();
  if (!text) {
    throw new Error(`Empty response (HTTP ${res.status}). The server may be starting up — try again in 30 seconds.`);
  }
  let data: { success: boolean; job_id?: string; error?: string; detail?: string };
  try {
    data = JSON.parse(text);
  } catch {
    throw new Error(`Server returned non-JSON (HTTP ${res.status}): ${text.slice(0, 200)}`);
  }
  if (!res.ok || !data.success || !data.job_id) {
    throw new Error(data.error ?? data.detail ?? `Processing failed (HTTP ${res.status})`);
  }
  return pollJob(data.job_id, onProgress);
}
