export interface CategoriesResponse {
  categories: string[];
  colors: Record<string, string>;
  orf_min_aa: number;
}

export async function fetchCategories(): Promise<CategoriesResponse> {
  const res = await fetch("/api/categories");
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
  const res = await fetch("/api/align", {
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

export async function submitProcess(formData: FormData): Promise<unknown[]> {
  const res = await fetch("/api/process", {
    method: "POST",
    body: formData,
  });
  const text = await res.text();
  if (!text) {
    throw new Error("Server returned an empty response — it may have timed out. Try fewer/smaller files.");
  }
  let data: { success: boolean; results?: unknown[]; error?: string; detail?: string };
  try {
    data = JSON.parse(text);
  } catch {
    throw new Error("Invalid response from server. It may have timed out or crashed.");
  }
  if (!res.ok || !data.success) {
    throw new Error(data.error ?? data.detail ?? "Processing failed");
  }
  return data.results ?? [];
}
