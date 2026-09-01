import type { DashboardPayload, DatasetRow, GenericRows } from "../types";

async function get<T>(url: string): Promise<T> {
  const response = await fetch(url);
  if (!response.ok) {
    const body = await response.text();
    throw new Error(body || `Request failed: ${response.status}`);
  }
  return response.json() as Promise<T>;
}

export const api = {
  health: () => get<{ ok: boolean; atlas_root: string }>("/api/health"),
  dashboard: () => get<DashboardPayload>("/api/dashboard"),
  datasets: () => get<{ rows: DatasetRow[]; columns: string[]; source_file?: string }>("/api/datasets"),
  signature: () => get<GenericRows>("/api/signature"),
  candidates: () => get<GenericRows>("/api/candidates"),
  cmap: () => get<GenericRows>("/api/cmap"),
  docking: () => get<GenericRows>("/api/docking")
};
