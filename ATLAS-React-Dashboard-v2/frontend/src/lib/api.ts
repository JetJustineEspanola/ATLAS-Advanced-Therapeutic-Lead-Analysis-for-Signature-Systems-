import type {
  DashboardPayload,
  DashboardSettings,
  DatasetRow,
  DeveloperStatistics,
  GenericRows,
  ResearchStatistics
} from "../types";

const API_BASE = (import.meta.env.VITE_API_URL ?? "").replace(/\/$/, "");

async function request<T>(url: string, init?: RequestInit): Promise<T> {
  const response = await fetch(`${API_BASE}${url}`, init);
  if (!response.ok) {
    let message = `Request failed: ${response.status}`;
    try {
      const data = await response.json();
      message = data.detail ?? JSON.stringify(data);
    } catch {
      const text = await response.text();
      if (text) message = text;
    }
    throw new Error(message);
  }
  return response.json() as Promise<T>;
}

export const api = {
  health: () => request<{ ok: boolean; atlas_root: string; api_version: string; settings_file: string }>("/api/health"),
  dashboard: () => request<DashboardPayload>("/api/dashboard"),
  datasets: () => request<{ rows: DatasetRow[]; columns: string[]; source_file?: string }>("/api/datasets"),
  signature: () => request<GenericRows>("/api/signature"),
  candidates: () => request<GenericRows>("/api/candidates"),
  cmap: () => request<GenericRows>("/api/cmap"),
  docking: () => request<GenericRows>("/api/docking"),
  settings: () => request<DashboardSettings>("/api/settings"),
  saveSettings: (settings: DashboardSettings) => request<{ ok: boolean; settings: DashboardSettings }>("/api/settings", {
    method: "PUT",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(settings)
  }),
  researchStatistics: () => request<ResearchStatistics>("/api/research/statistics"),
  developerStatistics: () => request<DeveloperStatistics>("/api/developer/statistics")
};
