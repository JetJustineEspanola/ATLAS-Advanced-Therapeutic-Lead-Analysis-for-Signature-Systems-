export type StatusTone = "success" | "warning" | "muted" | "danger";

export interface Metric {
  key: string;
  label: string;
  value: string | number;
  suffix?: string;
  status: string;
  tone: StatusTone;
  progress?: number | null;
}

export interface FunnelItem {
  label: string;
  value: string | number;
  ratio: number;
  tone: "primary" | "secondary" | "warning" | "muted";
}

export interface Candidate {
  name: string;
  connectivity_score?: number | string | null;
  docking_score?: number | string | null;
  status?: string;
  target?: string | null;
  final_score?: number | string | null;
  priority?: string | null;
}

export interface Activity {
  message: string;
  timestamp: string;
  tone: StatusTone;
}

export interface DashboardPayload {
  project: { name: string; atlas_root: string; last_updated: string | null };
  metrics: Metric[];
  funnel: FunnelItem[];
  top_candidates: Candidate[];
  activity: Activity[];
  primary_validation: string[];
  warnings: string[];
}

export interface DatasetRow {
  dataset_id: string;
  source?: string;
  title?: string;
  category?: string;
  score?: number | string;
  modality?: string;
  phenotype_confidence?: string;
  relationship_role?: string;
  sample_count?: number | string;
}

export interface GenericRows {
  rows: Record<string, unknown>[];
  columns: string[];
  source_file?: string | null;
}

export interface DashboardSettings {
  atlas_root: string;
  auto_refresh: boolean;
  refresh_seconds: number;
  table_row_limit: number;
  developer_mode: boolean;
  show_scientific_guardrails: boolean;
  dense_tables: boolean;
}

export interface ResearchStatistics {
  summary: Record<string, number>;
  deg: {
    available: boolean;
    source_file?: string;
    points: Array<{ gene: string; log2fc: number; neglog10_padj: number; padj: number; class: string }>;
    stats: Record<string, number>;
  };
  dataset_categories: { available: boolean; rows: Array<{ category: string; count: number }> };
  pathways: { available: boolean; rows: Array<{ pathway: string; nes: number; fdr?: number | null }>; source_file?: string | null };
  candidates: { available: boolean; rows: Array<{ candidate: string; score?: number | null; tau?: number | null; docking?: number | null; target?: string | null }>; source_file?: string | null };
  tgfb: { available: boolean; rows: Array<{ dataset: string; nes: number; fdr?: number | null }>; source_file?: string | null };
  notes: string[];
}

export interface DeveloperStatistics {
  runtime: { python: string; platform: string; pid: number; cwd: string; atlas_root: string };
  filesystem: { total_gb: number; free_gb: number; used_gb: number; free_percent: number | null };
  queue_service: { available: boolean; active_state: string; sub_state: string; main_pid: string };
  dataset_queue: { source_file: string; exists: boolean; rows: number; status_counts: Record<string, number>; recent: Record<string, string>[]; error?: string };
  outputs: Array<{ path: string; exists: boolean; size_mb: number | null; modified: string | null }>;
  api: string[];
}
