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
