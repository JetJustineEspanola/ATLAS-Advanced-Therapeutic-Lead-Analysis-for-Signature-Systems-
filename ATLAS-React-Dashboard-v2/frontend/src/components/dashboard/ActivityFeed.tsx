import { CheckCircle2, CircleAlert, Clock3, RefreshCw } from "lucide-react";
import type { Activity } from "../../types";

const icons = { success: CheckCircle2, warning: RefreshCw, danger: CircleAlert, muted: Clock3 };

export function ActivityFeed({ rows }: { rows: Activity[] }) {
  return <section className="panel p-4"><h3 className="mb-4 text-lg font-semibold text-primary dark:text-white">Recent Activity</h3><div className="space-y-1">{rows.map((row, i) => { const Icon = icons[row.tone]; return <div key={`${row.timestamp}-${i}`} className="flex gap-3 rounded-md p-3 hover:bg-surface-low dark:hover:bg-slate-900"><Icon className="mt-0.5 h-4 w-4 shrink-0 text-secondary" /><div><p className="text-sm text-primary dark:text-slate-100">{row.message}</p><p className="mt-1 font-mono text-[11px] text-text-muted dark:text-slate-500">{row.timestamp}</p></div></div>; })}</div></section>;
}
