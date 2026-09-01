import { ChevronRight } from "lucide-react";
import type { Candidate } from "../../types";
import { StatusBadge } from "../common/StatusBadge";

function tone(status = "") {
  const s = status.toLowerCase();
  if (s.includes("priority") || s.includes("pass") || s.includes("complete")) return "success" as const;
  if (s.includes("caution") || s.includes("queue") || s.includes("progress")) return "warning" as const;
  if (s.includes("fail") || s.includes("risk")) return "danger" as const;
  return "muted" as const;
}

export function TopCandidatesTable({ rows }: { rows: Candidate[] }) {
  return (
    <section className="panel overflow-hidden">
      <div className="flex items-center justify-between border-b border-outline-variant bg-surface p-4 dark:border-slate-700 dark:bg-slate-900"><h3 className="text-lg font-semibold text-primary dark:text-white">Top Final Candidates</h3><span className="label-caps text-secondary">{rows.length} shown</span></div>
      <div className="overflow-x-auto">
        <table className="w-full min-w-[760px] border-collapse text-left">
          <thead className="border-b border-outline-variant bg-surface dark:border-slate-700 dark:bg-slate-900"><tr><th className="p-3 label-caps">Compound</th><th className="p-3 text-right label-caps">CMap</th><th className="p-3 text-right label-caps">Docking</th><th className="p-3 label-caps">Target</th><th className="p-3 label-caps">Status</th><th className="w-10 p-3" /></tr></thead>
          <tbody className="divide-y divide-outline-variant dark:divide-slate-800">
            {rows.length === 0 && <tr><td colSpan={6} className="p-8 text-center text-sm text-text-muted dark:text-slate-400">No candidate rows are available yet.</td></tr>}
            {rows.map((row, index) => <tr key={`${row.name}-${index}`} className="hover:bg-surface-low dark:hover:bg-slate-900"><td className="p-3 font-mono text-sm font-medium text-primary dark:text-white">{row.name}</td><td className="p-3 text-right font-mono text-sm">{row.connectivity_score ?? "—"}</td><td className="p-3 text-right font-mono text-sm">{row.docking_score ?? "—"}</td><td className="p-3 font-mono text-xs">{row.target ?? "—"}</td><td className="p-3"><StatusBadge text={row.status ?? row.priority ?? "Available"} tone={tone(row.status ?? row.priority)} /></td><td className="p-3 text-right text-text-muted"><ChevronRight className="h-4 w-4" /></td></tr>)}
          </tbody>
        </table>
      </div>
    </section>
  );
}
