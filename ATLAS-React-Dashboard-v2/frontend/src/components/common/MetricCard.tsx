import type { Metric } from "../../types";
import { StatusBadge } from "./StatusBadge";

export function MetricCard({ metric }: { metric: Metric }) {
  return (
    <section className="panel flex min-h-32 flex-col p-4">
      <div className="mb-3 flex items-start justify-between gap-3">
        <span className="label-caps text-primary dark:text-slate-100">{metric.label}</span>
        <StatusBadge text={metric.status} tone={metric.tone} />
      </div>
      <div className="mt-auto flex items-end justify-between gap-3">
        <span className="text-3xl font-semibold tracking-tight text-primary dark:text-white">{metric.value}</span>
        {metric.suffix && <span className="pb-1 text-xs text-text-muted dark:text-slate-400">{metric.suffix}</span>}
      </div>
      {typeof metric.progress === "number" && (
        <div className="mt-3 h-1 overflow-hidden rounded-full bg-surface-container dark:bg-slate-700">
          <div className="h-full bg-warning" style={{ width: `${Math.max(0, Math.min(100, metric.progress))}%` }} />
        </div>
      )}
    </section>
  );
}
