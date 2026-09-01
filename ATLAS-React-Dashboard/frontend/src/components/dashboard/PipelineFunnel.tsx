import type { FunnelItem } from "../../types";

const toneClasses = {
  primary: "border-primary/20 bg-primary/5 text-primary dark:border-slate-600 dark:bg-slate-800 dark:text-white",
  secondary: "border-secondary/30 bg-secondary/10 text-secondary dark:text-emerald-300",
  warning: "border-warning/40 bg-warning/15 text-amber-800 dark:text-amber-300",
  muted: "border-outline/30 bg-surface-container text-text-muted dark:border-slate-600 dark:bg-slate-800 dark:text-slate-300"
};

export function PipelineFunnel({ items }: { items: FunnelItem[] }) {
  return (
    <section className="panel overflow-hidden">
      <div className="border-b border-outline-variant p-4 dark:border-slate-700"><h3 className="text-lg font-semibold text-primary dark:text-white">Pipeline Attrition</h3></div>
      <div className="flex min-h-80 flex-col items-center justify-center gap-2 p-6">
        {items.map(item => <div key={item.label} className={`flex h-12 items-center justify-between rounded-sm border px-4 font-mono text-sm ${toneClasses[item.tone]}`} style={{ width: `${Math.max(34, Math.min(100, item.ratio * 100))}%` }}><span>{item.label}</span><strong>{item.value}</strong></div>)}
      </div>
    </section>
  );
}
