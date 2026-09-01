import type { StatusTone } from "../../types";

const tones: Record<StatusTone, string> = {
  success: "bg-secondary/10 text-secondary dark:bg-emerald-400/10 dark:text-emerald-300",
  warning: "bg-warning/10 text-amber-800 dark:text-amber-300",
  muted: "bg-surface-container text-text-muted dark:bg-slate-800 dark:text-slate-300",
  danger: "bg-red-50 text-danger dark:bg-red-950/40 dark:text-red-300"
};

export function StatusBadge({ text, tone = "muted" }: { text: string; tone?: StatusTone }) {
  return <span className={`inline-flex items-center rounded-full px-2 py-1 text-xs font-medium ${tones[tone]}`}>{text}</span>;
}
