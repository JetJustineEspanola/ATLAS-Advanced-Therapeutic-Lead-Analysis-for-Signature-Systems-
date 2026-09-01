import { useQuery } from "@tanstack/react-query";
import type { GenericRows } from "../types";
import { ErrorState, LoadingState } from "../components/common/PageState";

export function GenericTablePage({ queryKey, title, subtitle, loader }: { queryKey: string; title: string; subtitle: string; loader: () => Promise<GenericRows> }) {
  const query = useQuery({ queryKey: [queryKey], queryFn: loader });
  if (query.isLoading) return <LoadingState />;
  if (query.error) return <ErrorState error={query.error as Error} />;
  const rows = query.data?.rows ?? [];
  const columns = (query.data?.columns ?? []).slice(0, 12);
  return <section className="panel overflow-hidden">
    <div className="border-b border-outline-variant p-4 dark:border-slate-700"><h3 className="text-lg font-semibold text-primary dark:text-white">{title}</h3><p className="mt-1 text-sm text-text-muted dark:text-slate-400">{subtitle}</p>{query.data?.source_file && <p className="mt-2 font-mono text-[11px] text-text-muted dark:text-slate-500">{query.data.source_file}</p>}</div>
    <div className="max-h-[76vh] overflow-auto"><table className="w-full min-w-[900px] text-left text-sm"><thead className="sticky top-0 bg-surface dark:bg-slate-900"><tr className="border-b border-outline-variant dark:border-slate-700">{columns.map(c => <th key={c} className="p-3 label-caps">{c}</th>)}</tr></thead><tbody className="divide-y divide-outline-variant dark:divide-slate-800">{rows.length === 0 && <tr><td colSpan={Math.max(columns.length, 1)} className="p-10 text-center text-text-muted">No data available yet.</td></tr>}{rows.slice(0, 1000).map((row, i) => <tr key={i} className="hover:bg-surface-low dark:hover:bg-slate-900">{columns.map(c => <td key={c} className="max-w-[320px] truncate p-3 font-mono text-xs" title={String(row[c] ?? "")}>{String(row[c] ?? "—")}</td>)}</tr>)}</tbody></table></div>
  </section>;
}
