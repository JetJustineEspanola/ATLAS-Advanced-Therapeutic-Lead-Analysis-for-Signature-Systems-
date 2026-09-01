import { useMemo, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import { Search } from "lucide-react";
import { api } from "../lib/api";
import { ErrorState, LoadingState } from "../components/common/PageState";

export function DatasetsPage() {
  const [search, setSearch] = useState("");
  const query = useQuery({ queryKey: ["datasets"], queryFn: api.datasets });
  const settings = useQuery({ queryKey: ["settings"], queryFn: api.settings, staleTime: 30_000 });

  const rows = useMemo(() => {
    const q = search.trim().toLowerCase();
    const base = query.data?.rows ?? [];
    if (!q) return base;
    return base.filter((row) => Object.values(row).some((v) => String(v ?? "").toLowerCase().includes(q)));
  }, [query.data, search]);

  if (query.isLoading || settings.isLoading) return <LoadingState />;
  if (query.error) return <ErrorState error={query.error as Error} />;
  if (settings.error) return <ErrorState error={settings.error as Error} />;
  const limit = settings.data?.table_row_limit ?? 1000;
  const py = settings.data?.dense_tables ? "py-1.5" : "py-3";

  return (
    <div className="space-y-4">
      <div className="panel flex items-center gap-3 p-3"><Search className="h-4 w-4 text-text-muted" /><input value={search} onChange={(e) => setSearch(e.target.value)} placeholder="Filter dataset id, title, source, category…" className="w-full border-0 bg-transparent text-sm outline-none ring-0 focus:ring-0 dark:text-white" /></div>
      <section className="panel overflow-hidden">
        <div className="border-b border-outline-variant p-4 dark:border-slate-700"><h3 className="text-lg font-semibold">Dataset Catalog</h3><p className="mt-1 text-xs text-text-muted dark:text-slate-400">{rows.length.toLocaleString()} matching rows · displaying up to {Math.min(limit, rows.length).toLocaleString()}</p></div>
        <div className="max-h-[72vh] overflow-auto">
          <table className="w-full min-w-[1000px] text-left text-sm">
            <thead className="sticky top-0 z-10 bg-surface dark:bg-slate-900"><tr className="border-b border-outline-variant dark:border-slate-700">{["Dataset", "Source", "Category", "Score", "Modality", "Phenotype", "Relationship", "Title"].map((h) => <th key={h} className="p-3 label-caps">{h}</th>)}</tr></thead>
            <tbody className="divide-y divide-outline-variant dark:divide-slate-800">
              {rows.slice(0, limit).map((row, i) => <tr key={`${row.dataset_id}-${i}`} className="hover:bg-surface-low dark:hover:bg-slate-900"><td className={`px-3 ${py} font-mono text-xs`}>{row.dataset_id}</td><td className={`px-3 ${py}`}>{row.source ?? "—"}</td><td className={`px-3 ${py}`}>{row.category ?? "—"}</td><td className={`px-3 ${py} font-mono`}>{row.score ?? "—"}</td><td className={`px-3 ${py}`}>{row.modality ?? "—"}</td><td className={`px-3 ${py}`}>{row.phenotype_confidence ?? "—"}</td><td className={`px-3 ${py}`}>{row.relationship_role ?? "—"}</td><td className={`max-w-[420px] px-3 ${py}`}>{row.title ?? "—"}</td></tr>)}
            </tbody>
          </table>
        </div>
      </section>
    </div>
  );
}
