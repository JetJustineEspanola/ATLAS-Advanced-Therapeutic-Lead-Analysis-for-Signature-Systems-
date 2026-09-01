import { useMemo, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import { Search } from "lucide-react";
import { api } from "../lib/api";
import { ErrorState, LoadingState } from "../components/common/PageState";

export function DatasetsPage() {
  const [search, setSearch] = useState("");
  const query = useQuery({ queryKey: ["datasets"], queryFn: api.datasets });
  const rows = useMemo(() => { const q = search.trim().toLowerCase(); if (!q) return query.data?.rows ?? []; return (query.data?.rows ?? []).filter(row => Object.values(row).some(v => String(v ?? "").toLowerCase().includes(q))); }, [query.data, search]);
  if (query.isLoading) return <LoadingState />;
  if (query.error) return <ErrorState error={query.error as Error} />;
  return <div className="space-y-4">
    <div className="panel flex items-center gap-3 p-3"><Search className="h-4 w-4 text-text-muted" /><input value={search} onChange={e => setSearch(e.target.value)} placeholder="Filter dataset id, title, source, category…" className="w-full border-0 bg-transparent text-sm outline-none ring-0 focus:ring-0 dark:text-white" /></div>
    <section className="panel overflow-hidden"><div className="border-b border-outline-variant p-4 dark:border-slate-700"><h3 className="text-lg font-semibold">Dataset Catalog</h3><p className="mt-1 text-xs text-text-muted dark:text-slate-400">{rows.length} rows</p></div><div className="max-h-[72vh] overflow-auto"><table className="w-full min-w-[1000px] text-left text-sm"><thead className="sticky top-0 bg-surface dark:bg-slate-900"><tr className="border-b border-outline-variant dark:border-slate-700">{["Dataset", "Source", "Category", "Score", "Modality", "Phenotype", "Relationship", "Title"].map(h => <th key={h} className="p-3 label-caps">{h}</th>)}</tr></thead><tbody className="divide-y divide-outline-variant dark:divide-slate-800">{rows.map((row, i) => <tr key={`${row.dataset_id}-${i}`} className="hover:bg-surface-low dark:hover:bg-slate-900"><td className="p-3 font-mono text-xs">{row.dataset_id}</td><td className="p-3">{row.source ?? "—"}</td><td className="p-3">{row.category ?? "—"}</td><td className="p-3 font-mono">{row.score ?? "—"}</td><td className="p-3">{row.modality ?? "—"}</td><td className="p-3">{row.phenotype_confidence ?? "—"}</td><td className="p-3">{row.relationship_role ?? "—"}</td><td className="max-w-[420px] p-3">{row.title ?? "—"}</td></tr>)}</tbody></table></div></section>
  </div>;
}
