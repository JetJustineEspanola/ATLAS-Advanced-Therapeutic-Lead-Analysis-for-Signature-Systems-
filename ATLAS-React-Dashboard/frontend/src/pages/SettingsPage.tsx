import { useQuery } from "@tanstack/react-query";
import { api } from "../lib/api";
import { ErrorState, LoadingState } from "../components/common/PageState";

export function SettingsPage() {
  const query = useQuery({ queryKey: ["health"], queryFn: api.health });
  if (query.isLoading) return <LoadingState />;
  if (query.error) return <ErrorState error={query.error as Error} />;
  return <div className="grid grid-cols-1 gap-6 xl:grid-cols-2">
    <section className="panel p-5"><p className="label-caps text-secondary">Backend</p><h3 className="mt-2 text-lg font-semibold">ATLAS data connection</h3><dl className="mt-5 space-y-4 text-sm"><div><dt className="text-text-muted dark:text-slate-400">Health</dt><dd className="mt-1 font-mono">{query.data?.ok ? "ONLINE" : "OFFLINE"}</dd></div><div><dt className="text-text-muted dark:text-slate-400">ATLAS_ROOT</dt><dd className="mt-1 break-all font-mono">{query.data?.atlas_root}</dd></div></dl></section>
    <section className="panel p-5"><p className="label-caps text-secondary">Architecture</p><h3 className="mt-2 text-lg font-semibold">Scalable by default</h3><p className="mt-3 text-sm leading-6 text-text-muted dark:text-slate-400">React handles presentation and routing. FastAPI provides a stable API boundary over ATLAS outputs, so local CSV readers can later be replaced by DuckDB, PostgreSQL, object storage, or remote compute without rewriting the UI.</p></section>
  </div>;
}
