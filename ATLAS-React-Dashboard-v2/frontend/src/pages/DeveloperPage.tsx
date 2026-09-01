import { useQuery } from "@tanstack/react-query";
import { Activity, Database, HardDrive, Server } from "lucide-react";
import { api } from "../lib/api";
import { ErrorState, LoadingState } from "../components/common/PageState";

function Kpi({ label, value, detail, icon: Icon }: { label: string; value: string; detail: string; icon: typeof Server }) {
  return (
    <section className="panel p-4">
      <div className="flex items-start justify-between">
        <div>
          <p className="label-caps text-text-muted dark:text-slate-400">{label}</p>
          <p className="mt-2 text-2xl font-semibold text-primary dark:text-white">{value}</p>
          <p className="mt-1 text-xs text-text-muted dark:text-slate-500">{detail}</p>
        </div>
        <Icon className="h-5 w-5 text-secondary" />
      </div>
    </section>
  );
}

export function DeveloperPage() {
  const query = useQuery({ queryKey: ["developer-statistics"], queryFn: api.developerStatistics, refetchInterval: 15_000 });
  if (query.isLoading) return <LoadingState />;
  if (query.error) return <ErrorState error={query.error as Error} />;
  if (!query.data) return null;

  const d = query.data;
  return (
    <div className="space-y-6">
      <div>
        <p className="label-caps text-secondary">Developer mode</p>
        <h1 className="mt-1 text-2xl font-semibold">Runtime & pipeline diagnostics</h1>
        <p className="mt-2 text-sm text-text-muted dark:text-slate-400">Read-only operational diagnostics for the dashboard, queue service, filesystem, and important ATLAS artifacts.</p>
      </div>

      <div className="grid grid-cols-1 gap-4 md:grid-cols-2 xl:grid-cols-4">
        <Kpi label="Queue Service" value={d.queue_service.active_state} detail={`${d.queue_service.sub_state} · PID ${d.queue_service.main_pid}`} icon={Activity} />
        <Kpi label="Filesystem Free" value={`${d.filesystem.free_gb} GB`} detail={`${d.filesystem.free_percent ?? "—"}% free`} icon={HardDrive} />
        <Kpi label="Queue Rows" value={String(d.dataset_queue.rows)} detail={d.dataset_queue.exists ? "queue state file available" : "queue file missing"} icon={Database} />
        <Kpi label="API Runtime" value={`Python ${d.runtime.python}`} detail={`PID ${d.runtime.pid}`} icon={Server} />
      </div>

      <div className="grid grid-cols-1 gap-6 xl:grid-cols-2">
        <section className="panel overflow-hidden">
          <div className="border-b border-outline-variant p-4 dark:border-slate-700"><h3 className="font-semibold">Important output registry</h3></div>
          <div className="overflow-x-auto">
            <table className="w-full min-w-[620px] text-left text-sm">
              <thead><tr className="border-b border-outline-variant dark:border-slate-700"><th className="p-3 label-caps">Path</th><th className="p-3 label-caps">Exists</th><th className="p-3 label-caps">Size</th><th className="p-3 label-caps">Modified</th></tr></thead>
              <tbody className="divide-y divide-outline-variant dark:divide-slate-800">
                {d.outputs.map((o) => <tr key={o.path}><td className="p-3 font-mono text-xs">{o.path}</td><td className="p-3">{o.exists ? "YES" : "NO"}</td><td className="p-3 font-mono">{o.size_mb == null ? "—" : `${o.size_mb} MB`}</td><td className="p-3 font-mono text-xs">{o.modified ?? "—"}</td></tr>)}
              </tbody>
            </table>
          </div>
        </section>

        <section className="panel p-4">
          <h3 className="font-semibold">Queue status counts</h3>
          <div className="mt-4 space-y-2">
            {Object.keys(d.dataset_queue.status_counts).length === 0 && <p className="text-sm text-text-muted">No status column detected or queue is empty.</p>}
            {Object.entries(d.dataset_queue.status_counts).map(([k, v]) => <div key={k} className="flex items-center justify-between rounded-md border border-outline-variant px-3 py-2 dark:border-slate-700"><span className="font-mono text-xs">{k}</span><strong>{v}</strong></div>)}
          </div>
          <h4 className="mt-6 label-caps">Runtime</h4>
          <pre className="mt-2 overflow-auto rounded-md bg-primary p-3 text-xs text-slate-100">{JSON.stringify(d.runtime, null, 2)}</pre>
        </section>
      </div>

      <section className="panel p-4">
        <h3 className="font-semibold">API surface</h3>
        <div className="mt-4 grid gap-2 md:grid-cols-2 xl:grid-cols-3">
          {d.api.map((e) => <code key={e} className="rounded-md border border-outline-variant bg-surface-low px-3 py-2 text-xs dark:border-slate-700 dark:bg-slate-900">{e}</code>)}
        </div>
      </section>
    </div>
  );
}
