import { useQuery } from "@tanstack/react-query";
import { api } from "../lib/api";
import { MetricCard } from "../components/common/MetricCard";
import { ErrorState, LoadingState } from "../components/common/PageState";
import { PipelineFunnel } from "../components/dashboard/PipelineFunnel";
import { TopCandidatesTable } from "../components/dashboard/TopCandidatesTable";
import { ActivityFeed } from "../components/dashboard/ActivityFeed";

export function DashboardPage() {
  const settings = useQuery({ queryKey: ["settings"], queryFn: api.settings, staleTime: 30_000 });
  const query = useQuery({
    queryKey: ["dashboard"],
    queryFn: api.dashboard,
    refetchInterval: settings.data?.auto_refresh ? (settings.data.refresh_seconds * 1000) : false
  });

  if (query.isLoading || settings.isLoading) return <LoadingState />;
  if (query.error) return <ErrorState error={query.error as Error} />;
  if (settings.error) return <ErrorState error={settings.error as Error} />;
  if (!query.data) return null;

  return (
    <div className="space-y-6">
      <div className="flex flex-col justify-between gap-3 md:flex-row md:items-end">
        <div>
          <p className="label-caps text-secondary">Live pipeline state</p>
          <h1 className="mt-1 text-2xl font-semibold tracking-tight text-primary dark:text-white">Computational evidence overview</h1>
          <p className="mt-2 text-sm text-text-muted dark:text-slate-400">Read-focused status of the current ATLAS outputs under the configured project root.</p>
        </div>
        <div className="text-xs text-text-muted dark:text-slate-400">
          Last updated: <span className="font-mono">{query.data.project.last_updated ?? "Unknown"}</span>
        </div>
      </div>

      {query.data.warnings.length > 0 && (
        <div className="rounded-md border border-warning/40 bg-warning/10 px-4 py-3 text-sm text-amber-900 dark:text-amber-200">{query.data.warnings.join(" • ")}</div>
      )}

      {settings.data?.show_scientific_guardrails && (
        <div className="rounded-md border border-secondary/30 bg-secondary/5 px-4 py-3 text-sm leading-6 text-text-muted dark:text-slate-300">
          <strong className="text-secondary">Interpretation:</strong> CMap, target/network, docking, and ADMET layers are prioritization evidence. They do not establish resistance reversal or clinical efficacy without experimental validation.
        </div>
      )}

      <div className="grid grid-cols-1 gap-4 md:grid-cols-2 xl:grid-cols-4">
        {query.data.metrics.map((metric) => <MetricCard key={metric.key} metric={metric} />)}
      </div>

      <div className="grid grid-cols-1 gap-6 xl:grid-cols-12">
        <div className="xl:col-span-5"><PipelineFunnel items={query.data.funnel} /></div>
        <div className="xl:col-span-7"><TopCandidatesTable rows={query.data.top_candidates} /></div>
      </div>

      <div className="grid grid-cols-1 gap-6 xl:grid-cols-12">
        <div className="xl:col-span-8"><ActivityFeed rows={query.data.activity} /></div>
        <section className="panel p-4 xl:col-span-4">
          <h3 className="text-lg font-semibold text-primary dark:text-white">Primary Validation</h3>
          <p className="mt-1 text-sm text-text-muted dark:text-slate-400">Cohorts currently accepted by the independence-aware gate.</p>
          <div className="mt-4 space-y-2">
            {query.data.primary_validation.length === 0 && <p className="text-sm text-text-muted">No primary-validation datasets are currently exposed by the catalog.</p>}
            {query.data.primary_validation.map((dataset) => <div key={dataset} className="rounded-md border border-outline-variant px-3 py-2 font-mono text-sm dark:border-slate-700">{dataset}</div>)}
          </div>
        </section>
      </div>
    </div>
  );
}
