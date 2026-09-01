import { useQuery } from "@tanstack/react-query";
import { api } from "../lib/api";
import { ErrorState, LoadingState } from "../components/common/PageState";
import { MetricCard } from "../components/common/MetricCard";
import { CandidateScoreChart, DatasetCompositionChart, PathwayChart, TgfbChart, VolcanoChart } from "../components/statistics/ResearchCharts";
import type { Metric } from "../types";

export function ResearchStatisticsPage() {
  const settings = useQuery({ queryKey: ["settings"], queryFn: api.settings, staleTime: 30_000 });
  const query = useQuery({
    queryKey: ["research-statistics"],
    queryFn: api.researchStatistics,
    refetchInterval: settings.data?.auto_refresh ? settings.data.refresh_seconds * 1000 : false
  });

  if (query.isLoading || settings.isLoading) return <LoadingState />;
  if (query.error) return <ErrorState error={query.error as Error} />;
  if (settings.error) return <ErrorState error={settings.error as Error} />;
  if (!query.data) return null;

  const s = query.data.summary;
  const deg = query.data.deg.stats;
  const cards: Metric[] = [
    { key: "datasets", label: "Catalog Datasets", value: s.datasets ?? 0, suffix: "datasets", status: "Current", tone: "success" },
    { key: "primary", label: "Primary Validation", value: s.primary_validation_datasets ?? 0, suffix: "datasets", status: "Gated", tone: "success" },
    { key: "deg", label: "Strict DEGs", value: deg.strict_abs_log2fc_ge_1 ?? 0, suffix: "FDR < 0.05, |LFC| ≥ 1", status: query.data.deg.available ? "Computed" : "Unavailable", tone: query.data.deg.available ? "success" : "muted" },
    { key: "candidates", label: "Integrated Evidence", value: s.integrated_candidates ?? 0, suffix: "candidate rows", status: query.data.candidates.available ? "Available" : "Waiting", tone: query.data.candidates.available ? "success" : "muted" }
  ];

  return (
    <div className="space-y-6">
      <div>
        <p className="label-caps text-secondary">Research analytics</p>
        <h1 className="mt-1 text-2xl font-semibold tracking-tight text-primary dark:text-white">Statistics & research charts</h1>
        <p className="mt-2 max-w-4xl text-sm leading-6 text-text-muted dark:text-slate-400">This page reads current ATLAS outputs and builds exploratory/diagnostic visualizations. It does not turn computational evidence into clinical claims; every chart should be interpreted within the evidence layer that produced it.</p>
      </div>

      <div className="grid grid-cols-1 gap-4 md:grid-cols-2 xl:grid-cols-4">{cards.map((m) => <MetricCard key={m.key} metric={m} />)}</div>

      <div className="grid grid-cols-1 gap-6 xl:grid-cols-12">
        <div className="xl:col-span-8"><VolcanoChart data={query.data.deg} /></div>
        <div className="xl:col-span-4"><DatasetCompositionChart data={query.data.dataset_categories} /></div>
      </div>

      <div className="grid grid-cols-1 gap-6 xl:grid-cols-2"><PathwayChart data={query.data.pathways} /><CandidateScoreChart data={query.data.candidates} /></div>
      <TgfbChart data={query.data.tgfb} />

      {settings.data?.show_scientific_guardrails && (
        <section className="panel p-5">
          <h3 className="text-lg font-semibold">Interpretation guardrails</h3>
          <div className="mt-4 grid gap-3 md:grid-cols-2">{query.data.notes.map((note) => <div key={note} className="rounded-md border border-outline-variant bg-surface-low p-3 text-sm leading-6 text-text-muted dark:border-slate-700 dark:bg-slate-900 dark:text-slate-300">{note}</div>)}</div>
        </section>
      )}
    </div>
  );
}
