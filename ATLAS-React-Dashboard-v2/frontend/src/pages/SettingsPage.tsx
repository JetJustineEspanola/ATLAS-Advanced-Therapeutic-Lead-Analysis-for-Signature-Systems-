import { useEffect, useState } from "react";
import type { FormEvent } from "react";
import { useMutation, useQuery, useQueryClient } from "@tanstack/react-query";
import { CheckCircle2, FlaskConical, Save, Server, Settings2, ShieldCheck } from "lucide-react";
import { api } from "../lib/api";
import type { DashboardSettings } from "../types";
import { ErrorState, LoadingState } from "../components/common/PageState";

function Toggle({ checked, onChange, label, description }: { checked: boolean; onChange: (v: boolean) => void; label: string; description: string }) {
  return (
    <label className="flex cursor-pointer items-start justify-between gap-4 rounded-md border border-outline-variant p-3 dark:border-slate-700">
      <div><p className="text-sm font-medium text-primary dark:text-white">{label}</p><p className="mt-1 text-xs leading-5 text-text-muted dark:text-slate-400">{description}</p></div>
      <input type="checkbox" checked={checked} onChange={(e) => onChange(e.target.checked)} className="mt-1 h-4 w-4 rounded border-outline-variant text-secondary focus:ring-secondary" />
    </label>
  );
}

export function SettingsPage() {
  const queryClient = useQueryClient();
  const settingsQuery = useQuery({ queryKey: ["settings"], queryFn: api.settings });
  const healthQuery = useQuery({ queryKey: ["health"], queryFn: api.health });
  const [form, setForm] = useState<DashboardSettings | null>(null);
  const [saved, setSaved] = useState(false);

  useEffect(() => { if (settingsQuery.data) setForm(settingsQuery.data); }, [settingsQuery.data]);

  const mutation = useMutation({
    mutationFn: api.saveSettings,
    onSuccess: (data) => {
      setForm(data.settings);
      setSaved(true);
      window.setTimeout(() => setSaved(false), 2500);
      queryClient.invalidateQueries();
    }
  });

  if (settingsQuery.isLoading || healthQuery.isLoading || !form) return <LoadingState />;
  if (settingsQuery.error) return <ErrorState error={settingsQuery.error as Error} />;
  if (healthQuery.error) return <ErrorState error={healthQuery.error as Error} />;

  function set<K extends keyof DashboardSettings>(key: K, value: DashboardSettings[K]) {
    setForm((prev) => prev ? { ...prev, [key]: value } : prev);
  }

  function submit(e: FormEvent) {
    e.preventDefault();
    setSaved(false);
    if (form) mutation.mutate(form);
  }

  return (
    <form onSubmit={submit} className="space-y-6">
      <div className="flex flex-col justify-between gap-3 md:flex-row md:items-end">
        <div>
          <p className="label-caps text-secondary">Configuration</p>
          <h1 className="mt-1 text-2xl font-semibold tracking-tight">Dashboard settings</h1>
          <p className="mt-2 text-sm text-text-muted dark:text-slate-400">Only non-secret dashboard configuration is editable here. API keys and credentials are never returned to the frontend.</p>
        </div>
        <button type="submit" disabled={mutation.isPending} className="flex items-center justify-center gap-2 rounded-md bg-secondary px-4 py-2.5 text-xs font-bold uppercase tracking-wider text-white disabled:opacity-50">
          <Save className="h-4 w-4" /> {mutation.isPending ? "Saving…" : "Save settings"}
        </button>
      </div>

      {saved && <div className="flex items-center gap-2 rounded-md border border-secondary/30 bg-secondary/10 px-4 py-3 text-sm text-secondary"><CheckCircle2 className="h-4 w-4" /> Settings saved and active.</div>}
      {mutation.error && <div className="rounded-md border border-red-300 bg-red-50 px-4 py-3 text-sm text-red-800 dark:border-red-900 dark:bg-red-950/30 dark:text-red-300">{(mutation.error as Error).message}</div>}

      <div className="grid grid-cols-1 gap-6 xl:grid-cols-2">
        <section className="panel p-5">
          <div className="flex items-center gap-2"><Server className="h-5 w-5 text-secondary" /><h2 className="text-lg font-semibold">ATLAS connection</h2></div>
          <p className="mt-2 text-sm leading-6 text-text-muted dark:text-slate-400">The backend reads research outputs from this directory. Changing it takes effect on subsequent API requests.</p>

          <label className="mt-5 block">
            <span className="label-caps">ATLAS_ROOT</span>
            <input value={form.atlas_root} onChange={(e) => set("atlas_root", e.target.value)} className="mt-2 w-full rounded-md border border-outline-variant bg-surface-lowest px-3 py-2 font-mono text-sm focus:border-secondary focus:ring-secondary dark:border-slate-700 dark:bg-slate-900" />
          </label>

          <div className="mt-4 rounded-md border border-outline-variant bg-surface-low p-3 text-sm dark:border-slate-700 dark:bg-slate-900">
            <div className="flex items-center justify-between"><span>API health</span><strong className={healthQuery.data?.ok ? "text-secondary" : "text-red-600"}>{healthQuery.data?.ok ? "ONLINE" : "ROOT NOT FOUND"}</strong></div>
            <div className="mt-2 flex items-center justify-between gap-4"><span>API version</span><code>{healthQuery.data?.api_version}</code></div>
            <div className="mt-2 break-all text-xs text-text-muted dark:text-slate-500">Current backend root: {healthQuery.data?.atlas_root}</div>
          </div>
        </section>

        <section className="panel p-5">
          <div className="flex items-center gap-2"><Settings2 className="h-5 w-5 text-secondary" /><h2 className="text-lg font-semibold">Performance & display</h2></div>
          <div className="mt-5 grid gap-4 sm:grid-cols-2">
            <label><span className="label-caps">Refresh interval</span><div className="mt-2 flex items-center gap-2"><input type="number" min={5} max={3600} value={form.refresh_seconds} onChange={(e) => set("refresh_seconds", Number(e.target.value))} className="w-full rounded-md border border-outline-variant bg-surface-lowest px-3 py-2 font-mono text-sm dark:border-slate-700 dark:bg-slate-900" /><span className="text-xs text-text-muted">sec</span></div></label>
            <label><span className="label-caps">Table / chart limit</span><input type="number" min={100} max={10000} step={100} value={form.table_row_limit} onChange={(e) => set("table_row_limit", Number(e.target.value))} className="mt-2 w-full rounded-md border border-outline-variant bg-surface-lowest px-3 py-2 font-mono text-sm dark:border-slate-700 dark:bg-slate-900" /></label>
          </div>
          <div className="mt-4 space-y-3">
            <Toggle checked={form.auto_refresh} onChange={(v) => set("auto_refresh", v)} label="Auto refresh" description="Allow pages to refresh automatically using the configured interval." />
            <Toggle checked={form.dense_tables} onChange={(v) => set("dense_tables", v)} label="Dense tables" description="Use tighter row spacing for large research tables. Pages can progressively adopt this preference." />
          </div>
        </section>

        <section className="panel p-5">
          <div className="flex items-center gap-2"><FlaskConical className="h-5 w-5 text-secondary" /><h2 className="text-lg font-semibold">Research behavior</h2></div>
          <div className="mt-5 space-y-3">
            <Toggle checked={form.show_scientific_guardrails} onChange={(v) => set("show_scientific_guardrails", v)} label="Show scientific guardrails" description="Keep interpretation warnings visible beside research statistics and evidence layers." />
            <Toggle checked={form.developer_mode} onChange={(v) => set("developer_mode", v)} label="Developer mode" description="Unlock runtime, filesystem, queue-service, output-registry, and API diagnostics. This remains read-only and does not expose shell execution." />
          </div>
        </section>

        <section className="panel p-5">
          <div className="flex items-center gap-2"><ShieldCheck className="h-5 w-5 text-secondary" /><h2 className="text-lg font-semibold">Safety boundary</h2></div>
          <p className="mt-3 text-sm leading-6 text-text-muted dark:text-slate-400">The web interface is intentionally read-focused. Settings can change dashboard behavior and the data root, but they cannot execute arbitrary terminal commands, alter scientific thresholds in pipeline scripts, or silently start expensive CMap/docking jobs.</p>
          <div className="mt-4 rounded-md border border-secondary/30 bg-secondary/5 p-3 text-sm"><strong>Recommended next control layer:</strong> a dedicated job API with an allowlist of pipeline actions, job IDs, log streaming, cancellation, disk guards, and explicit confirmations for expensive stages.</div>
        </section>
      </div>
    </form>
  );
}
