import {
  Bar,
  BarChart,
  CartesianGrid,
  Cell,
  Legend,
  Pie,
  PieChart,
  ReferenceLine,
  ResponsiveContainer,
  Scatter,
  ScatterChart,
  Tooltip,
  XAxis,
  YAxis
} from "recharts";
import type { ResearchStatistics } from "../../types";

const PALETTE = ["#006a61", "#1e293b", "#f59e0b", "#6b7280", "#8b5cf6", "#0ea5e9", "#be123c"];

function EmptyChart({ text }: { text: string }) {
  return <div className="flex h-72 items-center justify-center px-6 text-center text-sm text-text-muted dark:text-slate-400">{text}</div>;
}

export function VolcanoChart({ data }: { data: ResearchStatistics["deg"] }) {
  return (
    <section className="panel overflow-hidden">
      <div className="border-b border-outline-variant p-4 dark:border-slate-700">
        <h3 className="text-lg font-semibold">Differential-expression volcano</h3>
        <p className="mt-1 text-xs text-text-muted dark:text-slate-400">
          x = log2 fold change; y = −log10 adjusted p-value. Display points are performance-limited, while summary counts use the full file.
        </p>
      </div>
      {!data.available ? <EmptyChart text="No compatible differential-expression file was found." /> : (
        <div className="h-[420px] p-3">
          <ResponsiveContainer width="100%" height="100%">
            <ScatterChart margin={{ top: 10, right: 20, bottom: 20, left: 10 }}>
              <CartesianGrid strokeDasharray="3 3" opacity={0.25} />
              <XAxis type="number" dataKey="log2fc" name="log2FC" tick={{ fontSize: 11 }} label={{ value: "log2 fold change", position: "insideBottom", offset: -10 }} />
              <YAxis type="number" dataKey="neglog10_padj" name="−log10(FDR)" tick={{ fontSize: 11 }} width={52} />
              <ReferenceLine x={-1} strokeDasharray="4 4" />
              <ReferenceLine x={1} strokeDasharray="4 4" />
              <Tooltip cursor={{ strokeDasharray: "3 3" }} formatter={(value) => Number(value).toFixed(3)} />
              <Scatter data={data.points.filter((p) => p.class === "other")} fill="#94a3b8" name="Other" />
              <Scatter data={data.points.filter((p) => p.class === "up")} fill="#006a61" name="Up" />
              <Scatter data={data.points.filter((p) => p.class === "down")} fill="#be123c" name="Down" />
              <Legend />
            </ScatterChart>
          </ResponsiveContainer>
        </div>
      )}
    </section>
  );
}

export function DatasetCompositionChart({ data }: { data: ResearchStatistics["dataset_categories"] }) {
  return (
    <section className="panel overflow-hidden">
      <div className="border-b border-outline-variant p-4 dark:border-slate-700">
        <h3 className="text-lg font-semibold">Dataset evidence roles</h3>
        <p className="mt-1 text-xs text-text-muted dark:text-slate-400">Independence-aware dataset classification currently present in the catalog.</p>
      </div>
      {!data.available ? <EmptyChart text="Dataset category information is not available." /> : (
        <div className="h-80 p-3">
          <ResponsiveContainer width="100%" height="100%">
            <PieChart>
              <Pie data={data.rows} dataKey="count" nameKey="category" cx="50%" cy="46%" outerRadius={95} label={({ name, value }) => `${name}: ${value}`}>
                {data.rows.map((_, i) => <Cell key={i} fill={PALETTE[i % PALETTE.length]} />)}
              </Pie>
              <Tooltip />
            </PieChart>
          </ResponsiveContainer>
        </div>
      )}
    </section>
  );
}

export function PathwayChart({ data }: { data: ResearchStatistics["pathways"] }) {
  const rows = [...data.rows].sort((a, b) => a.nes - b.nes);
  return (
    <section className="panel overflow-hidden">
      <div className="border-b border-outline-variant p-4 dark:border-slate-700">
        <h3 className="text-lg font-semibold">Pathway enrichment</h3>
        <p className="mt-1 text-xs text-text-muted dark:text-slate-400">Top pathways ranked by absolute normalized enrichment score from a compatible current output.</p>
      </div>
      {!data.available ? <EmptyChart text="No compatible GSEA/pathway CSV with an NES column was found." /> : (
        <div className="h-[520px] p-3">
          <ResponsiveContainer width="100%" height="100%">
            <BarChart data={rows} layout="vertical" margin={{ top: 5, right: 20, bottom: 5, left: 145 }}>
              <CartesianGrid strokeDasharray="3 3" opacity={0.25} />
              <XAxis type="number" />
              <YAxis type="category" dataKey="pathway" width={140} tick={{ fontSize: 10 }} />
              <Tooltip />
              <ReferenceLine x={0} />
              <Bar dataKey="nes" name="NES">
                {rows.map((row, i) => <Cell key={i} fill={row.nes >= 0 ? "#006a61" : "#be123c"} />)}
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        </div>
      )}
    </section>
  );
}

export function CandidateScoreChart({ data }: { data: ResearchStatistics["candidates"] }) {
  const usable = data.rows.filter((r) => typeof r.score === "number").slice(0, 12).reverse();
  return (
    <section className="panel overflow-hidden">
      <div className="border-b border-outline-variant p-4 dark:border-slate-700">
        <h3 className="text-lg font-semibold">Integrated candidate ranking</h3>
        <p className="mt-1 text-xs text-text-muted dark:text-slate-400">Integrated evidence score where the current output exposes a compatible score column.</p>
      </div>
      {!data.available || usable.length === 0 ? <EmptyChart text="No numeric integrated candidate score was detected in the current evidence matrix." /> : (
        <div className="h-[420px] p-3">
          <ResponsiveContainer width="100%" height="100%">
            <BarChart data={usable} layout="vertical" margin={{ top: 5, right: 20, bottom: 5, left: 105 }}>
              <CartesianGrid strokeDasharray="3 3" opacity={0.25} />
              <XAxis type="number" />
              <YAxis type="category" dataKey="candidate" width={100} tick={{ fontSize: 10 }} />
              <Tooltip />
              <Bar dataKey="score" fill="#006a61" name="Integrated score" />
            </BarChart>
          </ResponsiveContainer>
        </div>
      )}
    </section>
  );
}

export function TgfbChart({ data }: { data: ResearchStatistics["tgfb"] }) {
  return (
    <section className="panel overflow-hidden">
      <div className="border-b border-outline-variant p-4 dark:border-slate-700">
        <h3 className="text-lg font-semibold">TGF-beta ranked validation</h3>
        <p className="mt-1 text-xs text-text-muted dark:text-slate-400">Direction and magnitude of NES values if a compatible TGF-beta validation table is available.</p>
      </div>
      {!data.available ? <EmptyChart text="No compatible TGF-beta ranked-validation CSV was auto-detected." /> : (
        <div className="h-80 p-3">
          <ResponsiveContainer width="100%" height="100%">
            <BarChart data={data.rows} margin={{ top: 5, right: 20, bottom: 45, left: 10 }}>
              <CartesianGrid strokeDasharray="3 3" opacity={0.25} />
              <XAxis dataKey="dataset" angle={-30} textAnchor="end" height={75} tick={{ fontSize: 10 }} />
              <YAxis />
              <ReferenceLine y={0} />
              <Tooltip />
              <Bar dataKey="nes" name="NES">
                {data.rows.map((row, i) => <Cell key={i} fill={row.nes >= 0 ? "#006a61" : "#be123c"} />)}
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        </div>
      )}
    </section>
  );
}
