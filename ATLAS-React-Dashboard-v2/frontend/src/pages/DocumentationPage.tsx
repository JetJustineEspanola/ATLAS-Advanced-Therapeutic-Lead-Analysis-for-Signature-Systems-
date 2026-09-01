import type { ReactNode } from "react";
import {
  BookOpen,
  Database,
  FlaskConical,
  GitBranch,
  Server,
  ShieldCheck,
  Terminal,
  Workflow,
} from "lucide-react";

const sources = [
  [
    "NCBI GEO",
    "Primary discovery/validation studies and processed gene-expression metadata.",
  ],
  [
    "NCBI SRA",
    "Raw sequencing runs and accession-level links for transcriptomic studies.",
  ],
  [
    "EMBL-EBI BioStudies / ArrayExpress",
    "External study discovery plus SDRF/sample metadata; E-GEOD records are treated as GEO mirrors, not independent evidence.",
  ],
  ["ENA", "Sequencing cross-references associated with EBI studies."],
  ["BioSamples", "Sample-level metadata context."],
  [
    "GDC / TCGA",
    "Patient/translational cancer-genomics context, not a direct acquired trastuzumab-resistance model by itself.",
  ],
  ["cBioPortal", "Cancer cohort/genomic context and translational support."],
  [
    "CLUE / LINCS CMap",
    "Perturbational signatures used to identify compounds whose transcriptional effects oppose the resistance signature.",
  ],
  ["DepMap", "Planned/optional model and target-dependency context."],
  ["CPTAC / PDC", "Planned/optional proteomic support."],
];

const stages = [
  [
    "00A–00A2",
    "Discovery",
    "Discover GEO/SRA and EBI/BioStudies leads; merge external leads while preserving source identity and cross-references.",
  ],
  [
    "00B–00C",
    "Catalog & eligibility",
    "Build metadata catalog and apply broad scientific eligibility rules.",
  ],
  [
    "00D / 00D-EBI",
    "Metadata enrichment",
    "Resolve samples, SDRF records, phenotype-relevant fields and study metadata.",
  ],
  [
    "00D1–00D3",
    "Phenotype & relationship audit",
    "Audit resistant/sensitive labels, curated mappings, HER2 context, duplicated/umbrella study relationships.",
  ],
  [
    "00C4",
    "Independence-aware scoring",
    "Assign PRIMARY_VALIDATION, SUPPORTING, EXPLORATORY and exclusion roles while preventing duplicate evidence counting.",
  ],
  [
    "00W",
    "Scientific automation gate",
    "Block progression when mandatory scientific conditions are not satisfied.",
  ],
  [
    "00E–00G",
    "Expression & DE",
    "Acquire expression matrices, audit design, and calculate primary-validation differential expression.",
  ],
  [
    "00H–00J",
    "Cross-study consensus",
    "Measure direction concordance and derive externally replicated resistance signatures.",
  ],
  [
    "00K–00R",
    "Pathway/TGF-beta validation",
    "Validate Hallmark/pathway behavior and gene/module-level TGF-beta remodeling across studies.",
  ],
  [
    "00S",
    "Validated evidence export",
    "Export the resistance evidence that is allowed to influence drug prioritization.",
  ],
  [
    "04K–04M",
    "CMap & identity",
    "Parse CMap outputs, rank perturbagens, and resolve compound identity.",
  ],
  [
    "04N–04Q",
    "Regulatory/safety/target/network",
    "Layer regulatory context, safety flags, target annotation, and network plausibility.",
  ],
  [
    "04R",
    "Final prioritization",
    "Rank candidates using transparent evidence layers.",
  ],
  [
    "04S",
    "Target-supported docking",
    "Dock only candidates with biologically justified targets and suitable experimental structures.",
  ],
  [
    "04T–04U",
    "ADMET & integrated evidence",
    "Combine structural/ADMET context with the full evidence matrix for experimental prioritization.",
  ],
];

const guardrails = [
  "Negative CMap tau/connectivity means transcriptional opposition. It is not proof that a compound reverses trastuzumab resistance in cells, animals, or patients.",
  "A PubChem record or ClinicalTrials.gov record does not establish FDA approval.",
  "Docking is a structural hypothesis and quality-control layer. It is not proof of binding or efficacy.",
  "PAINS flags indicate possible assay interference; they are not equivalent to toxicity.",
  "Network support means biological plausibility, not causal proof.",
  "Absence of regulatory evidence should remain missing/unknown evidence rather than being converted into a negative score by default.",
  "TCGA/GDC patient data are translational context and should not be mislabeled as a direct acquired-resistance experiment.",
  "Umbrella datasets and repository mirrors must not be counted as independent biological replication.",
];

function Section({
  id,
  title,
  icon: Icon,
  children,
}: {
  id: string;
  title: string;
  icon: typeof BookOpen;
  children: ReactNode;
}) {
  return (
    <section id={id} className="panel scroll-mt-24 p-5 md:p-6">
      <div className="flex items-center gap-3">
        <div className="rounded-md bg-secondary/10 p-2 text-secondary">
          <Icon className="h-5 w-5" />
        </div>
        <h2 className="text-xl font-semibold text-primary dark:text-white">
          {title}
        </h2>
      </div>
      <div className="mt-5 text-sm leading-7 text-text-muted dark:text-slate-300">
        {children}
      </div>
    </section>
  );
}

export function DocumentationPage() {
  return (
    <div className="mx-auto max-w-[1500px] space-y-6">
      <div className="grid gap-6 xl:grid-cols-[280px_minmax(0,1fr)]">
        <aside className="panel h-fit p-4 xl:sticky xl:top-20">
          <p className="label-caps text-secondary">Documentation</p>
          <h1 className="mt-2 text-xl font-semibold">
            ATLAS technical & research guide
          </h1>
          <nav className="mt-4 space-y-1 text-sm">
            {[
              ["overview", "1. System overview"],
              ["architecture", "2. Architecture"],
              ["sources", "3. Data sources"],
              ["pipeline", "4. Pipeline stages"],
              ["statistics", "5. Statistics & charts"],
              ["settings-doc", "6. Settings"],
              ["guardrails", "7. Scientific guardrails"],
              ["operations", "8. Operations"],
              ["api-doc", "9. API reference"],
              ["development", "10. Development guide"],
              ["troubleshooting", "11. Troubleshooting"],
              ["handoff", "12. Experimental handoff"],
            ].map(([id, label]) => (
              <a
                key={id}
                href={`#${id}`}
                className="block rounded px-2 py-1.5 hover:bg-surface-container dark:hover:bg-slate-800"
              >
                {label}
              </a>
            ))}
          </nav>
        </aside>

        <div className="space-y-6">
          <div>
            <p className="label-caps text-secondary">Living documentation</p>
            <h1 className="mt-1 text-3xl font-semibold tracking-tight text-primary dark:text-white">
              ATLAS dashboard, pipeline, evidence and operations
            </h1>
            <p className="mt-3 max-w-4xl text-sm leading-7 text-text-muted dark:text-slate-400">
              This documentation is intentionally built into the frontend so the
              research interface explains what every data layer means, how it is
              generated, and what conclusions it does and does not support.
            </p>
          </div>

          <Section id="overview" title="1. System overview" icon={BookOpen}>
            <p>
              <strong className="text-primary dark:text-white">ATLAS</strong> is
              a computational research workflow for studying trastuzumab
              resistance in HER2-positive breast cancer. The dashboard is a
              read-focused interface over the pipeline rather than a replacement
              for the scientific workflow.
            </p>
            <p className="mt-3">
              Its main responsibilities are: discover and qualify external
              datasets; construct and validate a resistance signature;
              characterize pathways and reproducible modules; use perturbational
              data to nominate compounds; add identity, safety, target, network
              and structural evidence; and export a transparent shortlist for
              experimental testing.
            </p>
          </Section>

          <Section
            id="architecture"
            title="2. Software & data architecture"
            icon={Server}
          >
            <div className="grid gap-4 lg:grid-cols-2">
              <div>
                <h3 className="font-semibold text-primary dark:text-white">
                  Frontend
                </h3>
                <p className="mt-1">
                  React + TypeScript + Vite. React Router separates research
                  domains into pages, TanStack Query handles server
                  state/refetching, Recharts renders scientific visualizations,
                  and Tailwind provides the theme.
                </p>
              </div>
              <div>
                <h3 className="font-semibold text-primary dark:text-white">
                  Backend adapter
                </h3>
                <p className="mt-1">
                  FastAPI provides a stable `/api/*` boundary over ATLAS output
                  files. This lets the frontend remain unchanged if CSV readers
                  are later replaced by DuckDB queries, PostgreSQL, object
                  storage, or a remote computation service.
                </p>
              </div>
              <div>
                <h3 className="font-semibold text-primary dark:text-white">
                  Research storage
                </h3>
                <p className="mt-1">
                  Large matrices belong in Parquet/Zarr/HDF5/GCTX-style storage;
                  DuckDB is appropriate for metadata/catalog queries. Raw
                  acquisitions and checksums should remain immutable where
                  possible.
                </p>
              </div>
              <div>
                <h3 className="font-semibold text-primary dark:text-white">
                  Execution boundary
                </h3>
                <p className="mt-1">
                  The dashboard currently exposes read and configuration
                  endpoints. Pipeline execution should remain behind an explicit
                  job-control layer rather than arbitrary shell execution from
                  the browser.
                </p>
              </div>
            </div>
            <pre className="mt-5 overflow-auto rounded-md bg-primary p-4 text-xs leading-6 text-slate-100">{`Browser (React)\n      ↓ /api\nFastAPI adapter\n      ↓\nATLAS_ROOT\n ├─ data/ metadata + matrices\n ├─ results/ DE, validation, CMap, docking, integrated evidence\n ├─ scripts/ pipeline stages\n └─ results/pipeline_state/ queue and automation state`}</pre>
          </Section>

          <Section
            id="sources"
            title="3. Data sources & evidence roles"
            icon={Database}
          >
            <div className="overflow-x-auto">
              <table className="w-full min-w-[760px] text-left">
                <thead>
                  <tr className="border-b border-outline-variant dark:border-slate-700">
                    <th className="p-2 label-caps">Source</th>
                    <th className="p-2 label-caps">Role</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-outline-variant dark:divide-slate-800">
                  {sources.map(([s, r]) => (
                    <tr key={s}>
                      <td className="p-2 font-mono text-xs text-primary dark:text-white">
                        {s}
                      </td>
                      <td className="p-2">{r}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            <p className="mt-4">
              <strong>Independence rule:</strong> cross-repository mirrors, the
              same biological samples deposited in multiple repositories, and
              umbrella studies containing already-counted sub-studies should be
              linked as relationships rather than treated as independent
              validation.
            </p>
          </Section>

          <Section id="pipeline" title="4. Pipeline stages" icon={Workflow}>
            <div className="overflow-x-auto">
              <table className="w-full min-w-[850px] text-left">
                <thead>
                  <tr className="border-b border-outline-variant dark:border-slate-700">
                    <th className="p-2 label-caps">Stage</th>
                    <th className="p-2 label-caps">Layer</th>
                    <th className="p-2 label-caps">Purpose</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-outline-variant dark:divide-slate-800">
                  {stages.map(([s, l, p]) => (
                    <tr key={s}>
                      <td className="p-2 font-mono text-xs">{s}</td>
                      <td className="p-2 font-medium text-primary dark:text-white">
                        {l}
                      </td>
                      <td className="p-2">{p}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </Section>

          <Section
            id="statistics"
            title="5. Statistics & research charts"
            icon={FlaskConical}
          >
            <p>
              The{" "}
              <strong className="text-primary dark:text-white">
                Statistics
              </strong>{" "}
              page is generated from current output files under `ATLAS_ROOT`.
              Missing data are displayed as unavailable instead of being
              replaced with demo values.
            </p>
            <div className="mt-4 space-y-3">
              <div>
                <strong className="text-primary dark:text-white">
                  Volcano plot:
                </strong>{" "}
                x-axis is log2 fold change and y-axis is −log10 adjusted
                p-value/FDR. The dashboard uses the complete file for counts but
                can limit rendered points for browser performance.
              </div>
              <div>
                <strong className="text-primary dark:text-white">
                  Dataset evidence-role chart:
                </strong>{" "}
                visualizes the current independence-aware classification, such
                as PRIMARY_VALIDATION, SUPPORTING, EXPLORATORY and exclusions
                when present.
              </div>
              <div>
                <strong className="text-primary dark:text-white">
                  Pathway chart:
                </strong>{" "}
                attempts to find a compatible GSEA/pathway result containing
                pathway names and NES. Positive and negative NES should be read
                as direction of enrichment within the analyzed ranking, not as a
                universal pathway state.
              </div>
              <div>
                <strong className="text-primary dark:text-white">
                  Candidate chart:
                </strong>{" "}
                reads integrated candidate scores if a compatible numeric score
                exists. Scores are ranking aids built from evidence layers, not
                calibrated probabilities.
              </div>
              <div>
                <strong className="text-primary dark:text-white">
                  TGF-beta validation chart:
                </strong>{" "}
                reads compatible ranked-validation outputs. Cross-study sign
                differences should be described as reproducible pathway
                remodeling when direction is not monotonic across studies.
              </div>
            </div>
          </Section>

          <Section
            id="settings-doc"
            title="6. Functional settings"
            icon={ShieldCheck}
          >
            <p>
              Settings are stored in `.atlas-dashboard-settings.json` inside the
              dashboard directory. No API keys are exposed in this file or
              returned through the settings API.
            </p>
            <div className="mt-4 space-y-2">
              <p>
                <strong>ATLAS_ROOT</strong> — filesystem location of the
                research project; must exist before it can be saved.
              </p>
              <p>
                <strong>Auto refresh</strong> — enables timed dashboard refresh.
              </p>
              <p>
                <strong>Refresh interval</strong> — controls query refresh
                frequency, bounded to avoid accidental request floods.
              </p>
              <p>
                <strong>Table row limit</strong> — UI/data-display limit used
                for large result tables and research plot display density.
              </p>
              <p>
                <strong>Developer mode</strong> — unlocks runtime, queue,
                filesystem, and API diagnostics. It does not grant shell access.
              </p>
              <p>
                <strong>Scientific guardrails</strong> — preserves
                interpretation warnings in the research UI.
              </p>
              <p>
                <strong>Dense tables</strong> — compact presentation preference
                for large datasets.
              </p>
            </div>
          </Section>

          <Section
            id="guardrails"
            title="7. Scientific interpretation guardrails"
            icon={ShieldCheck}
          >
            <div className="grid gap-3 lg:grid-cols-2">
              {guardrails.map((g) => (
                <div
                  key={g}
                  className="rounded-md border border-outline-variant bg-surface-low p-3 dark:border-slate-700 dark:bg-slate-900"
                >
                  {g}
                </div>
              ))}
            </div>
          </Section>

          <Section
            id="operations"
            title="8. Operations & runbook"
            icon={Terminal}
          >
            <p>
              Run backend and frontend in separate terminals during development:
            </p>
            <pre className="mt-3 overflow-auto rounded-md bg-primary p-4 text-xs leading-6 text-slate-100">{`# Terminal 1\ncd /home/regulus/Documents/ATLAS/ATLAS-React-Dashboard-v2\n./run_backend.sh\n\n# Terminal 2\ncd /home/regulus/Documents/ATLAS/ATLAS-React-Dashboard-v2\n./run_frontend.sh\n\n# Browser\nhttp://localhost:5173`}</pre>
            <p className="mt-4">
              The queue worker is intentionally separate from the website.
              Developer mode can read `atlas-dataset-queue.service` state and
              queue CSV status but does not start, stop or mutate the research
              pipeline.
            </p>
          </Section>

          <Section id="api-doc" title="9. API reference" icon={GitBranch}>
            <div className="grid gap-2 md:grid-cols-2">
              {[
                [
                  "GET /api/health",
                  "Backend/root connectivity and API version.",
                ],
                [
                  "GET /api/dashboard",
                  "Dashboard metrics, funnel, candidates and activity.",
                ],
                ["GET /api/datasets", "Independence-aware dataset catalog."],
                ["GET /api/signature", "Differential-expression rows."],
                ["GET /api/cmap", "CMap/integrated evidence table."],
                [
                  "GET /api/docking",
                  "Auto-detected docking/structural output.",
                ],
                ["GET /api/candidates", "Integrated candidate evidence."],
                [
                  "GET /api/research/statistics",
                  "Derived research-chart payloads.",
                ],
                ["GET /api/settings", "Non-secret dashboard settings."],
                [
                  "PUT /api/settings",
                  "Validate and persist dashboard settings.",
                ],
                [
                  "GET /api/developer/statistics",
                  "Developer diagnostics; requires developer mode.",
                ],
              ].map(([e, d]) => (
                <div
                  key={e}
                  className="rounded-md border border-outline-variant p-3 dark:border-slate-700"
                >
                  <code className="text-xs text-secondary">{e}</code>
                  <p className="mt-1">{d}</p>
                </div>
              ))}
            </div>
          </Section>

          <Section
            id="development"
            title="10. Development guide"
            icon={GitBranch}
          >
            <p>
              <strong>Adding a research page:</strong> create a page component
              under `frontend/src/pages`, add reusable components under
              `components`, expose a backend endpoint if new data are required,
              register the route in `App.tsx`, then add the navigation entry in
              `Sidebar.tsx`.
            </p>
            <p className="mt-3">
              <strong>Adding a chart:</strong> calculate statistical values on
              the backend from full datasets where practical, return only the
              fields necessary for visualization, and document threshold
              definitions. Avoid deriving scientific meaning from visual
              formatting alone.
            </p>
            <p className="mt-3">
              <strong>Scaling:</strong> when outputs grow beyond convenient CSV
              reads, keep the API contract but implement readers against
              DuckDB/Parquet or a database. This avoids coupling React
              components to storage internals.
            </p>
          </Section>

          <Section
            id="troubleshooting"
            title="11. Troubleshooting"
            icon={Terminal}
          >
            <div className="space-y-3">
              <p>
                <strong>Backend says root missing:</strong> open Settings and
                verify `ATLAS_ROOT`, or check `/api/health`.
              </p>
              <p>
                <strong>Chart says unavailable:</strong> the reader did not find
                a compatible current output. This is preferable to plotting
                fabricated values. Check the source path shown on the
                corresponding results page.
              </p>
              <p>
                <strong>Developer endpoint returns 403:</strong> enable
                Developer mode in Settings and save.
              </p>
              <p>
                <strong>Frontend cannot reach API:</strong> confirm Uvicorn is
                on `127.0.0.1:8000`; Vite proxies `/api` there during
                development.
              </p>
              <p>
                <strong>Settings do not save:</strong> ensure the dashboard
                directory is writable and the requested `ATLAS_ROOT` exists.
              </p>
              <p>
                <strong>Port already in use:</strong> stop the previous dev
                server or change the port consistently in Uvicorn and
                `vite.config.ts`.
              </p>
            </div>
          </Section>

          <Section
            id="handoff"
            title="12. Experimental handoff"
            icon={FlaskConical}
          >
            <p>
              The computational endpoint is a ranked set of candidates and
              mechanisms for laboratory validation. A minimal experimental
              design should compare a resistant HER2-positive model under
              vehicle/control, trastuzumab, candidate alone, and candidate +
              trastuzumab, with biological replication and an assay appropriate
              to the biological claim.
            </p>
            <p className="mt-3">
              Where a candidate is linked to a target mechanistically,
              orthogonal validation such as target knockdown/perturbation can
              help distinguish target-mediated effects from nonspecific drug
              effects. The website should continue to display computational
              evidence separately from future wet-lab evidence rather than
              merging them into one opaque score.
            </p>
          </Section>
        </div>
      </div>
    </div>
  );
}
