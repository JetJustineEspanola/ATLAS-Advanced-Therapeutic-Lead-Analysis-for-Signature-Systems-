import { createBrowserRouter, RouterProvider } from "react-router-dom";
import { Shell } from "./components/layout/Shell";
import { DashboardPage } from "./pages/DashboardPage";
import { DatasetsPage } from "./pages/DatasetsPage";
import { GenericTablePage } from "./pages/GenericTablePage";
import { SettingsPage } from "./pages/SettingsPage";
import { ResearchStatisticsPage } from "./pages/ResearchStatisticsPage";
import { DeveloperPage } from "./pages/DeveloperPage";
import { DocumentationPage } from "./pages/DocumentationPage";
import { api } from "./lib/api";

const router = createBrowserRouter([{ element: <Shell />, children: [
  { path: "/", element: <DashboardPage /> },
  { path: "/statistics", element: <ResearchStatisticsPage /> },
  { path: "/datasets", element: <DatasetsPage /> },
  { path: "/signature", element: <GenericTablePage queryKey="signature" title="Signature Discovery" subtitle="Differential-expression evidence from the current discovery analysis." loader={api.signature} /> },
  { path: "/cmap", element: <GenericTablePage queryKey="cmap" title="CMap Results" subtitle="Connectivity-map evidence exposed from current ATLAS outputs. Negative connectivity is transcriptional opposition, not efficacy proof." loader={api.cmap} /> },
  { path: "/docking", element: <GenericTablePage queryKey="docking" title="Docking Results" subtitle="Target-supported structural evidence. Docking is not efficacy or binding proof." loader={api.docking} /> },
  { path: "/candidates", element: <GenericTablePage queryKey="candidates" title="Final Candidates" subtitle="Integrated computational evidence used for experimental prioritization." loader={api.candidates} /> },
  { path: "/documentation", element: <DocumentationPage /> },
  { path: "/developer", element: <DeveloperPage /> },
  { path: "/settings", element: <SettingsPage /> }
]}]);

export default function App() { return <RouterProvider router={router} />; }
