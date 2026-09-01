import { useEffect, useState } from "react";
import { Outlet, useLocation } from "react-router-dom";
import { useQueryClient } from "@tanstack/react-query";
import { Sidebar } from "./Sidebar";
import { Topbar } from "./Topbar";

const titles: Record<string, string> = {
  "/": "Overview", "/datasets": "Datasets", "/signature": "Signature Discovery", "/cmap": "CMap Results", "/docking": "Docking Results", "/candidates": "Final Candidates", "/settings": "Settings"
};

export function Shell() {
  const [navOpen, setNavOpen] = useState(false);
  const [dark, setDark] = useState(() => localStorage.getItem("atlas-theme") === "dark");
  const location = useLocation();
  const queryClient = useQueryClient();

  useEffect(() => {
    document.documentElement.classList.toggle("dark", dark);
    localStorage.setItem("atlas-theme", dark ? "dark" : "light");
  }, [dark]);

  return (
    <div className="min-h-screen bg-background text-text-main dark:bg-slate-950 dark:text-slate-100">
      <Sidebar open={navOpen} onClose={() => setNavOpen(false)} />
      <div className="min-h-screen lg:ml-64">
        <Topbar title={titles[location.pathname] ?? "ATLAS"} onMenu={() => setNavOpen(true)} dark={dark} onToggleDark={() => setDark(v => !v)} onRefresh={() => queryClient.invalidateQueries()} />
        <main className="p-4 md:p-6"><Outlet /></main>
      </div>
    </div>
  );
}
