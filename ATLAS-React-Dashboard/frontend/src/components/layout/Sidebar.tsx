import { Activity, BarChart3, Database, Dna, FlaskConical, LayoutDashboard, Settings, X } from "lucide-react";
import { NavLink } from "react-router-dom";

const items = [
  { to: "/", label: "Dashboard", icon: LayoutDashboard },
  { to: "/datasets", label: "Datasets", icon: Database },
  { to: "/signature", label: "Signature Discovery", icon: Dna },
  { to: "/cmap", label: "CMap Results", icon: BarChart3 },
  { to: "/docking", label: "Docking Results", icon: FlaskConical },
  { to: "/candidates", label: "Final Candidates", icon: Activity }
];

export function Sidebar({ open, onClose }: { open: boolean; onClose: () => void }) {
  return (
    <>
      {open && <button className="fixed inset-0 z-30 bg-black/20 lg:hidden" aria-label="Close navigation overlay" onClick={onClose} />}
      <aside className={["fixed left-0 top-0 z-40 flex h-screen w-64 flex-col border-r border-outline-variant bg-surface px-4 py-8 transition-transform dark:border-slate-700 dark:bg-slate-950", open ? "translate-x-0" : "-translate-x-full lg:translate-x-0"].join(" ")}>
        <div className="mb-8 flex items-start justify-between">
          <div>
            <h1 className="text-2xl font-bold tracking-tight text-primary dark:text-white">ATLAS</h1>
            <p className="mt-2 text-xs leading-5 text-text-muted dark:text-slate-400">Trastuzumab resistance computational pipeline</p>
          </div>
          <button className="lg:hidden" onClick={onClose} aria-label="Close navigation"><X className="h-5 w-5" /></button>
        </div>
        <nav className="flex-1 space-y-1">
          {items.map(({ to, label, icon: Icon }) => (
            <NavLink key={to} to={to} onClick={onClose} className={({ isActive }) => ["flex items-center gap-3 rounded-md px-3 py-2 text-sm font-medium transition-colors", isActive ? "border-r-2 border-secondary bg-secondary/10 text-secondary dark:text-emerald-300" : "text-text-muted hover:bg-surface-container dark:text-slate-300 dark:hover:bg-slate-900"].join(" ")}>
              <Icon className="h-4 w-4" />{label}
            </NavLink>
          ))}
        </nav>
        <NavLink to="/settings" onClick={onClose} className="flex items-center gap-3 border-t border-outline-variant px-3 pt-4 text-sm text-text-muted dark:border-slate-700 dark:text-slate-300"><Settings className="h-4 w-4" />Settings</NavLink>
      </aside>
    </>
  );
}
