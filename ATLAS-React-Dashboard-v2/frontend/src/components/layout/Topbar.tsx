import { Menu, Moon, RefreshCw, Sun } from "lucide-react";

export function Topbar({ title, onMenu, dark, onToggleDark, onRefresh }: { title: string; onMenu: () => void; dark: boolean; onToggleDark: () => void; onRefresh: () => void }) {
  return (
    <header className="sticky top-0 z-20 flex h-16 items-center justify-between border-b border-outline-variant bg-surface/95 px-4 backdrop-blur md:px-6 dark:border-slate-700 dark:bg-slate-950/95">
      <div className="flex items-center gap-3">
        <button className="lg:hidden" onClick={onMenu} aria-label="Open navigation"><Menu className="h-5 w-5" /></button>
        <h2 className="text-xl font-semibold text-primary dark:text-white">{title}</h2>
      </div>
      <div className="flex items-center gap-2">
        <button onClick={onToggleDark} className="rounded-md border border-outline-variant p-2 text-text-muted hover:bg-surface-container dark:border-slate-700 dark:text-slate-300 dark:hover:bg-slate-900" aria-label="Toggle theme">{dark ? <Sun className="h-4 w-4" /> : <Moon className="h-4 w-4" />}</button>
        <button onClick={onRefresh} className="flex items-center gap-2 rounded-md bg-secondary px-3 py-2 text-xs font-bold uppercase tracking-wider text-white hover:bg-teal-800"><RefreshCw className="h-4 w-4" />Refresh</button>
      </div>
    </header>
  );
}
