import { AlertTriangle, LoaderCircle } from "lucide-react";

export function LoadingState() {
  return <div className="panel flex min-h-48 items-center justify-center gap-3 p-8 text-text-muted dark:text-slate-300"><LoaderCircle className="h-5 w-5 animate-spin" />Loading ATLAS data…</div>;
}

export function ErrorState({ error }: { error: Error }) {
  return <div className="panel flex min-h-48 items-center justify-center gap-3 p-8 text-danger dark:text-red-300"><AlertTriangle className="h-5 w-5" /><span>{error.message}</span></div>;
}
