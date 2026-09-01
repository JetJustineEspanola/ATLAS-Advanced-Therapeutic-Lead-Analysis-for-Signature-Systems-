import type { Config } from "tailwindcss";

export default {
  content: ["./index.html", "./src/**/*.{ts,tsx}"],
  darkMode: "class",
  theme: {
    extend: {
      colors: {
        background: "#f7f9fb",
        surface: "#f7f9fb",
        "surface-lowest": "#ffffff",
        "surface-low": "#f2f4f6",
        "surface-container": "#eceef0",
        "surface-high": "#e6e8ea",
        primary: "#091426",
        "primary-container": "#1e293b",
        secondary: "#006a61",
        "secondary-container": "#86f2e4",
        outline: "#75777d",
        "outline-variant": "#c5c6cd",
        "text-main": "#191c1e",
        "text-muted": "#45474c",
        warning: "#f59e0b",
        danger: "#ba1a1a"
      },
      fontFamily: {
        sans: ["Inter", "sans-serif"],
        mono: ["JetBrains Mono", "monospace"]
      },
      boxShadow: {
        panel: "0 1px 2px rgba(9,20,38,.05)"
      }
    }
  },
  plugins: []
} satisfies Config;
