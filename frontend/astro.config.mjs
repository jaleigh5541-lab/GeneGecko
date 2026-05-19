import { defineConfig } from "astro/config";
import react from "@astrojs/react";
import tailwind from "@astrojs/tailwind";

export default defineConfig({
  integrations: [react(), tailwind({ applyBaseStyles: false })],
  server: { port: 5173 },
  vite: {
    server: {
      proxy: {
        "/api": "http://localhost:5000",
      },
      watch: {
        usePolling: true,
      },
    },
  },
});
