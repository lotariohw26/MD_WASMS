// Live-reload script taken from https://github.com/evanw/esbuild/issues/802#issuecomment-819578182
import esbuild from "esbuild";
import http from "http";
import { spawn } from "child_process";
import buildExamples from "../examples/build_examples_json.mjs";
import process from "process";

const OUT_DIR = "./shinylive";

const clients = [];

let watch = false;
let serve = false;
if (process.argv.some((x) => x === "--watch")) {
  watch = true;
}
if (process.argv.some((x) => x === "--serve")) {
  watch = true;
  serve = true;
}

const onRebuild = (error, result) => {
  clients.forEach((res) => res.write("data: update\n\n"));
  clients.length = 0;

  console.log(
    `[${new Date().toISOString()}]` +
      (error ? error : " Rebuilding JS files...")
  );
};

let watchProp = {};
if (watch) {
  watchProp = { watch: { onRebuild } };
}

esbuild
  .build({
    bundle: true,
    entryPoints: ["src/Components/App.tsx"],
    outfile: `${OUT_DIR}/shinylive.js`,
    format: "esm",
    target: "es2020",
    ...watchProp,
    plugins: [
      {
        name: "example-builder",
        setup(build) {
          build.onStart(() => {
            // On every rebuild make sure the examples are up to date.
            // One issue is this won't force esbuild to watch for the changes
            // of the example files themselves so live-reloading won't work
            buildExamples();
          });
        },
      },
    ],
  })
  .catch(() => process.exit(1));

esbuild
  .build({
    bundle: true,
    entryPoints: ["src/pyodide-worker.ts", "src/inject-socket.ts"],
    outdir: `${OUT_DIR}`,
    format: "esm",
    target: "es2020",
    ...watchProp,
  })
  .catch(() => process.exit(1));

esbuild
  .build({
    bundle: false,
    entryPoints: ["src/run-python-blocks.ts"],
    outdir: `${OUT_DIR}`,
    format: "esm",
    target: "es2020",
    ...watchProp,
  })
  .catch(() => process.exit(1));

esbuild
  .build({
    bundle: true,
    entryPoints: ["src/serviceworker.ts"],
    outdir: "./",
    format: "iife",
    target: "es2020",
    ...watchProp,
  })
  .catch(() => process.exit(1));

if (serve) {
  esbuild.serve({ servedir: "site/", port: 3001 }, {}).then(() => {
    http
      .createServer((req, res) => {
        const { url, method, headers } = req;

        if (req.url === "/esbuild")
          return clients.push(
            res.writeHead(200, {
              "Content-Type": "text/event-stream",
              "Cache-Control": "no-cache",
              Connection: "keep-alive",
            })
          );

        req.pipe(
          http.request(
            { hostname: "0.0.0.0", port: 3001, path: url, method, headers },
            (proxyRes) => {
              if (url === "/shinylive/shinylive.js") {
                // JS code for does auto-reloading. We'll inject it into
                // shinylive.js as it's sent.
                const jsReloadCode = `(() => {
                  if (window.location.host.includes("localhost")) {
                    console.log('%c~~~~~ Live Reload Enabled ~~~~~~', 'font-weight:bold;font-size:20px;color:white;display:block;background-color:green;padding:4px;border-radius:5px;');
                    new EventSource("/esbuild").onmessage = () => location.reload();
                  }
                })();`;

                const newHeaders = {
                  ...proxyRes.headers,
                  "content-length":
                    parseInt(proxyRes.headers["content-length"], 10) +
                    jsReloadCode.length,
                };

                res.writeHead(proxyRes.statusCode, newHeaders);
                res.write(jsReloadCode);
              } else {
                res.writeHead(proxyRes.statusCode, proxyRes.headers);
              }
              proxyRes.pipe(res, { end: true });
            }
          ),
          { end: true }
        );
      })
      .listen(3000);

    setTimeout(() => {
      const op = {
        darwin: ["open"],
        linux: ["xdg-open"],
        win32: ["cmd", "/c", "start"],
      };
      const platform = process.platform;
      if (clients.length === 0)
        spawn(op[platform][0], [`http://localhost:3000/examples`]);
    }, 1000); //open the default browser only if it is not opened yet
  });
}
