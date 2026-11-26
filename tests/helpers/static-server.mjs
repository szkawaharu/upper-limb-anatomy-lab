// Each browser audit serves its own checkout, avoiding stale or absent previews.
import { createServer } from 'node:http';
import { readFile } from 'node:fs/promises';
import { resolve, extname, sep } from 'node:path';
export async function staticServer(root = process.cwd()) {
  root = resolve(root);
  const server = createServer(async (req, res) => {
    try {
      const path = resolve(root, '.' + decodeURIComponent(new URL(req.url, 'http://local').pathname));
      if (path !== root && !path.startsWith(root + sep)) { res.writeHead(403).end(); return; }
      const file = path === root ? resolve(root, 'index.html') : path;
      const data = await readFile(file);
      res.writeHead(200, { 'content-type': ({'.html':'text/html','.js':'text/javascript','.json':'application/json','.png':'image/png'})[extname(file)] || 'application/octet-stream', 'cache-control':'no-store' }).end(data);
    } catch { res.writeHead(404).end(); }
  });
  await new Promise((resolve, reject) => { server.once('error', reject); server.listen(0, '127.0.0.1', resolve); });
  return { url: `http://127.0.0.1:${server.address().port}/index.html`, close: () => new Promise(resolve => server.close(resolve)) };
}
