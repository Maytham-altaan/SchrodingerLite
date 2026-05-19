"""Maestro — molecular visualization (Schrödinger Maestro equivalent).

Two viewing modes:
  - PyMOL native GUI  (offline, full features)
  - Py3Dmol HTML page (auto-opens in your browser)

The HTML viewer:
  - Embeds 3Dmol.js *inline* (downloads once, caches under share/) so the page
    works in Chrome under file:// without needing CDN access at view time.
  - Embeds each structure as base64 → bulletproof against any file content.
  - Picks a sensible default style per model (cartoon for proteins, sticks for
    small molecules / docked poses) so a pose-only file isn't blank.
"""
from __future__ import annotations
import argparse
import base64
import platform
import subprocess
import urllib.request
import webbrowser
from pathlib import Path
from typing import Iterable

from .utils import banner, get_logger, home, projects_dir, require


# ---------------------------------------------------------------------------
# 3Dmol.js bootstrap — download once, then embed inline forever

_3DMOL_URLS = [
    "https://cdn.jsdelivr.net/npm/3dmol@2.4.0/build/3Dmol-min.js",
    "https://cdn.jsdelivr.net/npm/3dmol/build/3Dmol-min.js",
    "https://unpkg.com/3dmol@2.4.0/build/3Dmol-min.js",
    "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.6/3Dmol-min.js",
    "https://3Dmol.csb.pitt.edu/build/3Dmol-min.js",
]


def _ensure_3dmol_js() -> str | None:
    """Return 3Dmol.js source. Cached at <home>/share/3Dmol-min.js.
    Returns None if all download mirrors fail."""
    cache = home() / "share" / "3Dmol-min.js"
    cache.parent.mkdir(parents=True, exist_ok=True)

    if cache.exists() and cache.stat().st_size > 200_000:
        return cache.read_text(encoding="utf-8", errors="replace")

    for url in _3DMOL_URLS:
        try:
            print(f"  Downloading 3Dmol.js from {url.split('/')[2]} …")
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=20) as r:
                data = r.read()
            if len(data) > 200_000:
                cache.write_bytes(data)
                print(f"  ✓ Cached 3Dmol.js ({len(data)//1024} KB) → {cache}")
                return data.decode("utf-8", errors="replace")
        except Exception as e:
            print(f"    failed ({e.__class__.__name__}: {e})")
            continue
    return None


# ---------------------------------------------------------------------------
# HTML template — 3Dmol.js is injected at __THREEDMOL_JS__ (inline) with a
# CDN <script> as a graceful fallback.

HTML_TEMPLATE = """<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>Maestro Viewer — SchrodingerLite</title>
<style>
  html, body { margin: 0; padding: 0; height: 100%; font-family: -apple-system, sans-serif; background: #0b0c10; color: #e8e8e8; overflow: hidden; }
  header { padding: 12px 18px; background: #1a1d23; border-bottom: 1px solid #2a2d33; display: flex; align-items: center; justify-content: space-between; }
  header h1 { margin: 0; font-size: 16px; font-weight: 600; color: #6cc7ff; }
  header .title { font-size: 12px; color: #aaa; margin-left: 10px; }
  #viewer { width: 100vw; height: calc(100vh - 110px); position: relative; }
  #status { position: absolute; top: 10px; right: 14px; font-size: 11px; color: #6cc7ff; z-index: 100; background: rgba(0,0,0,0.6); padding: 4px 8px; border-radius: 4px; }
  #legend { position: absolute; top: 10px; left: 14px; font-size: 11px; color: #ddd; z-index: 100; background: rgba(0,0,0,0.6); padding: 6px 10px; border-radius: 4px; max-width: 320px; }
  #legend .row { display:flex; align-items:center; margin: 2px 0; }
  #legend .sw { width: 10px; height: 10px; border-radius: 50%; margin-right: 6px; display:inline-block; }
  #controls { padding: 10px 18px; background: #1a1d23; border-top: 1px solid #2a2d33; font-size: 13px; }
  #err { padding: 24px; color: #ff8a8a; font-size: 14px; line-height: 1.6; }
  #err code { background: #1a1d23; padding: 2px 6px; border-radius: 3px; }
  button { background:#2a2d33; color:#e8e8e8; border:1px solid #3a3d43; padding:5px 12px; border-radius:4px; cursor:pointer; margin-right:6px; }
  button:hover { background:#3a3d43; }
</style>
</head>
<body>
<header>
  <div><h1 style="display:inline">Maestro Viewer</h1><span class="title">__TITLE__</span></div>
  <div style="font-size:11px;color:#666">3Dmol.js · SchrodingerLite</div>
</header>
<div id="viewer">
  <div id="status">loading…</div>
  <div id="legend"></div>
</div>
<div id="controls">
  <button onclick="applyStyle('default')">Default</button>
  <button onclick="applyStyle('cartoon')">Cartoon</button>
  <button onclick="applyStyle('stick')">Sticks</button>
  <button onclick="applyStyle('sphere')">Spheres</button>
  <button onclick="applyStyle('surface')">Surface + Sticks</button>
  <button onclick="if(window.viewer){viewer.zoomTo();viewer.render();}">Reset View</button>
  <button onclick="toggleSpin()">Toggle Spin</button>
</div>

__MODEL_TAGS__

<script>__THREEDMOL_JS__</script>
<!-- CDN fallback (only used if the inlined library above is missing) -->
<script>
if (typeof $3Dmol === 'undefined') {
  var s = document.createElement('script');
  s.src = 'https://cdn.jsdelivr.net/npm/3dmol@2.4.0/build/3Dmol-min.js';
  s.onerror = function(){ window._3dmolFail = true; };
  document.head.appendChild(s);
}
</script>

<script>
function showError(msg) {
  var v = document.getElementById('viewer');
  v.innerHTML = '<div id="err">' + msg + '</div>';
}
function setStatus(msg) {
  var s = document.getElementById('status');
  if (s) s.textContent = msg;
}
function b64decodeUtf8(b64) {
  var bin = atob(b64);
  var bytes = new Uint8Array(bin.length);
  for (var i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  return new TextDecoder('utf-8').decode(bytes);
}

// Per-model color palette for ligands
var LIG_COLORS = ['cyanCarbon','magentaCarbon','yellowCarbon','greenCarbon',
                  'orangeCarbon','purpleCarbon','pinkCarbon','blueCarbon'];

var MODELS = [];   // [{name, role, color, model}]
var SPINNING = false;

function bootstrap() {
  if (window._3dmolFail || typeof $3Dmol === 'undefined') {
    // Wait a beat for the CDN fallback, then give up
    setTimeout(function(){
      if (typeof $3Dmol === 'undefined') {
        showError(
          '<b>3Dmol library failed to load.</b><br><br>' +
          'The inlined library was missing and the CDN fallback also failed.<br>' +
          'Re-run the viewer with internet access, or load the page in a browser ' +
          'that allows CDN scripts under <code>file://</code> (Chrome / Firefox).'
        );
      } else {
        renderAll();
      }
    }, 1500);
    return;
  }
  renderAll();
}

function renderAll() {
  try {
    setStatus('creating viewer…');
    var element = document.getElementById('viewer');
    window.viewer = $3Dmol.createViewer(element, {backgroundColor: '#0b0c10'});

    var tags = document.querySelectorAll('script[type="chemical/x-base64"]');
    if (tags.length === 0) { showError('No model data found in HTML.'); return; }

    setStatus('loading ' + tags.length + ' model(s)…');
    var ligIdx = 0;
    for (var i = 0; i < tags.length; i++) {
      var t = tags[i];
      var fmt  = t.getAttribute('data-format') || 'pdb';
      var name = t.getAttribute('data-name') || ('model ' + (i+1));
      var role = t.getAttribute('data-role') || 'auto';
      var text = b64decodeUtf8(t.textContent.trim());
      var m = viewer.addModel(text, fmt);

      // Auto-detect role from atom count if not specified
      if (role === 'auto') {
        var n = m.selectedAtoms({}).length;
        role = (n > 200) ? 'receptor' : 'ligand';
      }
      var color = null;
      if (role === 'ligand') {
        color = LIG_COLORS[ligIdx % LIG_COLORS.length];
        ligIdx++;
      }
      MODELS.push({name: name, role: role, color: color, model: m, modelIdx: i});
    }

    buildLegend();
    setStatus('rendering…');
    applyStyle('default');
    setTimeout(function(){ setStatus(''); document.getElementById('status').style.display='none'; }, 800);
  } catch (e) {
    showError('<b>Render error:</b> ' + (e.message || e) + '<br><br>' +
              '<pre style="white-space:pre-wrap;font-size:11px">' + (e.stack || '') + '</pre>');
  }
}

function buildLegend() {
  var el = document.getElementById('legend');
  if (!el || MODELS.length === 0) return;
  var html = '<div style="font-weight:600;color:#6cc7ff;margin-bottom:4px">Models</div>';
  for (var i = 0; i < MODELS.length; i++) {
    var m = MODELS[i];
    var label = m.name + '  (' + m.role + ')';
    var swColor = m.role === 'receptor' ? '#6cc7ff'
                : m.color ? m.color.replace('Carbon','').toLowerCase()
                : '#aaa';
    html += '<div class="row"><span class="sw" style="background:' + swColor + '"></span>' + label + '</div>';
  }
  el.innerHTML = html;
}

function applyStyle(style) {
  if (typeof viewer === 'undefined') return;
  viewer.setStyle({}, {});  // clear
  viewer.removeAllSurfaces();

  for (var i = 0; i < MODELS.length; i++) {
    var m = MODELS[i];
    var sel = {model: m.modelIdx};
    var s;
    if (style === 'default') {
      s = (m.role === 'receptor')
        ? {cartoon: {color: 'spectrum'}, stick: {hidden: true}}
        : {stick: {colorscheme: m.color || 'orangeCarbon', radius: 0.18}};
    } else if (style === 'cartoon') {
      s = (m.role === 'receptor') ? {cartoon: {color: 'spectrum'}}
                                   : {stick: {colorscheme: m.color || 'orangeCarbon'}};
    } else if (style === 'stick') {
      s = {stick: {colorscheme: m.role === 'receptor' ? 'whiteCarbon' : (m.color || 'orangeCarbon')}};
    } else if (style === 'sphere') {
      s = {sphere: {scale: 0.35, colorscheme: m.color || 'default'}};
    } else if (style === 'surface') {
      if (m.role === 'receptor') {
        s = {cartoon: {color: 'spectrum'}};
      } else {
        s = {stick: {colorscheme: m.color || 'orangeCarbon'}};
      }
    }
    viewer.setStyle(sel, s);
  }

  if (style === 'surface') {
    // Add a transparent surface only on receptor models
    for (var j = 0; j < MODELS.length; j++) {
      if (MODELS[j].role === 'receptor') {
        viewer.addSurface($3Dmol.SurfaceType.VDW,
          {opacity: 0.55, color: 'white'},
          {model: MODELS[j].modelIdx});
      }
    }
  }

  viewer.zoomTo();
  viewer.render();
}

function toggleSpin() {
  if (typeof viewer === 'undefined') return;
  SPINNING = !SPINNING;
  viewer.spin(SPINNING ? 'y' : false);
  viewer.render();
}

if (document.readyState === 'loading') {
  window.addEventListener('load', bootstrap);
} else {
  bootstrap();
}
</script>
</body></html>
"""


# ---------------------------------------------------------------------------

_EXT_TO_FMT = {"pdb": "pdb", "ent": "pdb", "sdf": "sdf", "mol": "sdf",
               "mol2": "mol2", "xyz": "xyz", "pdbqt": "pdbqt"}


def _model_tag(idx: int, path: Path, role: str = "auto") -> str:
    """Encode a structure file into a base64 script tag — bulletproof against
    template-literal-breaking content."""
    ext = path.suffix.lstrip(".").lower()
    fmt = _EXT_TO_FMT.get(ext, "pdb")
    raw = path.read_bytes()
    b64 = base64.b64encode(raw).decode("ascii")
    name = path.stem
    return (f'<script type="chemical/x-base64" '
            f'data-format="{fmt}" data-name="{name}" data-role="{role}" '
            f'id="model{idx}">{b64}</script>')


def html_view(files: Iterable[Path | tuple[Path, str]],
              out_html: Path | None = None,
              title: str = "",
              open_browser: bool = True) -> Path:
    """Render a stand-alone HTML viewer for one or more structures.

    Each entry in *files* is either a Path (role autodetected) or
    a (Path, role) tuple where role is "receptor", "ligand", or "auto".
    """
    entries: list[tuple[Path, str]] = []
    for f in files:
        if isinstance(f, tuple):
            entries.append((Path(f[0]), f[1]))
        else:
            entries.append((Path(f), "auto"))

    if out_html is None:
        out_html = projects_dir() / "maestro_view.html"
    out_html.parent.mkdir(parents=True, exist_ok=True)

    model_tags = "\n".join(_model_tag(i, p, role)
                           for i, (p, role) in enumerate(entries))

    if not title:
        title = ", ".join(p.name for p, _ in entries[:3])
        if len(entries) > 3:
            title += f" (+{len(entries)-3} more)"

    # Inline 3Dmol.js (downloads + caches on first run)
    js = _ensure_3dmol_js() or ""
    if not js:
        print("  ! 3Dmol.js could not be downloaded — the page will try a CDN at view time.")

    html = (HTML_TEMPLATE
            .replace("__TITLE__", title)
            .replace("__MODEL_TAGS__", model_tags)
            .replace("__THREEDMOL_JS__", js))
    out_html.write_text(html, encoding="utf-8")

    if open_browser:
        _open_in_chrome_preferred(out_html)
    return out_html


def _open_in_chrome_preferred(path: Path) -> None:
    """Open in Chrome / Firefox first (Safari blocks CDN under file://)."""
    opened = False
    if platform.system() == "Darwin":
        for app in ("Google Chrome", "Firefox", "Brave Browser",
                    "Microsoft Edge", "Arc"):
            try:
                subprocess.run(["open", "-a", app, str(path.resolve())],
                               check=True, capture_output=True)
                print(f"  Opened in {app}")
                opened = True
                break
            except subprocess.CalledProcessError:
                continue
    if not opened:
        webbrowser.open(f"file://{path.resolve()}")
        print("  Opened in default browser.")
        print("  If 3D is blank, try Chrome:")
        print(f"    open -a 'Google Chrome' \"{path.resolve()}\"")


def pymol_view(files: Iterable[Path]) -> None:
    """Open files in PyMOL (open-source)."""
    require("pymol")
    cmd = ["pymol"] + [str(f) for f in files]
    subprocess.Popen(cmd)


# ---------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        prog="maestro",
        description="Maestro — view proteins, ligands, and docking poses."
    )
    p.add_argument("files", nargs="+", type=Path,
                   help="One or more structure files (.pdb .sdf .mol2 .pdbqt)")
    p.add_argument("--engine", choices=["html", "pymol"], default="html",
                   help="Viewer engine (default: html, opens in browser)")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="HTML output path (html mode only)")
    p.add_argument("--no-open", action="store_true",
                   help="Write file but don't auto-open browser")
    p.add_argument("--title", default="", help="Window title")
    args = p.parse_args(argv)

    banner("Maestro — molecular viewer",
           f"engine={args.engine}  files={len(args.files)}")

    if args.engine == "html":
        out = html_view(args.files, out_html=args.output,
                        title=args.title, open_browser=not args.no_open)
        print(f"✓ Viewer written: {out}")
    else:
        pymol_view(args.files)
        print("✓ Launched PyMOL")


if __name__ == "__main__":
    main()
