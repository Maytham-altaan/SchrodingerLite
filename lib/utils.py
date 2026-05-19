"""Shared utilities: paths, project dirs, logging, file helpers."""
from __future__ import annotations
import os
import sys
import time
import shutil
import logging
from pathlib import Path
from datetime import datetime

try:
    from rich.console import Console
    from rich.logging import RichHandler
    _HAS_RICH = True
except Exception:
    _HAS_RICH = False

# ---- Paths -----------------------------------------------------------------

def home() -> Path:
    """Root of the SchrodingerLite install (AI DOCKING/SchrodingerLite)."""
    env = os.environ.get("SCHRODINGER_LITE_HOME")
    if env:
        return Path(env)
    # Fall back to walking up from this file
    return Path(__file__).resolve().parent.parent

def projects_dir() -> Path:
    """All user docking projects live here."""
    env = os.environ.get("SCHRODINGER_LITE_PROJECTS")
    if env:
        p = Path(env)
    else:
        p = home() / "projects"
    p.mkdir(parents=True, exist_ok=True)
    return p

def new_project(name: str | None = None) -> Path:
    """Create a timestamped project folder inside projects/."""
    stamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    folder = f"{stamp}_{name}" if name else stamp
    p = projects_dir() / folder
    (p / "input").mkdir(parents=True, exist_ok=True)
    (p / "output").mkdir(parents=True, exist_ok=True)
    (p / "logs").mkdir(parents=True, exist_ok=True)
    return p

# ---- Console / logging -----------------------------------------------------

def console():
    if _HAS_RICH:
        return Console()
    return None

def banner(title: str, subtitle: str = ""):
    bar = "=" * 70
    if _HAS_RICH:
        c = Console()
        c.print(f"[bold cyan]{bar}[/]")
        c.print(f"[bold white]  {title}[/]")
        if subtitle:
            c.print(f"[dim]  {subtitle}[/]")
        c.print(f"[bold cyan]{bar}[/]")
    else:
        print(bar)
        print(f"  {title}")
        if subtitle:
            print(f"  {subtitle}")
        print(bar)

def get_logger(name: str, logfile: Path | None = None) -> logging.Logger:
    log = logging.getLogger(name)
    if log.handlers:
        return log
    log.setLevel(logging.INFO)
    if _HAS_RICH:
        log.addHandler(RichHandler(rich_tracebacks=True, show_time=True, show_path=False))
    else:
        h = logging.StreamHandler()
        h.setFormatter(logging.Formatter("[%(asctime)s] %(levelname)s %(message)s"))
        log.addHandler(h)
    if logfile:
        fh = logging.FileHandler(logfile)
        fh.setFormatter(logging.Formatter("[%(asctime)s] %(levelname)s %(message)s"))
        log.addHandler(fh)
    return log

# ---- Misc ------------------------------------------------------------------

def require(tool: str):
    """Make sure a CLI binary is on PATH."""
    if shutil.which(tool) is None:
        sys.stderr.write(
            f"ERROR: required executable '{tool}' not found on PATH.\n"
            f"  Did you `conda activate schrodinger-lite`?\n"
        )
        sys.exit(1)

def time_it(fn):
    def wrapper(*a, **kw):
        t0 = time.time()
        r = fn(*a, **kw)
        dt = time.time() - t0
        print(f"  [done in {dt:.1f}s]")
        return r
    return wrapper
