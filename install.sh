#!/usr/bin/env bash
# SchrodingerLite installer — sets up everything needed to run from anywhere on this Mac.
#
#   • Creates conda env 'schrodinger-lite' with RDKit, AutoDock Vina, OpenBabel,
#     PyMOL, OpenMM, Psi4, PDBFixer, scikit-learn, etc.
#   • Installs the schrodingerlite Python package in editable mode.
#   • Symlinks all launchers (schrodinger, glide, ligprep, …) into a PATH dir,
#     so you can run them from ANY terminal, in ANY folder.
#
# Re-running is safe; it will skip steps already done.

set -e
trap 'echo; echo "✗ Install failed at line $LINENO. See messages above." >&2' ERR

ROOT="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ENV_NAME="schrodinger-lite"
CYAN="\033[36m"; GREEN="\033[32m"; YELLOW="\033[33m"; BOLD="\033[1m"; NC="\033[0m"

echo -e "${BOLD}${CYAN}"
echo "================================================================="
echo "  SchrodingerLite — Local install"
echo "  Install root: $ROOT"
echo "=================================================================${NC}"

# --- 1. Find conda --------------------------------------------------------
CONDA_BIN=""
for c in "$HOME/miniforge3/bin/conda" "$HOME/miniconda3/bin/conda" \
         "$HOME/anaconda3/bin/conda" "/opt/homebrew/Caskroom/miniforge/base/bin/conda" \
         "/opt/homebrew/anaconda3/bin/conda" "/usr/local/anaconda3/bin/conda" \
         "$(command -v conda 2>/dev/null)"; do
  [ -x "$c" ] && CONDA_BIN="$c" && break
done

if [ -z "$CONDA_BIN" ]; then
  echo -e "${YELLOW}conda not found.${NC} Installing Miniforge via Homebrew…"
  if ! command -v brew >/dev/null; then
    echo "Homebrew is required. Install from https://brew.sh and re-run."
    exit 1
  fi
  brew install --cask miniforge
  # Initialize
  CONDA_BIN="/opt/homebrew/Caskroom/miniforge/base/bin/conda"
  [ -x "$CONDA_BIN" ] || CONDA_BIN="$(command -v conda)"
fi
echo -e "${GREEN}✓${NC} conda at $CONDA_BIN"

# shellcheck disable=SC1091
source "$(dirname "$CONDA_BIN")/../etc/profile.d/conda.sh"

# --- 2. Enable libmamba solver (much faster + better at conflicts) --------
conda config --set solver libmamba 2>/dev/null || {
  echo "Installing libmamba solver (faster dependency resolution)…"
  conda install -n base conda-libmamba-solver -y -c conda-forge >/dev/null 2>&1 || true
  conda config --set solver libmamba 2>/dev/null || true
}

# --- 3. Create env --------------------------------------------------------
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo -e "${YELLOW}!${NC} env '$ENV_NAME' exists. Removing any partial install and rebuilding…"
  conda env remove -n "$ENV_NAME" -y
fi
echo -e "${CYAN}→${NC} Creating conda env '$ENV_NAME' (5-10 min with libmamba)…"
conda env create -n "$ENV_NAME" -f "$ROOT/environment.yml"
conda activate "$ENV_NAME"
echo -e "${GREEN}✓${NC} env '$ENV_NAME' active ($(python --version))"

# --- 3. Install the schrodingerlite Python package ------------------------
echo -e "${CYAN}→${NC} Installing schrodingerlite (editable)…"
pip install -e "$ROOT" --quiet
echo -e "${GREEN}✓${NC} schrodingerlite installed"

# --- 4. Symlink launchers into a PATH dir ---------------------------------
# Prefer /usr/local/bin if writable, else ~/.local/bin (auto-add to PATH)
TARGET="/usr/local/bin"
if [ ! -w "$TARGET" ]; then
  TARGET="$HOME/.local/bin"
  mkdir -p "$TARGET"
fi
echo -e "${CYAN}→${NC} Symlinking launchers into $TARGET"
for cmd in schrodinger maestro prepwizard ligprep glide prime desmond jaguar macromodel canvas strike aidiscover; do
  ln -sf "$ROOT/bin/$cmd" "$TARGET/$cmd"
done

# Ensure ~/.local/bin is on PATH (zsh by default on modern macOS)
if [ "$TARGET" = "$HOME/.local/bin" ]; then
  for rc in "$HOME/.zshrc" "$HOME/.bash_profile" "$HOME/.bashrc"; do
    [ -f "$rc" ] || continue
    if ! grep -q 'SchrodingerLite PATH' "$rc"; then
      {
        echo ""
        echo "# SchrodingerLite PATH"
        echo "export PATH=\"\$HOME/.local/bin:\$PATH\""
      } >> "$rc"
      echo -e "${GREEN}✓${NC} added ~/.local/bin to PATH in $(basename "$rc")"
    fi
  done
fi

# --- 5. Done --------------------------------------------------------------
echo -e "${GREEN}"
echo "================================================================="
echo "  ✓ SchrodingerLite is installed."
echo "================================================================="
echo -e "${NC}"
echo "Open a NEW terminal window and try:"
echo "  schrodinger version"
echo "  schrodinger --help"
echo ""
echo "Run the full demo:"
echo "  schrodinger workflow run \"$ROOT/examples/hiv_protease_demo.yml\""
echo ""
echo "All outputs will be saved under:"
echo "  $ROOT/projects/"
