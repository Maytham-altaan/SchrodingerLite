"""AI Design — Claude-powered de novo molecule generation.

Calls the Anthropic Claude API to invent novel drug-like SMILES for a target,
with a med-chem rationale attached to each molecule.

Requires ANTHROPIC_API_KEY in the environment.
"""
from __future__ import annotations
import json
import os
import re
import time
from pathlib import Path
from typing import Iterable

from .utils import get_logger


CLAUDE_MODEL = os.environ.get("SCHRODINGER_LITE_CLAUDE_MODEL",
                              "claude-sonnet-4-6")

SYSTEM_PROMPT = """You are an expert medicinal chemist designing novel drug candidates.

When asked to design molecules, you MUST output strict JSON of the form:

{
  "molecules": [
    {
      "name": "short_unique_name",
      "smiles": "<valid canonical SMILES>",
      "rationale": "<one sentence: which interaction, scaffold choice, why drug-like>",
      "expected_advantages": ["<2-4 short bullet phrases>"]
    },
    ...
  ]
}

Rules:
- Output ONLY JSON. No markdown fences, no commentary outside JSON.
- All SMILES must be VALID and drug-like (MW 200-500, ClogP -1 to 5, HBD≤5, HBA≤10, rotB≤10).
- Avoid PAINS, reactive groups, metal centers, peptide-like chains >3 residues.
- Each molecule must be NOVEL — not a known drug, not a trivial atom swap of one.
- Diversify scaffolds: different ring systems, different vectors into the binding site.
- Prefer synthesizable chemistry (avoid 5+ stereocenters, unusual heterocycles).
"""


# ---------------------------------------------------------------------------

def _claude_client():
    try:
        import anthropic
    except ImportError as e:
        raise RuntimeError(
            "anthropic SDK not installed. Run:\n"
            "  conda activate schrodinger-lite && pip install anthropic"
        ) from e
    if not os.environ.get("ANTHROPIC_API_KEY"):
        raise RuntimeError(
            "ANTHROPIC_API_KEY not set. Get a key at https://console.anthropic.com\n"
            "Then in your terminal:\n"
            "  export ANTHROPIC_API_KEY='sk-ant-…'\n"
            "(add it to ~/.zshrc to persist)."
        )
    return anthropic.Anthropic()


def _parse_json(text: str) -> dict:
    """Extract JSON from Claude's reply, tolerating accidental fences."""
    # Strip markdown fences if present
    m = re.search(r"```(?:json)?\s*(\{.*?\})\s*```", text, re.DOTALL)
    if m:
        text = m.group(1)
    # Find the first complete JSON object
    start = text.find("{")
    end = text.rfind("}")
    if start < 0 or end < 0:
        raise ValueError(f"No JSON found in Claude response:\n{text[:500]}")
    return json.loads(text[start:end+1])


# ---------------------------------------------------------------------------

def generate_round(target_description: str,
                   n_molecules: int = 30,
                   feedback: list[dict] | None = None,
                   round_num: int = 1,
                   max_tokens: int = 8000) -> list[dict]:
    """Ask Claude to design molecules for a target.

    Args:
        target_description: free-text description of the target / binding site.
        n_molecules: how many to ask for in this round.
        feedback: list of dicts with prior results, each like
                  {"smiles": ..., "score": ..., "notes": ...}
        round_num: iteration number, used in the prompt.

    Returns:
        List of dicts: [{name, smiles, rationale, expected_advantages, round}, …]
    """
    log = get_logger("ai_design")
    client = _claude_client()

    if feedback:
        feedback_block = "\n".join(
            f"- {f['smiles']}  (score {f.get('score','?')} kcal/mol)  "
            f"{f.get('notes','')}"
            for f in feedback
        )
        user_msg = (
            f"Round {round_num}. TARGET:\n{target_description}\n\n"
            f"PREVIOUS RESULTS (top hits from last round):\n{feedback_block}\n\n"
            f"Now design {n_molecules} NEW molecules. Build on what worked, "
            f"avoid what didn't. Diversify scaffolds — do not just analogue the "
            f"top hit. Output the JSON schema described in the system prompt."
        )
    else:
        user_msg = (
            f"Round {round_num}. TARGET:\n{target_description}\n\n"
            f"Design {n_molecules} novel drug candidates. Diversify scaffolds — "
            f"include at least 4 distinct chemotypes. Output the JSON schema."
        )

    log.info(f"Calling Claude ({CLAUDE_MODEL}) for round {round_num}, "
             f"{n_molecules} molecules…")
    t0 = time.time()
    resp = client.messages.create(
        model=CLAUDE_MODEL,
        max_tokens=max_tokens,
        system=SYSTEM_PROMPT,
        messages=[{"role": "user", "content": user_msg}],
    )
    dt = time.time() - t0
    log.info(f"Claude responded in {dt:.1f}s "
             f"({resp.usage.input_tokens} in / {resp.usage.output_tokens} out)")

    text = "".join(b.text for b in resp.content if b.type == "text")
    try:
        data = _parse_json(text)
    except Exception as e:
        log.error(f"Failed to parse Claude response: {e}")
        log.error(f"Raw text:\n{text[:2000]}")
        return []

    out = []
    for m in data.get("molecules", []):
        if "smiles" not in m:
            continue
        m["round"] = round_num
        out.append(m)
    log.info(f"Round {round_num}: Claude proposed {len(out)} molecules")
    return out


# ---------------------------------------------------------------------------

def explain_results(top_hits: list[dict], target_description: str,
                    max_tokens: int = 2000) -> str:
    """Ask Claude to write a med-chem summary of the top hits."""
    client = _claude_client()
    listing = "\n".join(
        f"{i+1}. {h.get('name','?')}  {h['smiles']}  "
        f"ΔG={h.get('score','?')} kcal/mol  "
        f"rationale: {h.get('rationale','')}"
        for i, h in enumerate(top_hits[:10])
    )
    msg = (
        f"TARGET:\n{target_description}\n\n"
        f"TOP HITS from de novo discovery:\n{listing}\n\n"
        f"Write a 2-3 paragraph medicinal-chemistry summary: which scaffolds "
        f"dominate the top, what binding-site features they exploit, what "
        f"liabilities to watch for, and what to test experimentally first. "
        f"Plain text, no markdown."
    )
    resp = client.messages.create(
        model=CLAUDE_MODEL,
        max_tokens=max_tokens,
        messages=[{"role": "user", "content": msg}],
    )
    return "".join(b.text for b in resp.content if b.type == "text")
