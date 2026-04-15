#!/usr/bin/env python3
"""
make_kinematics_compendium.py

Generates a LaTeX compendium of photon kinematic distributions
(pT and |eta|, normalised and unnormalised) for all scenarios produced
by run_kinematics.sh.

Output structure expected:
  output/kinematics/{loose,iso_tight}/llgamma/{baseline,converted,all_conv}/plots/
    kinematics_pt.pdf, kinematics_eta.pdf         — unnormalised
    kinematics_pt_norm.pdf, kinematics_eta_norm.pdf — normalised

The compendium is written to:
  output/kinematics/kinematics_compendium.pdf

Usage:
    python3 make_kinematics_compendium.py
    python3 make_kinematics_compendium.py --base-dir path/to/kinematics
"""
import argparse
import os
import subprocess
import sys

# ── Configuration ──────────────────────────────────────────────────────
DEFAULT_BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "..", "..", "output", "kinematics")

CHANNEL = "llgamma"

# Ordered list: (isolation label, base subdir, conversion, selectionLabel)
SCENARIOS = [
    ("Loose isolation", "loose", "baseline",  "Unconverted $\\gamma$, loose iso"),
    ("Loose isolation", "loose", "converted", "Converted $\\gamma$, loose iso"),
    ("Loose isolation", "loose", "all_conv",  "All $\\gamma$, loose iso"),
    ("Tight isolation", "iso_tight", "baseline",  "Unconverted $\\gamma$, tight iso"),
    ("Tight isolation", "iso_tight", "converted", "Converted $\\gamma$, tight iso"),
    ("Tight isolation", "iso_tight", "all_conv",  "All $\\gamma$, tight iso"),
]


def generate_tex(base_dir):
    """Generate the full LaTeX source for the kinematics compendium."""
    lines = []
    lines.append(r"""\documentclass[11pt,a4paper]{article}
\usepackage[margin=1.8cm]{geometry}
\usepackage{graphicx}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{hyperref}
\usepackage{grffile}
\usepackage{float}

\title{Kinematics Compendium\\[0.3em]
  \large $Z \to \ell\ell\gamma$ --- Photon $p_{\mathrm{T}}$ and $|\eta|$ Distributions}
\author{Cell-Energy Reweighting Pipeline}
\date{\today}

\begin{document}
\maketitle
\tableofcontents
\clearpage
""")

    current_iso = None

    for iso_label, iso_dir, conv, scen_label in SCENARIOS:
        plots_path = os.path.join(base_dir, iso_dir, CHANNEL, conv, "plots")
        # graphicspath relative to tex file location (base_dir/)
        rel_plots = os.path.relpath(plots_path, base_dir)

        # Start new section for each isolation group
        if iso_label != current_iso:
            current_iso = iso_label
            sec_label = iso_label.lower().replace(" ", "-")
            lines.append(r"""
%% ====================================================================
\section{%(iso)s}
\label{sec:%(label)s}
%% ====================================================================
""" % {"iso": iso_label, "label": sec_label})

        # Subsection per conversion type
        conv_tex = {"baseline": "Unconverted",
                    "converted": "Converted",
                    "all_conv": "All conversions"}[conv]

        lines.append(r"""
\subsection{%(conv)s}

\begin{figure}[H]
\centering
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{%(rel)s/kinematics_pt.pdf}
  \caption{$p_{\mathrm{T}}$ (events)}
\end{subfigure}\hfill
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{%(rel)s/kinematics_pt_norm.pdf}
  \caption{$p_{\mathrm{T}}$ (normalised)}
\end{subfigure}\\[1em]
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{%(rel)s/kinematics_eta.pdf}
  \caption{$|\eta|$ (events)}
\end{subfigure}\hfill
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{%(rel)s/kinematics_eta_norm.pdf}
  \caption{$|\eta|$ (normalised)}
\end{subfigure}
\caption{%(scen)s.  Top: photon $p_{\mathrm{T}}$; bottom: photon $|\eta|$.
  Left: event counts; right: area-normalised.}
\end{figure}
\clearpage
""" % {"conv": conv_tex, "rel": rel_plots, "scen": scen_label})

    lines.append(r"\end{document}" + "\n")
    return "".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Build a kinematics compendium PDF.")
    parser.add_argument("--base-dir", default=DEFAULT_BASE,
                        help="Base kinematics output directory")
    args = parser.parse_args()

    base_dir = os.path.abspath(args.base_dir)
    tex_name = "kinematics_compendium.tex"
    pdf_name = "kinematics_compendium.pdf"
    tex_path = os.path.join(base_dir, tex_name)

    print(f"\n{'='*60}")
    print(f"  Building kinematics compendium")
    print(f"  Base: {base_dir}")
    print(f"{'='*60}")

    os.makedirs(base_dir, exist_ok=True)
    tex_src = generate_tex(base_dir)
    with open(tex_path, "w") as f:
        f.write(tex_src)
    print(f"  Written: {tex_path}")

    # Compile (twice for ToC)
    for pass_num in (1, 2):
        result = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", tex_name],
            cwd=base_dir, capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0 and pass_num == 2:
            print(f"  WARNING: pdflatex returned non-zero on pass {pass_num}")
            log_file = os.path.join(base_dir, "kinematics_compendium.log")
            if os.path.exists(log_file):
                with open(log_file) as lf:
                    log_lines = lf.readlines()
                for line in log_lines[-30:]:
                    line = line.rstrip()
                    if "!" in line or "Error" in line:
                        print(f"    {line}")

    # Clean auxiliary files
    for ext in ("aux", "log", "out", "toc"):
        aux = os.path.join(base_dir, f"kinematics_compendium.{ext}")
        if os.path.exists(aux):
            os.remove(aux)

    pdf_path = os.path.join(base_dir, pdf_name)
    if os.path.exists(pdf_path):
        size_mb = os.path.getsize(pdf_path) / 1024 / 1024
        print(f"  OK: {pdf_path} ({size_mb:.1f} MB)")
    else:
        print(f"  FAILED: {pdf_path} not produced")

    print(f"{'='*60}")


if __name__ == "__main__":
    main()
