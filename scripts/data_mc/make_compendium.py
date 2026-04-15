#!/usr/bin/env python3
"""
make_compendium.py

Generates a LaTeX result compendium PDF for each channel × scenario combination.
Each compendium includes all plots and the per-eta chi-squared table.

Usage:
    cd scripts/data_mc
    python3 make_compendium.py            # builds all 9
    python3 make_compendium.py eegamma baseline  # builds one
"""
import os
import re
import subprocess
import sys

# ── Configuration ──────────────────────────────────────────────────────
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "..", "..", "output",
                    "cell_energy_reweighting_Francisco_method", "data24")
REPORT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "..", "report")

CHANNELS  = ["eegamma", "mumugamma", "llgamma"]
SCENARIOS = ["baseline", "converted", "all_conv"]

CHAN_TEX  = {
    "eegamma":   r"$Z \to ee\gamma$",
    "mumugamma": r"$Z \to \mu\mu\gamma$",
    "llgamma":   r"$Z \to \ell\ell\gamma$",
}
SCEN_TEX = {
    "baseline":  "Unconverted (baseline)",
    "converted": "Converted",
    "all_conv":  "All conversions",
}

N_PAGES = 13  # pages per multi-page eta-binned PDF (14 bins minus crack)

ETA_LABELS = [
    r"[0.00, 0.20)", r"[0.20, 0.40)", r"[0.40, 0.60)", r"[0.60, 0.80)",
    r"[0.80, 1.00)", r"[1.00, 1.20)", r"[1.20, 1.30)", r"[1.30, 1.37)",
    r"[1.52, 1.60)", r"[1.60, 1.80)", r"[1.80, 2.00)", r"[2.00, 2.20)",
    r"[2.20, 2.40)",
]


# ── Chi-squared table extraction ──────────────────────────────────────
def extract_chi2_table(channel, scenario):
    """Extract the per-eta chi2 table for a specific channel/scenario
    from the combined chi2_tables.tex file."""
    chi2_file = os.path.join(REPORT_DIR, "chi2_tables.tex")
    if not os.path.exists(chi2_file):
        return "% chi2\\_tables.tex not found\n"

    with open(chi2_file) as f:
        content = f.read()

    # Split into individual table environments
    tables = re.findall(r"(\\begin\{table\}.*?\\end\{table\})", content, re.DOTALL)

    label_target = f"tab:chi2-{channel}-{scenario}"
    for table in tables:
        if label_target in table:
            return table

    return f"% Table for {channel}/{scenario} not found\n"


# ── LaTeX document generation ────────────────────────────────────────
def generate_tex(channel, scenario):
    """Generate the full LaTeX source for one compendium."""
    plots_dir = os.path.join(BASE, channel, scenario, "plots")
    # graphicspath relative to the .tex file location (channel/scenario/)
    rel_plots = "plots"

    chan_label = CHAN_TEX[channel]
    scen_label = SCEN_TEX[scenario]

    chi2_table = extract_chi2_table(channel, scenario)
    # Replace \chisq with inline definition since we define it in preamble
    # No need — we define \chisq in preamble

    lines = []
    lines.append(r"""\documentclass[11pt,a4paper]{article}
\usepackage[margin=1.8cm]{geometry}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{hyperref}
\usepackage{grffile}
\usepackage{float}

\newcommand{\abseta}{|\eta|}
\newcommand{\chisq}{\chi^2/n_{\mathrm{df}}}

\graphicspath{{%(plots_dir)s/}}

\title{Result Compendium\\[0.3em]
  \large %(chan)s\ --- \ %(scen)s}
\author{Cell-Energy Reweighting Pipeline}
\date{\today}

\begin{document}
\maketitle
\tableofcontents
\clearpage
""" % {"plots_dir": rel_plots, "chan": chan_label, "scen": scen_label})

    # ── Section 1: Chi-squared table ──────────────────────────────
    lines.append(r"""
%% ====================================================================
\section{Chi-Squared Summary}
\label{sec:chi2}
%% ====================================================================

""")
    lines.append(chi2_table)
    lines.append("\n\\clearpage\n")

    # ── Section: Integrated shower shape comparisons ──────────────
    lines.append(r"""
%% ====================================================================
\section{Integrated Shower Shape Comparisons}
\label{sec:shower-integrated}
%% ====================================================================

Data versus uncorrected MC, M1, and M2 corrections (all $|\eta|$ bins combined).

\begin{figure}[H]
\centering
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth,page=1]{rew_integrated.pdf}
  \caption{$R_{\eta}$}
\end{subfigure}\hfill
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth,page=2]{rew_integrated.pdf}
  \caption{$R_{\phi}$}
\end{subfigure}\hfill
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth,page=3]{rew_integrated.pdf}
  \caption{$w_{\eta 2}$}
\end{subfigure}
\caption{Integrated shower shape comparisons: Data vs uncorrected MC, M1 and M2.}
\end{figure}
\clearpage
""")

    # ── Section: Per-eta shower shape comparisons ─────────────────
    lines.append(r"""
%% ====================================================================
\section{Per-$|\eta|$ Shower Shape Comparisons}
\label{sec:shower}
%% ====================================================================

Data versus uncorrected MC, M1, and M2 corrections for each $|\eta|$ bin.
""")

    for page_idx in range(1, N_PAGES + 1):
        eta = ETA_LABELS[page_idx - 1]
        lines.append(r"""
\subsection*{$|\eta| \in %s$}
\begin{figure}[H]
\centering
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth,page=%d]{rew_reta.pdf}
  \caption{$R_{\eta}$}
\end{subfigure}\hfill
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth,page=%d]{rew_rphi.pdf}
  \caption{$R_{\phi}$}
\end{subfigure}\hfill
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth,page=%d]{rew_weta2.pdf}
  \caption{$w_{\eta 2}$}
\end{subfigure}
\caption{Shower shape comparisons, $|\eta| \in %s$.}
\end{figure}
""" % (eta, page_idx, page_idx, page_idx, eta))

    lines.append("\\clearpage\n")

    # ── Section: Fudge factor benchmark (integrated) ─────────────
    lines.append(r"""
%% ====================================================================
\section{Fudge Factor Comparison (Integrated)}
\label{sec:fudge}
%% ====================================================================

Stored (fudged) branch values versus unfudged MC and data.

\begin{figure}[H]
\centering
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth,page=1]{fudge_factors.pdf}
  \caption{$R_{\eta}$}
\end{subfigure}\hfill
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth,page=2]{fudge_factors.pdf}
  \caption{$R_{\phi}$}
\end{subfigure}\hfill
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth,page=3]{fudge_factors.pdf}
  \caption{$w_{\eta 2}$}
\end{subfigure}
\caption{Fudge factor comparison: Data (stored) vs MC (fudged) vs MC (unfudged).}
\end{figure}
\clearpage
""")

    var_names = [r"$R_{\eta}$", r"$R_{\phi}$", r"$w_{\eta 2}$"]

    # ── Section: Fudge factor per-eta ────────────────────────────
    lines.append(r"""
%% ====================================================================
\section{Fudge Factor Comparison (Per $|\eta|$ Bin)}
\label{sec:fudge-eta}
%% ====================================================================

Data (stored) vs MC (fudged) vs MC (unfudged) per $|\eta|$ bin.
""")

    for page_idx in range(1, N_PAGES + 1):
        eta = ETA_LABELS[page_idx - 1]
        pg_reta  = page_idx
        pg_rphi  = page_idx + N_PAGES
        pg_weta2 = page_idx + 2 * N_PAGES
        lines.append(r"""
\subsection*{$|\eta| \in %s$}
\begin{figure}[H]
\centering
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth,page=%d]{fudge_factors_eta.pdf}
  \caption{$R_{\eta}$}
\end{subfigure}\hfill
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth,page=%d]{fudge_factors_eta.pdf}
  \caption{$R_{\phi}$}
\end{subfigure}\hfill
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth,page=%d]{fudge_factors_eta.pdf}
  \caption{$w_{\eta 2}$}
\end{subfigure}
\caption{Fudge factors, $|\eta| \in %s$.}
\end{figure}
""" % (eta, pg_reta, pg_rphi, pg_weta2, eta))

    lines.append("\\clearpage\n")

    # ── Section: Computed vs Stored (integrated) ─────────────────
    lines.append(r"""
%% ====================================================================
\section{Computed vs.\ Stored Validation (Integrated)}
\label{sec:comp-vs-stor}
%% ====================================================================

Shower shapes computed from cell energies versus the branch-stored values.

\begin{figure}[H]
\centering
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth,page=1]{computed_vs_stored.pdf}
  \caption{$R_{\eta}$}
\end{subfigure}\hfill
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth,page=2]{computed_vs_stored.pdf}
  \caption{$R_{\phi}$}
\end{subfigure}\hfill
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth,page=3]{computed_vs_stored.pdf}
  \caption{$w_{\eta 2}$}
\end{subfigure}
\caption{Computed (from cells) vs stored (branch) shower shapes.}
\end{figure}
\clearpage
""")

    # ── Section: Computed vs Stored per-eta ──────────────────────
    lines.append(r"""
%% ====================================================================
\section{Computed vs.\ Stored (Per $|\eta|$ Bin)}
\label{sec:comp-vs-stor-eta}
%% ====================================================================

Cell-computed vs branch-stored data per $|\eta|$ bin.
""")

    for page_idx in range(1, N_PAGES + 1):
        eta = ETA_LABELS[page_idx - 1]
        pg_reta  = page_idx
        pg_rphi  = page_idx + N_PAGES
        pg_weta2 = page_idx + 2 * N_PAGES
        lines.append(r"""
\subsection*{$|\eta| \in %s$}
\begin{figure}[H]
\centering
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth,page=%d]{computed_vs_stored_eta.pdf}
  \caption{$R_{\eta}$}
\end{subfigure}\hfill
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth,page=%d]{computed_vs_stored_eta.pdf}
  \caption{$R_{\phi}$}
\end{subfigure}\hfill
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth,page=%d]{computed_vs_stored_eta.pdf}
  \caption{$w_{\eta 2}$}
\end{subfigure}
\caption{Computed vs.\ stored, $|\eta| \in %s$.}
\end{figure}
""" % (eta, pg_reta, pg_rphi, pg_weta2, eta))

    lines.append("\\clearpage\n")

    # ── Section 7: Cell energy profiles ───────────────────────────
    lines.append(r"""
%% ====================================================================
\section{Cell Energy Fraction Profiles}
\label{sec:cell-profiles}
%% ====================================================================

Mean cell energy fraction $\langle f_k \rangle$ maps in the $7 \times 11$
$(\eta \times \phi)$ grid.  Each page shows one $|\eta|$ bin.
""")

    for page_idx in range(1, N_PAGES + 1):
        eta = ETA_LABELS[page_idx - 1]
        lines.append(r"""
\subsection*{$|\eta| \in %s$}
\begin{figure}[H]
\centering
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth,page=%d]{cell_data.pdf}
  \caption{Data}
\end{subfigure}\hfill
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth,page=%d]{cell_mc.pdf}
  \caption{MC (uncorrected)}
\end{subfigure}\\[1em]
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth,page=%d]{cell_mc_m1.pdf}
  \caption{MC after M1}
\end{subfigure}\hfill
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth,page=%d]{cell_mc_m2.pdf}
  \caption{MC after M2}
\end{subfigure}
\caption{Cell energy fraction maps, $|\eta| \in %s$.}
\end{figure}
""" % (eta, page_idx, page_idx, page_idx, page_idx, eta))

    lines.append("\\clearpage\n")

    # ── Section 8: Cell corrections ───────────────────────────────
    lines.append(r"""
%% ====================================================================
\section{Cell-Level Corrections}
\label{sec:cell-corrections}
%% ====================================================================

M1 shift ($\Delta_k$) and M2 stretch ($s_k$) maps per $|\eta|$ bin.
""")

    for page_idx in range(1, N_PAGES + 1):
        eta = ETA_LABELS[page_idx - 1]
        lines.append(r"""
\subsection*{$|\eta| \in %s$}
\begin{figure}[H]
\centering
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth,page=%d]{cell_shift.pdf}
  \caption{M2 shift ($a_k$)}
\end{subfigure}\hfill
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth,page=%d]{cell_stretch.pdf}
  \caption{M2 stretch ($s_k$)}
\end{subfigure}
\caption{Correction maps, $|\eta| \in %s$.}
\end{figure}
""" % (eta, page_idx, page_idx, eta))

    lines.append(r"""

%% ====================================================================
\section{Kinematics}
\label{sec:kinematics}
%% ====================================================================

\begin{figure}[H]
\centering
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{kinematics_pt.pdf}
  \caption{Photon $p_{\mathrm{T}}$}
\end{subfigure}\hfill
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{kinematics_eta.pdf}
  \caption{Photon $|\eta|$}
\end{subfigure}
\caption{Normalised kinematic distributions: Data vs MC.}
\end{figure}

\end{document}
""")

    return "".join(lines)


# ── Build one compendium ─────────────────────────────────────────────
def build_compendium(channel, scenario):
    """Generate and compile one compendium PDF."""
    out_dir = os.path.join(BASE, channel, scenario)
    tex_name = f"result_compendium_{channel}_{scenario}.tex"
    pdf_name = f"result_compendium_{channel}_{scenario}.pdf"
    tex_path = os.path.join(out_dir, tex_name)

    print(f"\n{'='*60}")
    print(f"  Building: {channel} / {scenario}")
    print(f"{'='*60}")

    # Generate LaTeX
    tex_src = generate_tex(channel, scenario)
    with open(tex_path, "w") as f:
        f.write(tex_src)
    print(f"  Written: {tex_path}")

    # Compile (twice for ToC)
    for pass_num in (1, 2):
        result = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", tex_name],
            cwd=out_dir, capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0 and pass_num == 2:
            print(f"  WARNING: pdflatex returned non-zero on pass {pass_num}")
            log_file = os.path.join(out_dir, f"result_compendium_{channel}_{scenario}.log")
            if os.path.exists(log_file):
                with open(log_file) as lf:
                    log_lines = lf.readlines()
                for line in log_lines[-30:]:
                    line = line.rstrip()
                    if "!" in line or "Error" in line:
                        print(f"    {line}")

    # Clean auxiliary files
    for ext in ("aux", "log", "out", "toc"):
        aux = os.path.join(out_dir, f"result_compendium_{channel}_{scenario}.{ext}")
        if os.path.exists(aux):
            os.remove(aux)

    pdf_path = os.path.join(out_dir, pdf_name)
    if os.path.exists(pdf_path):
        size_mb = os.path.getsize(pdf_path) / 1024 / 1024
        print(f"  OK: {pdf_path} ({size_mb:.1f} MB)")
        return True
    else:
        print(f"  FAILED: {pdf_path} not produced")
        return False


# ── Main ─────────────────────────────────────────────────────────────
def main():
    if len(sys.argv) == 3:
        channels = [sys.argv[1]]
        scenarios = [sys.argv[2]]
    else:
        channels = CHANNELS
        scenarios = SCENARIOS

    ok, fail = 0, 0
    for ch in channels:
        for sc in scenarios:
            if build_compendium(ch, sc):
                ok += 1
            else:
                fail += 1

    print(f"\n{'='*60}")
    print(f"  Done: {ok} OK, {fail} failed")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
