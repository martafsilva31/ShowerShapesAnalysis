#!/usr/bin/env python3
"""
Generate and compile result_compendium PDFs for hist_sigma_cellrange and hist_sigma_wide.

Generates the full LaTeX from scratch (no external template needed).
Handles fudge_factors_eta.pdf having 1 page (summary) vs 39 (3x13 per-bin).
Extracts the correct chi2 table per channel/scenario from report/chi2_tables.tex.
"""

import os, re, subprocess, sys


def pdf_ok(path):
    """Return True if the PDF file exists and is not corrupted."""
    if not os.path.isfile(path):
        return False
    r = subprocess.run(["pdfinfo", path], capture_output=True, text=True)
    return r.returncode == 0


def pdf_pages(path):
    """Return the number of pages in a PDF, or 0 if invalid."""
    if not pdf_ok(path):
        return 0
    r = subprocess.run(["pdfinfo", path], capture_output=True, text=True)
    for line in r.stdout.splitlines():
        if line.startswith("Pages:"):
            return int(line.split()[-1])
    return 0

BASE     = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../output/Layer_2"))
PDFLATEX = "/cvmfs/sft.cern.ch/lcg/external/texlive/latest/bin/x86_64-linux/pdflatex"

VARIANTS = {
    "eta_loose":     r"$\eta$-only binning, loose isolation",
    "eta_tight":     r"$\eta$-only binning, tight isolation",
    "eta_pt_loose":  r"$\eta \times p_{\mathrm{T}}$ binning, loose isolation",
    "eta_pt_tight":  r"$\eta \times p_{\mathrm{T}}$ binning, tight isolation",
    "eta_mu_loose": r"$\eta \times \langle\mu\rangle$ binning, loose isolation",
}

SCENARIOS = ["converted", "inclusive", "unconverted"]

SCENARIO_LABELS = {
    "converted":   r"Converted $\gamma$",
    "unconverted": r"Unconverted $\gamma$",
    "inclusive":   r"Converted and Unconverted $\gamma$",
}

ETA_BINS = [
    (0.00, 0.20), (0.20, 0.40), (0.40, 0.60), (0.60, 0.80),
    (0.80, 1.00), (1.00, 1.20), (1.20, 1.30), (1.30, 1.37),
    (1.52, 1.60), (1.60, 1.80), (1.80, 2.00), (2.00, 2.20),
    (2.20, 2.40),
]
NBINS = len(ETA_BINS)  # 13

PT_BINS = [
    (10, 15), (15, 20), (20, 25), (25, 30), (30, 40), (40, 1000),
]
NPT = len(PT_BINS)  # 6

def pt_label(lo, hi):
    if hi >= 1000:
        return rf"$p_{{\mathrm{{T}}}} > {lo}$~GeV"
    return rf"${lo} < p_{{\mathrm{{T}}}} < {hi}$~GeV"

MU_BINS = [
    (0, 40), (40, 55), (55, 70),
]
NMU = len(MU_BINS)  # 3

def mu_label(lo, hi):
    if hi >= 119:
        return rf"$\langle\mu\rangle > {lo}$"
    return rf"${lo} < \langle\mu\rangle < {hi}$"

CHANNEL_LATEX = {
    "llgamma": r"$Z \to \ell\ell\gamma$",
    "eegamma": r"$Z \to ee\gamma$",
    "mumugamma": r"$Z \to \mu\mu\gamma$",
}

# ── helpers ───────────────────────────────────────────────────────────────────

def eta_label(lo, hi):
    return rf"$|\eta| \in [{lo:.2f}, {hi:.2f})$"


def extract_chi2_table(chi2_tables_path, channel, scenario):
    """Extract the channel/scenario table block from chi2_tables.tex."""
    label = f"tab:chi2-{channel}-{scenario}"
    with open(chi2_tables_path) as f:
        text = f.read()
    pattern = r"(\\begin\{table\}.*?\\end\{table\})"
    for m in re.finditer(pattern, text, re.DOTALL):
        if label in m.group(0):
            return m.group(0)
    return None


def extract_chi2_pt_tables(chi2_pt_path, scenario):
    """Extract all per-(|eta|,pT) chi2 table blocks for a given scenario
    from chi2_tables_eta_pt.tex (written by extract_chi2.C for eta_pt variants).
    Returns a list of LaTeX table strings (one per shower-shape variable).
    """
    if not os.path.isfile(chi2_pt_path):
        return []
    with open(chi2_pt_path) as f:
        text = f.read()
    pattern = r"(\\begin\{table\}.*?\\end\{table\})"
    result = []
    for m in re.finditer(pattern, text, re.DOTALL):
        if f"tab:chi2-etapt-{scenario}-" in m.group(0):
            result.append(m.group(0))
    return result


def three_panel(pdf, p1, p2, p3, caption, width="0.32"):
    """Three side-by-side subfigures from pages p1, p2, p3 of a multi-page PDF."""
    return rf"""
\begin{{figure}}[H]
\centering
\begin{{subfigure}}[b]{{{width}\textwidth}}
  \includegraphics[width=\textwidth,page={p1}]{{{pdf}}}
  \caption{{$R_{{\eta}}$}}
\end{{subfigure}}\hfill
\begin{{subfigure}}[b]{{{width}\textwidth}}
  \includegraphics[width=\textwidth,page={p2}]{{{pdf}}}
  \caption{{$R_{{\phi}}$}}
\end{{subfigure}}\hfill
\begin{{subfigure}}[b]{{{width}\textwidth}}
  \includegraphics[width=\textwidth,page={p3}]{{{pdf}}}
  \caption{{$w_{{\eta 2}}$}}
\end{{subfigure}}
\caption{{{caption}}}
\end{{figure}}
"""


def three_panel_3pdfs(pdf1, pdf2, pdf3, page, caption, width="0.32"):
    """Three side-by-side subfigures from the same page of three different PDFs."""
    return rf"""
\begin{{figure}}[H]
\centering
\begin{{subfigure}}[b]{{{width}\textwidth}}
  \includegraphics[width=\textwidth,page={page}]{{{pdf1}}}
  \caption{{$R_{{\eta}}$}}
\end{{subfigure}}\hfill
\begin{{subfigure}}[b]{{{width}\textwidth}}
  \includegraphics[width=\textwidth,page={page}]{{{pdf2}}}
  \caption{{$R_{{\phi}}$}}
\end{{subfigure}}\hfill
\begin{{subfigure}}[b]{{{width}\textwidth}}
  \includegraphics[width=\textwidth,page={page}]{{{pdf3}}}
  \caption{{$w_{{\eta 2}}$}}
\end{{subfigure}}
\caption{{{caption}}}
\end{{figure}}
"""


def two_panel(pdf1, pdf2, page, cap1, cap2, caption, width="0.48"):
    """Two side-by-side subfigures from the same page of two PDFs."""
    return rf"""
\begin{{figure}}[H]
\centering
\begin{{subfigure}}[b]{{{width}\textwidth}}
  \includegraphics[width=\textwidth,page={page}]{{{pdf1}}}
  \caption{{{cap1}}}
\end{{subfigure}}\hfill
\begin{{subfigure}}[b]{{{width}\textwidth}}
  \includegraphics[width=\textwidth,page={page}]{{{pdf2}}}
  \caption{{{cap2}}}
\end{{subfigure}}
\caption{{{caption}}}
\end{{figure}}
"""


def four_panel(pdf1, pdf2, pdf3, pdf4, page, caps, caption, width="0.48"):
    """2x2 grid from the same page of four PDFs."""
    return rf"""
\begin{{figure}}[H]
\centering
\begin{{subfigure}}[b]{{{width}\textwidth}}
  \includegraphics[width=\textwidth,page={page}]{{{pdf1}}}
  \caption{{{caps[0]}}}
\end{{subfigure}}\hfill
\begin{{subfigure}}[b]{{{width}\textwidth}}
  \includegraphics[width=\textwidth,page={page}]{{{pdf2}}}
  \caption{{{caps[1]}}}
\end{{subfigure}}\\[1em]
\begin{{subfigure}}[b]{{{width}\textwidth}}
  \includegraphics[width=\textwidth,page={page}]{{{pdf3}}}
  \caption{{{caps[2]}}}
\end{{subfigure}}\hfill
\begin{{subfigure}}[b]{{{width}\textwidth}}
  \includegraphics[width=\textwidth,page={page}]{{{pdf4}}}
  \caption{{{caps[3]}}}
\end{{subfigure}}
\caption{{{caption}}}
\end{{figure}}
"""


# ── build full tex ────────────────────────────────────────────────────────────

def build_tex(channel, scenario, variant_desc, chi2_table, plots_dir, report_dir="", variant=""):
    ch_label = CHANNEL_LATEX.get(channel, channel)
    sc_title = SCENARIO_LABELS.get(scenario, scenario.capitalize())

    # Pre-check which plot PDFs are valid
    def have(name):
        return pdf_ok(os.path.join(plots_dir, name))

    L = []
    L.append(r"""\documentclass[11pt,a4paper]{article}
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

\graphicspath{{plots/}}
""")

    L.append(rf"""\title{{Result Compendium\\[0.3em]
  \large {ch_label}\ --- \ {sc_title}\\[0.2em]
  \small {variant_desc}}}
\author{{Cell-Energy Reweighting Pipeline}}
\date{{\today}}

\begin{{document}}
\maketitle
\tableofcontents
\clearpage
""")

    # ── Section 1: Chi-Squared Summary ──
    L.append(r"""
%% ====================================================================
\section{Chi-Squared Summary}
\label{sec:chi2}
%% ====================================================================
""")
    if chi2_table:
        L.append(chi2_table)
    else:
        L.append(r"\emph{Chi-squared table not available for this combination.}")

    # Per-(|eta|,pT) chi2 tables — only for eta_pt variants
    if report_dir:
        chi2_pt_path = os.path.join(report_dir, "chi2_tables_eta_pt.tex")
        pt_tables = extract_chi2_pt_tables(chi2_pt_path, scenario)
        if pt_tables:
            L.append(r"""
%% Per-(eta,pT) chi-squared tables
\subsection*{Per-$(|\eta|, p_{\mathrm{T}})$ $\chi^2/n_{\mathrm{df}}$ Tables}
""")
            for t in pt_tables:
                L.append(t)

    L.append(r"\clearpage")

    # ── Section 2: Integrated Shower Shape Comparisons ──
    if have("rew_integrated.pdf"):
        L.append(r"""
%% ====================================================================
\section{Integrated Shower Shape Comparisons}
\label{sec:shower-integrated}
%% ====================================================================

Data versus uncorrected MC, M1, and M2 corrections (all $|\eta|$ bins combined).
""")
        L.append(three_panel("rew_integrated.pdf", 1, 2, 3,
                             "Integrated shower shape comparisons: Data vs uncorrected MC, M1 and M2."))
        L.append(r"\clearpage")

    # ── Section 3: Per-|eta| Shower Shape Comparisons ──
    L.append(r"""
%% ====================================================================
\section{Per-$|\eta|$ Shower Shape Comparisons}
\label{sec:shower}
%% ====================================================================

Data versus uncorrected MC, M1, and M2 corrections for each $|\eta|$ bin.
""")
    for i, (lo, hi) in enumerate(ETA_BINS):
        page = i + 1
        L.append(rf"\subsection*{{{eta_label(lo, hi)}}}")
        L.append(three_panel_3pdfs("rew_reta.pdf", "rew_rphi.pdf", "rew_weta2.pdf",
                                   page,
                                   rf"Shower shape comparisons, {eta_label(lo, hi)}."))
    L.append(r"\clearpage")

    # ── Section 3b: Per-mu Shower Shape Comparisons (eta_mu only) ──
    if have("rew_reta_mu00.pdf"):
        L.append(r"""
%% ====================================================================
\section{Per-$\langle\mu\rangle$ Shower Shape Comparisons}
\label{sec:shower-mu}
%% ====================================================================

Data versus uncorrected MC, M1, and M2 corrections for each $\langle\mu\rangle$ bin.
""")
        for imu, (mulo, muhi) in enumerate(MU_BINS):
            muf = f"mu{imu:02d}"
            reta_mu = f"rew_reta_{muf}.pdf"
            rphi_mu = f"rew_rphi_{muf}.pdf"
            weta_mu = f"rew_weta2_{muf}.pdf"
            if not have(reta_mu):
                continue
            n_pages_mu = pdf_pages(os.path.join(plots_dir, reta_mu))
            L.append(rf"\subsection{{{mu_label(mulo, muhi)}}}")
            for i, (lo, hi) in enumerate(ETA_BINS):
                page = i + 1
                if page > n_pages_mu:
                    break
                L.append(rf"\subsubsection*{{{eta_label(lo, hi)}}}")
                L.append(three_panel_3pdfs(reta_mu, rphi_mu, weta_mu,
                                           page,
                                           rf"Shower shape comparisons, {mu_label(mulo, muhi)}, {eta_label(lo, hi)}."))
            L.append(r"\clearpage")

    # ── Section 3c: Per-pT Shower Shape Comparisons (eta_pt only) ──
    if "eta_pt" in variant and have("rew_reta_pt00.pdf"):
        L.append(r"""
%% ====================================================================
\section{Per-$p_{\mathrm{T}}$ Shower Shape Comparisons}
\label{sec:shower-pt}
%% ====================================================================

Data versus uncorrected MC, M1, and M2 corrections for each $p_{\mathrm{T}}$ bin.
Each $p_{\mathrm{T}}$ bin contains all $|\eta|$ bins.
""")
        for ip, (ptlo, pthi) in enumerate(PT_BINS):
            ptf = f"pt{ip:02d}"
            reta_pt = f"rew_reta_{ptf}.pdf"
            rphi_pt = f"rew_rphi_{ptf}.pdf"
            weta_pt = f"rew_weta2_{ptf}.pdf"
            if not have(reta_pt):
                continue
            n_pages_pt = pdf_pages(os.path.join(plots_dir, reta_pt))
            L.append(rf"\subsection{{{pt_label(ptlo, pthi)}}}")
            for i, (lo, hi) in enumerate(ETA_BINS):
                page = i + 1
                if page > n_pages_pt:
                    break
                L.append(rf"\subsubsection*{{{eta_label(lo, hi)}}}")
                L.append(three_panel_3pdfs(reta_pt, rphi_pt, weta_pt,
                                           page,
                                           rf"Shower shape comparisons, {pt_label(ptlo, pthi)}, {eta_label(lo, hi)}."))
            L.append(r"\clearpage")

    # ── Section 3d: Shower Shapes vs ATLAS Fudge Factors ──
    if have("rew_vs_fudge_reta.pdf"):
        L.append(r"""
%% ====================================================================
\section{Shower Shape Comparisons vs.\ ATLAS Fudge Factors}
\label{sec:shower-vs-fudge}
%% ====================================================================

Data (cell-computed) versus original MC, M2 cell reweighting, and ATLAS fudge factors
(stored-branch fudged MC).  Provides a direct comparison of the cell-reweighting
correction against the existing ATLAS shower-shape corrections.
""")
        if have("rew_vs_fudge_integrated.pdf"):
            L.append(r"\subsection*{Integrated (all $|\eta|$ bins)}")
            L.append(three_panel("rew_vs_fudge_integrated.pdf", 1, 2, 3,
                                 r"Integrated: Data vs.\ original MC, M2 reweighting, and ATLAS fudged MC."))

        L.append(r"\subsection*{Per $|\eta|$ bin}")
        for i, (lo, hi) in enumerate(ETA_BINS):
            page = i + 1
            L.append(rf"\subsubsection*{{{eta_label(lo, hi)}}}")
            L.append(three_panel_3pdfs("rew_vs_fudge_reta.pdf",
                                       "rew_vs_fudge_rphi.pdf",
                                       "rew_vs_fudge_weta2.pdf",
                                       page,
                                       rf"Shower shapes vs.\ fudge factors, {eta_label(lo, hi)}."))

        if "eta_pt" in variant and have("rew_vs_fudge_reta_pt00.pdf"):
            L.append(r"\subsection*{Per $p_{\mathrm{T}}$ bin}")
            for ip, (ptlo, pthi) in enumerate(PT_BINS):
                ptf = f"pt{ip:02d}"
                reta_vf = f"rew_vs_fudge_reta_{ptf}.pdf"
                rphi_vf = f"rew_vs_fudge_rphi_{ptf}.pdf"
                weta_vf = f"rew_vs_fudge_weta2_{ptf}.pdf"
                if not have(reta_vf):
                    continue
                n_pages_vf = pdf_pages(os.path.join(plots_dir, reta_vf))
                L.append(rf"\subsubsection{{{pt_label(ptlo, pthi)}}}")
                for i, (lo, hi) in enumerate(ETA_BINS):
                    page = i + 1
                    if page > n_pages_vf:
                        break
                    L.append(rf"\paragraph*{{{eta_label(lo, hi)}}}")
                    L.append(three_panel_3pdfs(reta_vf, rphi_vf, weta_vf,
                                               page,
                                               rf"Shower shapes vs.\ fudge, {pt_label(ptlo, pthi)}, {eta_label(lo, hi)}."))
        L.append(r"\clearpage")

    # ── Section 4: Fudge Factor Comparison (Integrated) ──
    L.append(r"""
%% ====================================================================
\section{Fudge Factor Comparison (Integrated)}
\label{sec:fudge}
%% ====================================================================

Stored (fudged) branch values versus unfudged MC and data.
""")
    L.append(three_panel("fudge_factors.pdf", 1, 2, 3,
                         "Fudge factor comparison: Data (stored) vs MC (fudged) vs MC (unfudged)."))
    L.append(r"\clearpage")

    # ── Section 5: Fudge Factor (All bins summary / Per-eta) ──
    n_fudge_eta = pdf_pages(os.path.join(plots_dir, "fudge_factors_eta.pdf")) if have("fudge_factors_eta.pdf") else 0
    if n_fudge_eta >= 3 * NBINS:   # 39 pages → per-eta breakdown
        L.append(r"""
%% ====================================================================
\section{Fudge Factor Comparison (Per $|\eta|$ Bin)}
\label{sec:fudge-eta}
%% ====================================================================

Fudge factor comparisons for each $|\eta|$ bin.
""")
        for i, (lo, hi) in enumerate(ETA_BINS):
            page_reta = i + 1
            page_rphi = i + 1 + NBINS
            page_weta = i + 1 + 2 * NBINS
            L.append(rf"\subsection*{{{eta_label(lo, hi)}}}")
            L.append(three_panel("fudge_factors_eta.pdf", page_reta, page_rphi, page_weta,
                                 rf"Fudge factor comparison, {eta_label(lo, hi)}."))
        L.append(r"\clearpage")
    elif n_fudge_eta >= 1:
        L.append(r"""
%% ====================================================================
\section{Fudge Factor Comparison (All $|\eta|$ Bins)}
\label{sec:fudge-eta}
%% ====================================================================

Fudge factor comparisons for all $|\eta|$ bins.
""")
        L.append(r"""
\begin{figure}[H]
\centering
\includegraphics[width=0.85\textwidth,page=1]{fudge_factors_eta.pdf}
\caption{Fudge factors for all $|\eta|$ bins overlaid.}
\end{figure}
""")
        L.append(r"\clearpage")

    # ── Section 6: Computed vs. Stored (Integrated) ──
    L.append(r"""
%% ====================================================================
\section{Computed vs.\ Stored Validation (Integrated)}
\label{sec:comp-vs-stor}
%% ====================================================================

Shower shapes computed from cell energies versus the branch-stored values.
""")
    L.append(three_panel("computed_vs_stored.pdf", 1, 2, 3,
                         r"Computed (from cells) vs stored (branch) shower shapes."))
    L.append(r"\clearpage")

    # ── Section 7: Computed vs. Stored (Per-eta) ──
    n_cvs_eta = pdf_pages(os.path.join(plots_dir, "computed_vs_stored_eta.pdf")) if have("computed_vs_stored_eta.pdf") else 0
    if n_cvs_eta >= 3 * NBINS:   # 39 pages → per-eta breakdown available
        L.append(r"""
%% ====================================================================
\section{Computed vs.\ Stored (Per $|\eta|$ Bin)}
\label{sec:comp-vs-stor-eta}
%% ====================================================================

Cell-computed vs branch-stored data per $|\eta|$ bin.
""")
        for i, (lo, hi) in enumerate(ETA_BINS):
            page_reta = i + 1
            page_rphi = i + 1 + NBINS
            page_weta = i + 1 + 2 * NBINS
            L.append(rf"\subsection*{{{eta_label(lo, hi)}}}")
            L.append(three_panel("computed_vs_stored_eta.pdf", page_reta, page_rphi, page_weta,
                                 rf"Computed vs.\ stored, {eta_label(lo, hi)}."))
        L.append(r"\clearpage")
    elif n_cvs_eta >= 1:
        L.append(r"""
%% ====================================================================
\section{Computed vs.\ Stored (Per $|\eta|$ Bin)}
\label{sec:comp-vs-stor-eta}
%% ====================================================================

\emph{computed\_vs\_stored\_eta.pdf has only %d page(s); per-$|\eta|$ breakdown requires %d.}
""" % (n_cvs_eta, 3 * NBINS))
        L.append(r"\clearpage")

    # ── Section 8: Cell Energy Fraction Profiles ──
    cell_pdfs = ["cell_data.pdf", "cell_mc.pdf", "cell_mc_m1.pdf", "cell_mc_m2.pdf"]
    valid_cell_pdfs = [p for p in cell_pdfs if have(p)]
    if valid_cell_pdfs:
        L.append(r"""
%% ====================================================================
\section{Cell Energy Fraction Profiles}
\label{sec:cell-profiles}
%% ====================================================================

Mean cell energy fraction $\langle f_k \rangle$ maps in the $7 \times 11$
$(\eta \times \phi)$ grid.  Each page shows one $|\eta|$ bin.
""")
        cell_caps = {"cell_data.pdf": "Data", "cell_mc.pdf": "MC (uncorrected)",
                     "cell_mc_m1.pdf": "MC after M1", "cell_mc_m2.pdf": "MC after M2"}
        for i, (lo, hi) in enumerate(ETA_BINS):
            page = i + 1
            L.append(rf"\subsection*{{{eta_label(lo, hi)}}}")
            if len(valid_cell_pdfs) >= 4:
                L.append(four_panel(valid_cell_pdfs[0], valid_cell_pdfs[1],
                                    valid_cell_pdfs[2], valid_cell_pdfs[3],
                                    page,
                                    [cell_caps[p] for p in valid_cell_pdfs[:4]],
                                    rf"Cell energy fraction maps, {eta_label(lo, hi)}."))
            else:
                # Fall back: include whatever is available
                L.append(r"\begin{figure}[H]")
                L.append(r"\centering")
                w = f"{round(0.96/len(valid_cell_pdfs), 2):.2f}"
                for p in valid_cell_pdfs:
                    L.append(rf"\begin{{subfigure}}[b]{{{w}\textwidth}}")
                    L.append(rf"  \includegraphics[width=\textwidth,page={page}]{{{p}}}")
                    L.append(rf"  \caption{{{cell_caps[p]}}}")
                    L.append(r"\end{subfigure}\hfill")
                L.append(rf"\caption{{Cell energy fraction maps, {eta_label(lo, hi)}.}}")
                L.append(r"\end{figure}")
        L.append(r"\clearpage")

    # ── Section 8b: Per-mu Cell Energy Fraction Profiles (eta_mu only) ──
    if have("cell_data_mu00.pdf"):
        mu_cell_pdfs = ["cell_data", "cell_mc", "cell_mc_m1", "cell_mc_m2"]
        mu_cell_caps = {"cell_data": "Data", "cell_mc": "MC (uncorrected)",
                        "cell_mc_m1": "MC after M1", "cell_mc_m2": "MC after M2"}
        L.append(r"""
%% ====================================================================
\section{Per-$\langle\mu\rangle$ Cell Energy Fraction Profiles}
\label{sec:cell-profiles-mu}
%% ====================================================================

Mean cell energy fraction maps and corrections per $\langle\mu\rangle$ bin.
""")
        for imu, (mulo, muhi) in enumerate(MU_BINS):
            muf = f"mu{imu:02d}"
            avail = [p for p in mu_cell_pdfs if have(f"{p}_{muf}.pdf")]
            if not avail:
                continue
            n_pages_mu = pdf_pages(os.path.join(plots_dir, f"{avail[0]}_{muf}.pdf"))
            sf  = f"cell_shift_{muf}.pdf"
            stf = f"cell_stretch_{muf}.pdf"
            have_corr = have(sf) and have(stf)
            sf_pages  = pdf_pages(os.path.join(plots_dir, sf))  if have_corr else 0
            stf_pages = pdf_pages(os.path.join(plots_dir, stf)) if have_corr else 0
            L.append(rf"\subsection{{{mu_label(mulo, muhi)}}}")
            for i, (lo, hi) in enumerate(ETA_BINS):
                page = i + 1
                if page > n_pages_mu:
                    break
                L.append(rf"\subsubsection*{{{eta_label(lo, hi)}}}")
                fnames = [f"{p}_{muf}.pdf" for p in avail]
                if len(fnames) >= 4:
                    L.append(four_panel(fnames[0], fnames[1], fnames[2], fnames[3],
                                        page,
                                        [mu_cell_caps[p] for p in avail[:4]],
                                        rf"Cell energy fraction maps, {mu_label(mulo, muhi)}, {eta_label(lo, hi)}."))
                else:
                    L.append(r"\begin{figure}[H]")
                    L.append(r"\centering")
                    w = f"{round(0.96/len(fnames), 2):.2f}"
                    for fn, p in zip(fnames, avail):
                        L.append(rf"\begin{{subfigure}}[b]{{{w}\textwidth}}")
                        L.append(rf"  \includegraphics[width=\textwidth,page={page}]{{{fn}}}")
                        L.append(rf"  \caption{{{mu_cell_caps[p]}}}")
                        L.append(r"\end{subfigure}\hfill")
                    L.append(rf"\caption{{Cell energy fraction maps, {mu_label(mulo, muhi)}, {eta_label(lo, hi)}.}}")
                    L.append(r"\end{figure}")
                if have_corr and page <= sf_pages and page <= stf_pages:
                    L.append(two_panel(sf, stf, page,
                                       r"M2 shift ($a_k$)", r"M2 stretch ($s_k$)",
                                       rf"Correction maps, {mu_label(mulo, muhi)}, {eta_label(lo, hi)}."))
            L.append(r"\clearpage")

    # ── Section 8c: Per-pT Cell Energy Fraction Profiles + Corrections (eta_pt only) ──
    if "eta_pt" in variant and have("cell_data_pt00.pdf"):
        pt_cell_pdfs = ["cell_data", "cell_mc", "cell_mc_m1", "cell_mc_m2"]
        pt_cell_caps = {"cell_data": "Data", "cell_mc": "MC (uncorrected)",
                        "cell_mc_m1": "MC after M1", "cell_mc_m2": "MC after M2"}
        L.append(r"""
%% ====================================================================
\section{Per-$p_{\mathrm{T}}$ Cell Energy Fraction Profiles}
\label{sec:cell-profiles-pt}
%% ====================================================================

Mean cell energy fraction maps and corrections per $p_{\mathrm{T}}$ bin.
""")
        for ip, (ptlo, pthi) in enumerate(PT_BINS):
            ptf = f"pt{ip:02d}"
            avail = [p for p in pt_cell_pdfs if have(f"{p}_{ptf}.pdf")]
            if not avail:
                continue
            # Determine max pages from the first available PDF
            n_pages_pt = pdf_pages(os.path.join(plots_dir, f"{avail[0]}_{ptf}.pdf"))
            sf = f"cell_shift_{ptf}.pdf"
            stf = f"cell_stretch_{ptf}.pdf"
            have_corr = have(sf) and have(stf)
            sf_pages = pdf_pages(os.path.join(plots_dir, sf)) if have_corr else 0
            stf_pages = pdf_pages(os.path.join(plots_dir, stf)) if have_corr else 0
            L.append(rf"\subsection{{{pt_label(ptlo, pthi)}}}")
            for i, (lo, hi) in enumerate(ETA_BINS):
                page = i + 1
                if page > n_pages_pt:
                    break
                L.append(rf"\subsubsection*{{{eta_label(lo, hi)}}}")
                fnames = [f"{p}_{ptf}.pdf" for p in avail]
                if len(fnames) >= 4:
                    L.append(four_panel(fnames[0], fnames[1], fnames[2], fnames[3],
                                        page,
                                        [pt_cell_caps[p] for p in avail[:4]],
                                        rf"Cell energy fraction maps, {pt_label(ptlo, pthi)}, {eta_label(lo, hi)}."))
                else:
                    L.append(r"\begin{figure}[H]")
                    L.append(r"\centering")
                    w = f"{round(0.96/len(fnames), 2):.2f}"
                    for fn, p in zip(fnames, avail):
                        L.append(rf"\begin{{subfigure}}[b]{{{w}\textwidth}}")
                        L.append(rf"  \includegraphics[width=\textwidth,page={page}]{{{fn}}}")
                        L.append(rf"  \caption{{{pt_cell_caps[p]}}}")
                        L.append(r"\end{subfigure}\hfill")
                    L.append(rf"\caption{{Cell energy fraction maps, {pt_label(ptlo, pthi)}, {eta_label(lo, hi)}.}}")
                    L.append(r"\end{figure}")
                # Shift + stretch for this pT/eta bin
                if have_corr and page <= sf_pages and page <= stf_pages:
                    L.append(two_panel(sf, stf, page,
                                       r"M2 shift ($a_k$)", r"M2 stretch ($s_k$)",
                                       rf"Correction maps, {pt_label(ptlo, pthi)}, {eta_label(lo, hi)}."))
            L.append(r"\clearpage")

    # ── Section 9: Cell-Level Corrections (eta-only variants only) ──
    if "eta_pt" not in variant and "eta_mu" not in variant:
        shift_pages = pdf_pages(os.path.join(plots_dir, "cell_shift.pdf"))
        stretch_pages = pdf_pages(os.path.join(plots_dir, "cell_stretch.pdf"))
        L.append(r"""
%% ====================================================================
\section{Cell-Level Corrections}
\label{sec:cell-corrections}
%% ====================================================================

M1 shift ($\Delta_k$) and M2 stretch ($s_k$) maps per $|\eta|$ bin.
""")
        if shift_pages >= NBINS and stretch_pages >= NBINS:
            for i, (lo, hi) in enumerate(ETA_BINS):
                page = i + 1
                L.append(rf"\subsection*{{{eta_label(lo, hi)}}}")
                L.append(two_panel("cell_shift.pdf", "cell_stretch.pdf", page,
                                   r"M2 shift ($a_k$)", r"M2 stretch ($s_k$)",
                                   rf"Correction maps, {eta_label(lo, hi)}."))

    L.append(r"""
\end{document}""")

    return "\n".join(L)


# ── compile ────────────────────────────────────────────────────────────────────

def compile_tex(tex_str, out_dir, tex_name):
    """Write tex to out_dir and compile in place (so graphicspath works)."""
    tex_path = os.path.join(out_dir, tex_name)
    pdf_name = tex_name.replace(".tex", ".pdf")

    with open(tex_path, "w") as f:
        f.write(tex_str)

    for run in range(2):
        result = subprocess.run(
            [PDFLATEX, "-interaction=nonstopmode", "-halt-on-error", tex_name],
            cwd=out_dir, capture_output=True, text=True,
        )
        if result.returncode != 0:
            print(f"FAILED (pass {run+1})")
            for line in result.stdout.splitlines()[-20:]:
                print(" ", line)
            return False

    size = os.path.getsize(os.path.join(out_dir, pdf_name)) // 1024
    print(f"OK  ({size} kB)")
    return True


# ── main ─────────────────────────────────────────────────────────────────────

ok, fail = [], []

for variant, variant_desc in VARIANTS.items():
    chi2_tables_path = os.path.join(BASE, variant, "report", "chi2_tables.tex")

    for scenario in SCENARIOS:
        channel = "llgamma"
        sc_dir = os.path.join(BASE, variant, channel, scenario)
        if not os.path.isdir(sc_dir):
            print(f"[skip] {variant}/{channel}/{scenario} not found")
            continue

        tex_name = f"result_compendium_{variant}_{channel}_{scenario}.tex"
        tag = f"{variant}/{scenario}"
        print(f"Building {tag} ...", end=" ", flush=True)

        chi2_table = extract_chi2_table(chi2_tables_path, channel, scenario) if os.path.isfile(chi2_tables_path) else None

        report_dir = os.path.join(BASE, variant, "report")
        plots_dir = os.path.join(sc_dir, "plots")
        tex = build_tex(channel, scenario, variant_desc, chi2_table, plots_dir, report_dir, variant)

        if compile_tex(tex, sc_dir, tex_name):
            ok.append(tag)
        else:
            fail.append(tag)

print(f"\n{len(ok)} PDF(s) written.")
if fail:
    print("FAILED:", ", ".join(fail))
    sys.exit(1)
