#!/usr/bin/env python3
"""Analyze the phase-space slicing study at reconstruction level.

Reads the Delphes ROOT files of the slicing runs (produced with
--shower Pythia8 --detector Delphes by study/run_generation.sh) and
repeats the LHE-level comparison at reconstruction level, for the two
key observables:

  * MET  : from the MissingET branch (MissingET.MET[0]);
  * ptl1 : pt of the leading lepton, from the Electron and Muon
           branches (highest Electron.PT / Muon.PT of the event);

per slice width, compared to the unsliced reference:

  * combined weighted spectra (closure);
  * per-bin relative statistical uncertainty sqrt(sum w^2)/sum w;
  * the per-slice weight steps at the borders, crossed by the detector
    resolution: events generated in one slice reconstruct in the
    neighbouring one (migration), quantified with a parton-MET vs
    reco-MET map;
  * sigma additivity of the slices.

Environment: ROOT must be sourced
(source /home/ammelsayed/root/bin/thisroot.sh); the Delphes libraries of
the MG5aMC installation that produced the files are loaded exactly like
in typeIIIseesaw_phenom/samples/make_dataset.py.

Usage:
  source /home/ammelsayed/root/bin/thisroot.sh
  python3 study/analyze_delphes.py [runs_dir] [outdir]
"""
import glob
import json
import math
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import ROOT

DELPHES_PATH = os.environ.get(
    "DELPHES_HOME", "/home/ammelsayed/softwares/MG5_aMC_v3_7_2/Delphes")
ROOT.gInterpreter.AddIncludePath(DELPHES_PATH)
ROOT.gInterpreter.AddIncludePath(f"{DELPHES_PATH}/classes")
ROOT.gInterpreter.AddIncludePath(f"{DELPHES_PATH}/external")
if not os.path.exists(os.path.join(DELPHES_PATH, "libDelphes.so")):
    sys.exit(f"libDelphes.so not found in {DELPHES_PATH}")
# load by full path: LD_LIBRARY_PATH does not contain the MG5 Delphes dir
if ROOT.gSystem.Load(os.path.join(DELPHES_PATH, "libDelphes.so")) != 0:
    sys.exit(f"could not load libDelphes from {DELPHES_PATH}")
ROOT.gROOT.SetBatch(True)

VARIANTS = [
    ("ref", "no slicing"),
    ("w200", "200 GeV slices"),
    ("w100", "100 GeV slices"),
    ("w50", "50 GeV slices"),
]
COLORS = {"ref": "k", "w200": "tab:green", "w100": "tab:blue", "w50": "tab:red"}

# observable name -> (histogram bins, axis label)
OBS = {
    "met": (np.arange(0.0, 1101.0, 10.0), "reco MET [GeV]"),
    "ptl1": (np.arange(100.0, 1101.0, 10.0), "reco $p_T$(leading lepton) [GeV]"),
}


def parse_log_sigmas(logpath):
    if not os.path.exists(logpath):
        return {}
    text = open(logpath, errors="ignore").read()
    out = {}
    for m in re.finditer(
            r"Results Summary for run: (\S+) tag.*?"
            r"Cross-section :\s+([0-9.eE+-]+) \+- ([0-9.eE+-]+)",
            text, re.S):
        out[m.group(1)] = (float(m.group(2)), float(m.group(3)))
    return out


def slice_bounds(name):
    m = re.match(r"misset_(\d+)_(\d+|inf)$", name)
    if not m:
        return 0.0, float("inf")
    return float(m.group(1)), (float("inf") if m.group(2) == "inf"
                               else float(m.group(2)))


def read_delphes(path):
    """Return (met_reco, ptl1_reco, weights, met_parton) of one ROOT file.

    met comes from the MissingET branch, ptl1 from the Electron and Muon
    branches, weights from Event.Weight. met_parton is the transverse
    momentum of the vector sum of the stable invisible particles of the
    hard event (the parton-level MET the slices were cut on).
    """
    rdf = ROOT.RDataFrame("Delphes", path)
    arrays = rdf.Define("met", "MissingET.MET[0]") \
                .Define("wgt", "Event.Weight") \
                .Define("ptl1",
                        "max(Electron.PT.size() > 0 ? Max(Electron.PT) : -1.0,"
                        "     Muon.PT.size() > 0 ? Max(Muon.PT) : -1.0)") \
                .AsNumpy(["met", "wgt", "ptl1"])
    met_reco = np.asarray(arrays["met"], dtype=float)
    ptl1_reco = np.asarray(arrays["ptl1"], dtype=float)
    # Event.Weight arrives as one RVec per event (size 1 for these files)
    ws = np.array([float(np.asarray(w).ravel()[0]) for w in arrays["wgt"]],
                  dtype=float)

    # parton-level MET from the stable invisible particles
    met_parton = np.zeros_like(met_reco)
    chain = ROOT.TChain("Delphes")
    chain.Add(path)
    part = ROOT.TClonesArray("GenParticle")
    chain.SetBranchAddress("Particle", ROOT.AddressOf(part))
    nev = min(chain.GetEntries(), met_reco.size)
    for i in range(nev):
        chain.GetEntry(i, 0)
        px = py = 0.0
        for p in part:
            if p.Status != 1:
                continue
            if abs(p.PID) in (12, 14, 16):
                px += p.Px
                py += p.Py
        met_parton[i] = math.hypot(px, py)
    return met_reco, ptl1_reco, ws, met_parton


def hist_relerr(x, w, bins):
    h, _ = np.histogram(x, bins=bins, weights=w)
    h2, _ = np.histogram(x, bins=bins, weights=w * w)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(h > 0, np.sqrt(h2) / np.abs(h), np.nan)
    return h, rel


def load_variant(runs_dir, tag):
    if tag == "ref":
        paths = sorted(glob.glob(
            os.path.join(runs_dir, "ref/Events/*/tag_1_delphes_events.root")))
    else:
        paths = sorted(glob.glob(
            os.path.join(runs_dir, f"{tag}/Events/misset_*/tag_1_delphes_events.root")))
    if not paths:
        return None
    sigmas = parse_log_sigmas(os.path.join(runs_dir, "LogFiles", f"{tag}.log"))
    slices = []
    for path in paths:
        name = os.path.basename(os.path.dirname(path))
        lo, hi = slice_bounds(name)
        met_reco, ptl1_reco, ws, met_parton = read_delphes(path)
        sigma, sigma_err = sigmas.get(name, (None, None))
        # rescale so sum(w) = sigma, making runs directly comparable in pb
        if sigma is not None and ws.sum() > 0:
            ws = ws * (sigma / ws.sum())
        slices.append({"name": name, "lo": lo, "hi": hi, "met_reco": met_reco,
                       "ptl1_reco": ptl1_reco, "met_parton": met_parton,
                       "ws": ws, "sigma": sigma, "sigma_err": sigma_err,
                       "path": path})
    return {"slices": slices}


def main():
    runs_dir = sys.argv[1] if len(sys.argv) > 1 else "runs"
    outdir = sys.argv[2] if len(sys.argv) > 2 else "study/plots"
    os.makedirs(outdir, exist_ok=True)

    loaded = {}
    for tag, label in VARIANTS:
        s = load_variant(runs_dir, tag)
        if s is None:
            print(f"[WARN] variant {tag} not found, skipped")
            continue
        loaded[tag] = s
    if "ref" not in loaded:
        sys.exit("Delphes reference run not found")

    ref_sl = loaded["ref"]["slices"][0]
    ref_sigma = ref_sl["sigma"]
    ref_sigma_e = ref_sl["sigma_err"]
    ref_mets = ref_sl["met_reco"]
    ref_ws = ref_sl["ws"]

    summary = {"variants": {}}
    for tag, label in VARIANTS:
        if tag not in loaded or tag == "ref":
            continue
        s = loaded[tag]
        ws = np.concatenate([sl["ws"] for sl in s["slices"]])

        sig = np.array([sl["sigma"] for sl in s["slices"] if sl["sigma"] is not None])
        sig_e = np.array([sl["sigma_err"] for sl in s["slices"] if sl["sigma"] is not None])
        sigma_sum = float(sig.sum())
        sigma_sum_e = float(np.sqrt((sig_e ** 2).sum()))
        pull = ((sigma_sum - ref_sigma)
                / math.sqrt(sigma_sum_e ** 2 + ref_sigma_e ** 2)
                if ref_sigma else float("nan"))

        borders = []
        for a, b in zip(s["slices"][:-1], s["slices"][1:]):
            bnd = a["hi"]
            below = a["ws"][(a["met_reco"] >= bnd - 10.0) & (a["met_reco"] < bnd + 10.0)]
            above = b["ws"][(b["met_reco"] >= bnd - 10.0) & (b["met_reco"] < bnd + 10.0)]
            borders.append({
                "border": bnd,
                "mean_w_below": float(below.mean()) if below.size else None,
                "mean_w_above": float(above.mean()) if above.size else None,
                "n_reco_above_border_from_below_slice": int((a["met_reco"] >= bnd).sum()),
                "n_reco_below_border_from_above_slice": int((b["met_reco"] < bnd).sum()),
            })

        summary["variants"][tag] = {
            "label": label, "n_runs": len(s["slices"]),
            "n_events_total": int(ws.size),
            "sigma_sum_pb": sigma_sum, "sigma_sum_err_pb": sigma_sum_e,
            "sigma_ref_pb": ref_sigma, "sigma_ref_err_pb": ref_sigma_e,
            "sigma_additivity_pull": pull,
            "n_eff": float(ws.sum() ** 2 / (ws ** 2).sum()),
            "borders": borders,
        }

        # per-observable closure + rel err figures
        for obs, (bins, xlabel) in OBS.items():
            key = f"{obs}_reco"
            ref_vals = ref_sl[key]
            ref_h, ref_rel = hist_relerr(ref_vals, ref_ws, bins)
            vals = np.concatenate([sl[key] for sl in s["slices"]])
            h, rel = hist_relerr(vals, ws, bins)
            centers = 0.5 * (bins[:-1] + bins[1:])
            summary["variants"][tag][f"closure_{obs}"] = {
                "worst_abs_log_ratio": float(np.nanmax(np.abs(np.log(
                    np.where((h > 0) & (ref_h > 0), h / ref_h, np.nan)[
                        (centers > 250) & (centers < 900)])))) if np.any(
                    (h > 0) & (ref_h > 0)) else None,
            }

            fig, (a1, a2, a3) = plt.subplots(
                3, 1, sharex=True, figsize=(8, 9),
                gridspec_kw={"height_ratios": [2.2, 1.3, 1.3]})
            a1.step(centers, ref_h, where="mid", color="k",
                    label="no slicing (ref)")
            a1.step(centers, h, where="mid", color=COLORS[tag], label=label)
            a1.set_yscale("log")
            a1.set_ylabel(r"$\Sigma w$ per bin [pb]")
            a1.legend(fontsize=9)
            a1.set_title(f"Slicing study, reconstruction level: {label}")
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where(ref_h > 0, h / ref_h, np.nan)
            a2.errorbar(centers, ratio,
                        yerr=np.where(ref_h > 0, ref_rel * ratio, np.nan),
                        fmt="o", ms=3, color=COLORS[tag])
            a2.axhline(1.0, color="gray", lw=0.8)
            a2.set_ylim(0.7, 1.3)
            a2.set_ylabel("combined / ref (both in pb)")
            for b in [sl["hi"] for sl in s["slices"][:-1]]:
                a2.axvline(b, color="gray", lw=0.4, alpha=0.6)
            a3.step(centers, ref_rel, where="mid", color="k",
                    label="no slicing (ref)")
            a3.step(centers, rel, where="mid", color=COLORS[tag], label=label)
            a3.set_yscale("log")
            a3.set_ylabel(r"$\delta = \sqrt{\Sigma w^2}/\Sigma w$ per bin")
            a3.set_xlabel(xlabel)
            a3.legend(fontsize=9)
            fig.tight_layout()
            fig.savefig(os.path.join(outdir, f"closure_{tag}_{obs}.png"), dpi=130)
            plt.close(fig)

    # rel-err profile summaries per observable
    for obs, (bins, xlabel) in OBS.items():
        key = f"{obs}_reco"
        ref_h, ref_rel = hist_relerr(ref_sl[key], ref_ws, bins)
        centers = 0.5 * (bins[:-1] + bins[1:])
        fig, ax = plt.subplots(figsize=(8, 5.5))
        ax.step(centers, ref_rel, where="mid", color="k", label="no slicing (ref)")
        for tag, lbl in VARIANTS:
            if tag == "ref" or tag not in loaded:
                continue
            vals = np.concatenate([sl[key] for sl in loaded[tag]["slices"]])
            ws = np.concatenate([sl["ws"] for sl in loaded[tag]["slices"]])
            _, rel = hist_relerr(vals, ws, bins)
            ax.step(centers, rel, where="mid", color=COLORS[tag], label=lbl)
        ax.set_yscale("log")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$\delta = \sqrt{\Sigma w^2}/\Sigma w$ per bin")
        ax.set_title(f"Per-bin relative statistical uncertainty "
                     f"(reconstruction level, {obs})")
        ax.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, f"relerr_profiles_delphes_{obs}.png"),
                    dpi=130)
        plt.close(fig)

    # per-slice weight steps (same weights for all observables)
    centers = 0.5 * (OBS["met"][0][:-1] + OBS["met"][0][1:])
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.axhline(ref_ws.mean(), color="k", lw=1.2, label="no slicing (ref)")
    for tag, lbl in VARIANTS:
        if tag == "ref" or tag not in loaded:
            continue
        for sl in loaded[tag]["slices"]:
            if not sl["ws"].size:
                continue
            hi = sl["hi"] if math.isfinite(sl["hi"]) else 1100.0
            ax.plot([max(sl["lo"], 150.0), hi], [sl["ws"].mean()] * 2,
                    color=COLORS[tag], lw=1.6)
        ax.plot([], [], color=COLORS[tag], label=lbl)
    ax.set_yscale("log")
    ax.set_xlabel("misset [GeV]")
    ax.set_ylabel("mean event weight [pb] (normalized: $\\Sigma w=\\sigma$)")
    ax.set_title("Per-slice event weight: steps at the slice borders")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "weight_steps_delphes.png"), dpi=130)
    plt.close(fig)

    # stacked per-slice contributions to the reco MET spectrum (w50):
    # shows the migration across borders directly
    if "w50" in loaded:
        s = loaded["w50"]
        bins, xlabel = OBS["met"]
        fig, ax = plt.subplots(figsize=(8, 5.5))
        colors = plt.cm.viridis(np.linspace(0, 0.95, len(s["slices"])))
        bottoms = np.zeros(len(bins) - 1)
        for sl, c in zip(s["slices"], colors):
            hh, _ = np.histogram(sl["met_reco"], bins=bins, weights=sl["ws"])
            ax.bar(centers, hh, bottom=bottoms, width=np.diff(bins),
                   color=c, alpha=0.85, label=sl["name"])
            bottoms += hh
        hh, _ = np.histogram(ref_mets, bins=bins, weights=ref_ws)
        ax.step(centers, hh, where="mid", color="k", lw=1.2,
                label="no slicing (ref)")
        for b in [sl["hi"] for sl in s["slices"][:-1]]:
            ax.axvline(b, color="gray", lw=0.5, ls=":")
        ax.set_yscale("log")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$\Sigma w$ per bin [pb]")
        ax.set_title("50 GeV slices: per-slice contributions to the reco MET "
                     "spectrum (migration across parton borders)")
        ax.legend(fontsize=6, ncol=2)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "migration_w50.png"), dpi=130)
        plt.close(fig)

        # parton vs reco MET map
        pm = np.concatenate([sl["met_parton"] for sl in s["slices"]])
        rm = np.concatenate([sl["met_reco"] for sl in s["slices"]])
        step = max(1, pm.size // 200000)
        fig, ax = plt.subplots(figsize=(6.5, 6))
        ax.hist2d(pm[::step], rm[::step], bins=[np.arange(100, 1101, 25),
                                                np.arange(0, 1101, 25)],
                   norm=matplotlib.colors.LogNorm(), cmap="viridis")
        for b in [sl["hi"] for sl in s["slices"][:-1]]:
            ax.axvline(b, color="red", lw=0.5, alpha=0.7)
            ax.axhline(b, color="red", lw=0.5, alpha=0.7)
        ax.plot([0, 1100], [0, 1100], color="white", lw=0.8, alpha=0.6)
        ax.set_xlabel("parton-level MET [GeV] (slice borders in red)")
        ax.set_ylabel("reco MET [GeV]")
        ax.set_title("50 GeV slices: parton vs reco MET")
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "parton_vs_reco_w50.png"), dpi=130)
        plt.close(fig)

    with open(os.path.join(outdir, "delphes_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)

    print("# Slicing study (reconstruction level)\n")
    for tag, d in summary["variants"].items():
        print(f"## {d['label']} ({d['n_runs']} runs, {d['n_events_total']} events)\n")
        print(f"- sigma additivity: sum(sigma_i) = {d['sigma_sum_pb']:.4f} "
              f"+- {d['sigma_sum_err_pb']:.4f} pb vs sigma_ref = "
              f"{d['sigma_ref_pb']:.4f} +- {d['sigma_ref_err_pb']:.4f} pb "
              f"(pull = {d['sigma_additivity_pull']:+.2f})")
        print(f"- N_eff = {d['n_eff']:.3g}")
        mig_up = sum(b["n_reco_above_border_from_below_slice"] for b in d["borders"])
        mig_dn = sum(b["n_reco_below_border_from_above_slice"] for b in d["borders"])
        print(f"- reco events migrated up across a border: {mig_up}; "
              f"down: {mig_dn}")
        for obs in OBS:
            c = d.get(f"closure_{obs}", {}).get("worst_abs_log_ratio")
            if c is not None:
                print(f"- worst |ln(combined/ref)| for {obs} in 250-900 GeV: {c:.3f}")
        print()


if __name__ == "__main__":
    main()
