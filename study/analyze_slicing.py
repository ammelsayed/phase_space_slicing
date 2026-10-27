#!/usr/bin/env python3
"""Analyze the phase-space slicing study at parton level (LHE files).

For an unsliced reference run and several slice-width variants (each a set
of MG5aMC runs covering misset in slices up to 1000 GeV, last one
unbounded), this script:

  * combines the slices of a variant into a single weighted sample
    (each LHE event keeps its stored weight ~ the sigma of its slice);
  * checks the closure of the combined weighted misset spectrum against
    the reference;
  * compares the per-bin relative statistical uncertainty
    sqrt(sum w^2)/sum w of the combined samples to the reference;
  * quantifies the border artifacts: the per-slice weight discontinuity
    at the slice boundaries, exact-border double counting (the MG5aMC
    misset/missetmax cuts are inclusive on both edges), and the
    additivity of the per-slice cross sections (sum sigma_i vs sigma_ref).

Usage: python3 study/analyze_slicing.py [runs_dir] [outdir]
"""
import glob
import json
import math
import os
import re
import sys

import numpy as np
import pylhe
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INVISIBLE = {12, 14, 16}

VARIANTS = [
    ("ref", "no slicing"),        # runs/ref: single unsliced run
    ("w200", "200 GeV slices"),   # runs/w200: 5 runs
    ("w100", "100 GeV slices"),   # runs/w100: 9 runs
    ("w50", "50 GeV slices"),     # runs/w50: 17 runs + unbounded
]
COLORS = {"ref": "k", "w200": "tab:green", "w100": "tab:blue", "w50": "tab:red"}

BINS = np.arange(150.0, 1101.0, 10.0)
WIDTH_GEVS = [200, 100, 50]


def read_lhe_mets(path):
    """Return (mets, weights) of all events in an LHE file (pylhe)."""
    mets, ws = [], []
    for ev in pylhe.LHEFile.fromfile(path).events:
        px = py = 0.0
        for p in ev.particles:
            if p.status == 1 and abs(p.id) in INVISIBLE:
                px += p.px
                py += p.py
        mets.append(math.hypot(px, py))
        ws.append(ev.eventinfo.weight)
    return np.array(mets), np.array(ws)


def parse_log_sigmas(logpath):
    """Map run name -> (sigma, sigma_err) from an MG5aMC launch log."""
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
    """misset_150_200 -> (150.0, 200.0); misset_1000_inf -> (1000.0, inf)."""
    m = re.match(r"misset_(\d+)_([\d]+|inf)$", name)
    lo = float(m.group(1))
    hi = float("inf") if m.group(2) == "inf" else float(m.group(2))
    return lo, hi


def load_variant(runs_dir, tag):
    """Load all slices of a variant. Returns None if not present yet."""
    if tag == "ref":
        paths = sorted(glob.glob(
            os.path.join(runs_dir, "ref/Events/*/unweighted_events.lhe.gz")))
    else:
        paths = sorted(glob.glob(
            os.path.join(runs_dir, f"{tag}/Events/misset_*/unweighted_events.lhe.gz")))
    if not paths:
        return None
    sigmas = parse_log_sigmas(os.path.join(runs_dir, "LogFiles", f"{tag}.log"))

    slices = []
    for path in paths:
        name = os.path.basename(os.path.dirname(path))
        if name.startswith("misset_"):
            lo, hi = slice_bounds(name)
        else:  # reference run (run_01 or ref): full range
            lo, hi = 150.0, float("inf")
        mets, ws = read_lhe_mets(path)
        sigma, sigma_err = sigmas.get(name, (None, None))
        # MG5aMC stores per-event weights ~ sigma of the run they come
        # from; rescale so that sum(w) = sigma exactly, making the bin
        # contents of different runs directly comparable in pb
        if sigma is not None and ws.sum() > 0:
            ws = ws * (sigma / ws.sum())
        slices.append({"name": name, "lo": lo, "hi": hi, "mets": mets,
                       "ws": ws, "sigma": sigma, "sigma_err": sigma_err,
                       "path": path})
    return {"slices": slices}


def combined(sample):
    mets = np.concatenate([s["mets"] for s in sample["slices"]])
    ws = np.concatenate([s["ws"] for s in sample["slices"]])
    return mets, ws


def hist_relerr(x, w, bins):
    h, _ = np.histogram(x, bins=bins, weights=w)
    h2, _ = np.histogram(x, bins=bins, weights=w * w)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(h > 0, np.sqrt(h2) / np.abs(h), np.nan)
    return h, rel


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
        sys.exit("reference run not found")

    ref_mets, ref_ws = combined(loaded["ref"])
    ref_h, ref_rel = hist_relerr(ref_mets, ref_ws, BINS)
    centers = 0.5 * (BINS[:-1] + BINS[1:])

    summary = {"bins": BINS.tolist(), "variants": {}}

    # ---------- per-variant figures ----------
    for tag, label in VARIANTS:
        if tag not in loaded or tag == "ref":
            continue
        s = loaded[tag]
        mets, ws = combined(s)
        h, rel = hist_relerr(mets, ws, BINS)

        # sigma additivity
        sig = np.array([sl["sigma"] for sl in s["slices"] if sl["sigma"] is not None])
        sig_e = np.array([sl["sigma_err"] for sl in s["slices"] if sl["sigma"] is not None])
        ref_sigma = loaded["ref"]["slices"][0]["sigma"]
        ref_sigma_e = loaded["ref"]["slices"][0]["sigma_err"]
        sigma_sum = float(sig.sum())
        sigma_sum_e = float(np.sqrt((sig_e ** 2).sum()))
        pull = ((sigma_sum - ref_sigma)
                / math.sqrt(sigma_sum_e ** 2 + ref_sigma_e ** 2)
                if ref_sigma else float("nan"))

        # border checks
        borders = []
        for a, b in zip(s["slices"][:-1], s["slices"][1:]):
            bnd = a["hi"]
            below = a["ws"][a["mets"] >= bnd - 10.0]
            above = b["ws"][b["mets"] <= bnd + 10.0]
            mean_below = float(below.mean()) if below.size else float("nan")
            mean_above = float(above.mean()) if above.size else float("nan")
            # exact-border double counting (inclusive edges)
            tol = 1e-6
            n_lo = int((np.abs(a["mets"] - bnd) < tol).sum())
            n_hi = int((np.abs(b["mets"] - bnd) < tol).sum())
            borders.append({"border": bnd, "mean_w_below": mean_below,
                            "mean_w_above": mean_above,
                            "w_jump": mean_above / mean_below if mean_below else None,
                            "n_at_border_below": n_lo, "n_at_border_above": n_hi})
        # slice bound validation
        for sl in s["slices"]:
            if sl["mets"].size:
                if sl["mets"].min() < sl["lo"] - 1e-6 or (
                        math.isfinite(sl["hi"]) and sl["mets"].max() > sl["hi"] + 1e-6):
                    print(f"[WARN] {sl['name']}: met outside "
                          f"[{sl['lo']:g}, {sl['hi']:g}]")

        # closure of each slice against the reference spectrum:
        # sigma_i (MG5 integration) vs the reference content in [lo, hi]
        ref_mask_sum = float(ref_ws.sum())
        closure = []
        for sl in s["slices"]:
            if sl["sigma"] is None or not math.isfinite(sl["hi"]):
                continue
            m = (ref_mets >= sl["lo"]) & (ref_mets < sl["hi"])
            ref_content = float(ref_ws[m].sum())  # already in pb
            n_ref = int(m.sum())
            ref_stat = ref_content * math.sqrt(max(1, n_ref)) / n_ref if n_ref else float("nan")
            c_pull = ((sl["sigma"] - ref_content)
                      / math.sqrt(sl["sigma_err"] ** 2 + ref_stat ** 2)
                      if n_ref else float("nan"))
            closure.append({"slice": sl["name"], "sigma_i": sl["sigma"],
                            "sigma_i_err": sl["sigma_err"],
                            "ref_content": ref_content,
                            "pull": c_pull})
        worst_pull = max((abs(c["pull"]) for c in closure), default=0)

        neff = float(ws.sum() ** 2 / (ws ** 2).sum())
        summary["variants"][tag] = {
            "label": label, "n_runs": len(s["slices"]),
            "n_events_total": int(mets.size),
            "sigma_sum_pb": sigma_sum, "sigma_sum_err_pb": sigma_sum_e,
            "sigma_ref_pb": ref_sigma, "sigma_ref_err_pb": ref_sigma_e,
            "sigma_additivity_pull": pull,
            "sum_w": float(ws.sum()), "n_eff": neff,
            "ref_n_eff": float(ref_ws.sum() ** 2 / (ref_ws ** 2).sum()),
            "worst_closure_pull": worst_pull,
            "closure": closure,
            "slices": [{"name": sl["name"], "lo": sl["lo"], "hi": sl["hi"],
                        "n": int(sl["mets"].size), "sigma": sl["sigma"],
                        "sigma_err": sl["sigma_err"],
                        "mean_w": float(sl["ws"].mean()) if sl["ws"].size else None}
                       for sl in s["slices"]],
            "borders": borders,
        }

        # closure + rel-err figure
        fig, (a1, a2, a3) = plt.subplots(
            3, 1, sharex=True, figsize=(8, 9),
            gridspec_kw={"height_ratios": [2.2, 1.3, 1.3]})
        a1.step(centers, ref_h, where="mid", color="k", label="no slicing (ref)")
        a1.step(centers, h, where="mid", color=COLORS[tag], label=label)
        a1.set_yscale("log")
        a1.set_ylabel(r"$\Sigma w$ per bin [pb]")
        a1.legend(fontsize=9)
        a1.set_title(f"Slicing study, parton level: {label}")
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(ref_h > 0, h / ref_h, np.nan)
        a2.errorbar(centers, ratio, yerr=np.where(ref_h > 0, ref_rel * ratio, np.nan),
                    fmt="o", ms=3, color=COLORS[tag])
        a2.axhline(1.0, color="gray", lw=0.8)
        a2.set_ylim(0.7, 1.3)
        a2.set_ylabel("combined / ref (both in pb)")
        for b in [sl["hi"] for sl in s["slices"][:-1]]:
            a2.axvline(b, color="gray", lw=0.4, alpha=0.6)
        a3.step(centers, ref_rel, where="mid", color="k", label="no slicing (ref)")
        a3.step(centers, rel, where="mid", color=COLORS[tag], label=label)
        a3.set_yscale("log")
        a3.set_ylabel(r"$\delta = \sqrt{\Sigma w^2}/\Sigma w$ per bin")
        a3.set_xlabel("misset [GeV]")
        a3.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, f"closure_{tag}.png"), dpi=130)
        plt.close(fig)

    # ---------- summary figures ----------
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.step(centers, ref_rel, where="mid", color="k", label="no slicing (ref)")
    for tag, lbl in VARIANTS:
        if tag == "ref" or tag not in loaded:
            continue
        mets, ws = combined(loaded[tag])
        _, rel = hist_relerr(mets, ws, BINS)
        ax.step(centers, rel, where="mid", color=COLORS[tag], label=lbl)
    ax.set_yscale("log")
    ax.set_xlabel("misset [GeV]")
    ax.set_ylabel(r"$\delta = \sqrt{\Sigma w^2}/\Sigma w$ per bin")
    ax.set_title("Per-bin relative statistical uncertainty (parton level)")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "relerr_profiles.png"), dpi=130)
    plt.close(fig)

    # mean event weight vs met: the per-slice weight steps at the borders
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
    fig.savefig(os.path.join(outdir, "weight_steps.png"), dpi=130)
    plt.close(fig)

    with open(os.path.join(outdir, "slicing_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)

    # ---------- markdown report ----------
    print("# Slicing study (parton level)\n")
    ref_info = summary["variants"].get("ref") or {}
    for tag, d in summary["variants"].items():
        if tag == "ref":
            continue
        print(f"## {d['label']} ({d['n_runs']} runs, "
              f"{d['n_events_total']} events)\n")
        print(f"- sigma additivity: sum(sigma_i) = {d['sigma_sum_pb']:.4f} "
              f"+- {d['sigma_sum_err_pb']:.4f} pb  vs  sigma_ref = "
              f"{d['sigma_ref_pb']:.4f} +- {d['sigma_ref_err_pb']:.4f} pb "
              f"(pull = {d['sigma_additivity_pull']:+.2f})")
        print(f"- effective sample size: N_eff = {d['n_eff']:.3g} "
              f"(reference {d['ref_n_eff']:.3g})")
        print(f"- per-slice closure vs reference spectrum: worst |pull| = "
              f"{d['worst_closure_pull']:.2f}")
        print(f"  (sigma_i from the slice integration vs reference content "
              f"in the same misset range)")
        worst = max((abs(math.log(b["w_jump"])) for b in d["borders"]
                     if b["w_jump"]), default=0)
        print(f"- largest per-slice weight step at a border: "
              f"x{math.exp(worst):.2f}")
        ndc = sum(b["n_at_border_below"] + b["n_at_border_above"]
                  for b in d["borders"])
        print(f"- events sitting exactly on a border value: {ndc}")
        print()
        print("| slice | N | sigma [pb] | mean weight [pb] |")
        print("|---|---|---|---|")
        for sl in d["slices"]:
            hi = "inf" if math.isinf(sl["hi"]) else f"{sl['hi']:g}"
            print(f"| {sl['name']} [{sl['lo']:g}, {hi}] | {sl['n']} | "
                  f"{sl['sigma']:.4f} +- {sl['sigma_err']:.4f} | "
                  f"{sl['mean_w']:.4f} |")
        print()


if __name__ == "__main__":
    main()
