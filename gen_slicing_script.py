#!/usr/bin/env python3
"""Generate the MG5aMC phase-space slicing script for the W -> l nu sample.

The base script is the w012j_MLM.mg5 of the phase_space_biasing repository
(generate p p > lep nu QCD=0 QED=2, ptl1min 150, 13.6 TeV) with:
  * 'set ptj 30.0' and 'set etaj 2.5' removed (no jets in the final state),
  * the '@0' multiplicity label dropped from the generate line,
  * one named 'launch' per slice of the slicing variable:
        launch w012j_MLM --name=misset_{lo}_{hi}
        set <min_cut> {lo}
        set <max_cut> {hi}
        done
    so that every slice is an independent MG5aMC run with its own cross
    section, written to Events/<name>/.

Run-card cut names verified in MG5aMC v3.7.2: the MET minimum cut is
'misset' (there is NO 'missetmin') and the maximum cut is 'missetmax'
(negative = no cut). Both edges are inclusive in the generated cuts.f.
The LAST slice is always unbounded on the upper edge ('missetmax = -1.0',
named e.g. misset_450_inf) so that no phase space is discarded.

Examples:
    python3 gen_slicing_script.py                          # misset 150-450 + [450,inf], 50 GeV slices
    python3 gen_slicing_script.py --width 25 --max 1000    # finer slices, wider range
    python3 gen_slicing_script.py --detector Delphes       # reconstruction-level stage

Note: the default lower edge is 150 GeV (= ptl1min): at LO the lepton and
the neutrino are back-to-back, so MET = ptl1 >= ptl1min for every event
and slices below 150 GeV contain no phase space (MG5aMC reports a zero
result for them).
"""
import argparse

# run_card parameter names of the slicing variable (min, max), verified
# for MG5aMC v3.7.2 (Template/LO/Cards/run_card.dat). Extend the dict for
# other variables, or override with --min-param/--max-param.
CUT_PARAMS = {
    "misset": ("misset", "missetmax"),  # min is 'misset', NOT 'missetmin'!
}


def launch_block(proc, name, shower, detector, analysis, first,
                 min_param, max_param, lo, hi, nevents, unbounded=False):
    lines = [
        f"launch {proc} --name={name}",
        f"shower = {shower}",
        f"detector = {detector}",
        f"analysis = {analysis}",
        "madspin = OFF",
        "reweight = OFF",
        "done",
    ]
    if first:
        # full run-card definition; identical values are re-used by the
        # later launch blocks (the process run_card persists)
        lines += [
            "set ebeam 6800",
            "# shower and detector are OFF: pythia8/delphes card lines commented out",
            "# ./Cards_edited/pythia8_card_TuneA14.dat",
            "# ./Cards_edited/delphes_card_ATLAS_FatJobPTMin50.tcl",
            "set pdlabel lhapdf",
            "set lhaid 331900",
            "# set ickkw 1",
            "# set xqcut 30.0",
            "# set maxjetflavor 5",
            "set use_syst False",
            "set cut_decays True",
            "set ptl1min 150.0",
            "set etal 2.5",
            "set sde_strategy 2",
            "set hard_survey 2",
        ]
    # every block is explicit about what defines the run: slice bounds,
    # seed and number of events
    if unbounded:
        max_line = f"set {max_param} -1.0  # last slice: no upper cut"
    else:
        max_line = f"set {max_param} {hi:.1f}"
    lines += [
        f"set {min_param} {lo:.1f}",
        max_line,
        "set iseed 0",
        f"set nevents {nevents}",
        "done",
        "",
    ]
    return "\n".join(lines)


def main():
    p = argparse.ArgumentParser(
        description="Generate the MG5aMC slicing script (one named launch "
                    "per slice of the slicing variable).")
    p.add_argument("--variable", default="misset",
                   help="slicing variable (default: misset)")
    p.add_argument("--min", type=float, default=150.0,
                   help="lower edge of the sliced range (default: 150.0 = "
                        "ptl1min; slices below it are empty because "
                        "MET = ptl1 >= ptl1min at LO)")
    p.add_argument("--max", type=float, default=500.0,
                   help="upper edge of the sliced range (default: 500.0)")
    p.add_argument("--width", type=float, default=50.0,
                   help="slice width in GeV (default: 50.0)")
    p.add_argument("--nevents", default="100_000",
                   help="events per slice (default: 100_000)")
    p.add_argument("--name-format", default="{var}_{lo:g}_{hi:g}",
                   help="run naming style (default: {var}_{lo:g}_{hi:g})")
    p.add_argument("--process", default="w012j_MLM",
                   help="MG5aMC process/output name (default: w012j_MLM)")
    p.add_argument("--output", default=None,
                   help="output script path (default: <process>.mg5)")
    p.add_argument("--shower", default="OFF",
                   help="shower switch (default: OFF)")
    p.add_argument("--detector", default="OFF",
                   help="detector switch (default: OFF; use Delphes for the "
                        "reconstruction-level stage)")
    p.add_argument("--analysis", default="MadAnalysis5",
                   help="analysis switch (default: MadAnalysis5)")
    p.add_argument("--min-param", default=None,
                   help="run_card name of the minimum cut "
                        "(default: from CUT_PARAMS)")
    p.add_argument("--max-param", default=None,
                   help="run_card name of the maximum cut "
                        "(default: from CUT_PARAMS)")
    args = p.parse_args()

    if args.variable not in CUT_PARAMS and (not args.min_param or
                                            not args.max_param):
        p.error(f"unknown variable '{args.variable}': pass --min-param and "
                f"--max-param explicitly (known: {list(CUT_PARAMS)})")
    min_param, max_param = CUT_PARAMS[args.variable]
    if args.min_param:
        min_param = args.min_param
    if args.max_param:
        max_param = args.max_param

    lo_min, hi_max, width = args.min, args.max, args.width
    if width <= 0:
        p.error("--width must be positive")
    if hi_max <= lo_min:
        p.error("--max must be greater than --min")

    n_slices = (hi_max - lo_min) / width
    n_full = round(n_slices)
    if abs(n_slices - n_full) > 1e-9:
        print(f"WARNING: range [{lo_min:g}, {hi_max:g}] is not a multiple of "
              f"the width {width:g}; the last slice is narrower.")

    blocks = []
    edges = [lo_min + i * width for i in range(n_full + 1)]
    if edges[-1] < hi_max - 1e-9:
        edges.append(hi_max)
    # the last slice is always unbounded on the upper edge: no phase space
    # may be discarded (missetmax = -1 disables the max cut in MG5aMC)
    pairs = list(zip(edges[:-1], edges[1:]))
    for i, (lo, hi) in enumerate(pairs):
        last = i == len(pairs) - 1
        name = args.name_format.format(
            var=args.variable, lo=lo, hi=float("inf") if last else hi)
        blocks.append(launch_block(args.process, name, args.shower,
                                   args.detector, args.analysis, i == 0,
                                   min_param, max_param, lo, hi, args.nevents,
                                   unbounded=last))

    header = f"""\
set automatic_html_opening False --no_save
set notification_center False --no_save
set zerowidth_tchannel False --no_save
set nb_core 40 --no_save
set nb_core_pythia8 40 --no_save
set nb_core_delphes 40 --no_save
import model sm
define p = g u c d s b u~ c~ d~ s~ b~
define j = g u c d s b u~ c~ d~ s~ b~
define l+ = e+ mu+ ta+
define l- = e- mu- ta-
define vl = ve vm vt
define vl~ = ve~ vm~ vt~
define lep = l+ l-
define nu = vl vl~
generate p p > lep nu QCD=0 QED=2
output {args.process} --f
"""

    script = header + "\n".join(blocks)
    output = args.output or f"{args.process}.mg5"
    with open(output, "w") as f:
        f.write(script)

    def edge_str(hi):
        return "inf" if hi == float("inf") else f"{hi:g}"

    print(f"wrote {output}: {len(blocks)} slices of {args.variable} "
          f"({edges[0]:g} to {edge_str(edges[-1])}, width {width:g}, "
          f"{args.nevents} events each; last slice unbounded)")
    for i, (lo, hi) in enumerate(pairs):
        print(f"  [{i}] {args.name_format.format(var=args.variable, lo=lo, hi=float('inf') if i == len(pairs) - 1 else hi)}"
              f"  ->  [{lo:g}, {edge_str(hi)}]")


if __name__ == "__main__":
    main()
