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
    python3 gen_slicing_script.py --ref                    # inclusive reference (no MET cut)
    python3 gen_slicing_script.py --detector Delphes       # with reconstruction

Note: the default lower edge is 150 GeV (= ptl1min): at LO the lepton and
the neutrino are back-to-back, so MET = ptl1 >= ptl1min for every event
and slices below 150 GeV contain no phase space (MG5aMC reports a zero
result for them).

With --pythia8-card/--delphes-card the (uncommented) card paths are
written into every launch block; they must be valid relative to the
directory the generated script is run from (runmg5.sh runs it from its
own directory).
"""
import argparse

# run_card parameter names of the slicing variable (min, max), verified
# for MG5aMC v3.7.2 (Template/LO/Cards/run_card.dat). Extend the dict for
# other variables, or override with --min-param/--max-param.
CUT_PARAMS = {
    "misset": ("misset", "missetmax"),  # min is 'misset', NOT 'missetmin'!
}


def launch_block(proc, name, shower, detector, analysis, first,
                 min_param, max_param, lo, hi, nevents, unbounded=False,
                 pythia8_card=None, delphes_card=None, bias_module=None,
                 bias_parameters=None):
    lines = [
        f"launch {proc} --name={name}",
        f"shower = {shower}",
        f"detector = {detector}",
        f"analysis = {analysis}",
        "madspin = OFF",
        "reweight = OFF",
        "done",
    ]
    # mirror the original working scripts: the first run_card line after
    # the switch stage is 'set ebeam', then the card paths, then the rest
    lines.append("set ebeam 6800")
    if bias_module:
        lines.append(f"set bias_module {bias_module}")
    if bias_parameters:
        lines.append(f"set bias_parameters {bias_parameters}")
    if pythia8_card:
        lines.append(pythia8_card)
    if delphes_card:
        lines.append(delphes_card)
    if first:
        lines += [
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
        max_line = f"set {max_param} -1.0  # no upper cut"
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
    p.add_argument("--ref", action="store_true",
                   help="generate the inclusive reference script instead of "
                        "slices: single run, no MET cut (misset 0.0, "
                        "missetmax -1.0), named 'ref'")
    p.add_argument("--pythia8-card", default=None,
                   help="path of the pythia8 card, uncommented in every "
                        "launch block (e.g. ../Cards/pythia8_card_TuneA14.dat "
                        "for a script in runs/)")
    p.add_argument("--delphes-card", default=None,
                   help="path of the delphes card, uncommented in every "
                        "launch block (e.g. "
                        "../Cards/delphes_card_ATLAS_FatJetPTMin50.tcl)")
    p.add_argument("--geometric", type=float, default=None, metavar="RATIO",
                   help="geometric slicing: slice edges at min*RATIO^k up to "
                        "--max instead of fixed --width steps (adaptive "
                        "widths for a steeply falling spectrum)")
    p.add_argument("--add-jets", action="store_true",
                   help="add the lep nu j and lep nu j j multiplicity "
                        "processes (no MLM matching)")
    p.add_argument("--bias-module", default=None,
                   help="absolute path of a bias module directory, set in "
                        "every launch block (for biasing-vs-slicing tests)")
    p.add_argument("--bias-parameters", default=None,
                   help="bias_parameters dictionary string, e.g. "
                        "\"{'a': 1.0}\"")
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

    generate_lines = "generate p p > lep nu QCD=0 QED=2"
    if args.add_jets:
        generate_lines += ("\nadd process p p > lep nu j QCD=1 QED=2 @1"
                           "\nadd process p p > lep nu j j QCD=2 QED=2 @2")
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
{generate_lines}
output {args.process} --f
"""

    if args.ref:
        blocks = [launch_block(args.process, "ref", args.shower, args.detector,
                               args.analysis, True, min_param, max_param,
                               0.0, None, args.nevents, unbounded=True,
                               pythia8_card=args.pythia8_card,
                               delphes_card=args.delphes_card,
                               bias_module=args.bias_module,
                               bias_parameters=args.bias_parameters)]
        script = header + "\n".join(blocks)
        output = args.output or f"{args.process}.mg5"
        with open(output, "w") as f:
            f.write(script)
        print(f"wrote {output}: inclusive reference (no MET cut), "
              f"{args.nevents} events, shower={args.shower}, "
              f"detector={args.detector}")
        return

    lo_min, hi_max, width = args.min, args.max, args.width
    if args.geometric is not None:
        if args.geometric <= 1.0:
            p.error("--geometric RATIO must be > 1")
        if hi_max <= lo_min:
            p.error("--max must be greater than --min")
        edges = [lo_min]
        while edges[-1] < hi_max:
            edges.append(min(edges[-1] * args.geometric, hi_max))
        edges = sorted(set(round(e) for e in edges))
        width = None
    else:
        if width <= 0:
            p.error("--width must be positive")
        if hi_max <= lo_min:
            p.error("--max must be greater than --min")
        n_slices = (hi_max - lo_min) / width
        n_full = round(n_slices)
        if abs(n_slices - n_full) > 1e-9:
            print(f"WARNING: range [{lo_min:g}, {hi_max:g}] is not a multiple "
                  f"of the width {width:g}; the last slice is narrower.")
        edges = [lo_min + i * width for i in range(n_full + 1)]
        if edges[-1] < hi_max - 1e-9:
            edges.append(hi_max)

    blocks = []
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
                                   unbounded=last,
                                   pythia8_card=args.pythia8_card,
                                   delphes_card=args.delphes_card,
                                   bias_module=args.bias_module,
                                   bias_parameters=args.bias_parameters))

    script = header + "\n".join(blocks)
    output = args.output or f"{args.process}.mg5"
    with open(output, "w") as f:
        f.write(script)

    def edge_str(hi):
        return "inf" if hi == float("inf") else f"{hi:g}"

    w_str = (f"width {width:g}" if width is not None
             else f"geometric ratio {args.geometric:g}")
    print(f"wrote {output}: {len(blocks)} slices of {args.variable} "
          f"({edges[0]:g} to {edge_str(edges[-1])}, {w_str}, "
          f"{args.nevents} events each; last slice unbounded)")
    for i, (lo, hi) in enumerate(pairs):
        print(f"  [{i}] {args.name_format.format(var=args.variable, lo=lo, hi=float('inf') if i == len(pairs) - 1 else hi)}"
              f"  ->  [{lo:g}, {edge_str(hi)}]")


if __name__ == "__main__":
    main()
