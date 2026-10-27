# phase_space_slicing

Study of **phase-space slicing** with MadGraph5_aMC@NLO (MG5aMC) for the
`pp → W → ℓν` sample at 13.6 TeV: the same process is generated in
several successive **slices of the missing-transverse-energy spectrum**
(one MG5aMC run per slice, selected through the run_card MET cuts), and
the slices are combined into a single weighted sample. The goal is to
quantify:

1. **the impact of the slice width** — one wide run vs. many narrow
   slices: tail statistics per unit CPU, unweighting efficiency per
   slice, accumulation of per-slice normalization uncertainties;
2. **the border issues** — what happens to events (and to the combined
   distributions) at the boundaries between two slices, both at parton
   level (LHE) and at reconstruction level (Delphes), where detector
   resolution migrates events across the parton-level borders.

This is the classic statistics-enhancement technique that predates the
event-generation biasing studied in the companion repository
[phase_space_biasing](https://github.com/ammelsayed/phase_space_biasing):
instead of biasing the sampling inside a single run, the phase space is
split into chunks that are generated independently and recombined with
their own cross sections.

## The generation script

`gen_slicing_script.py` generates the MG5aMC script
`w012j_MLM.mg5` from a small set of knobs:

```bash
python3 gen_slicing_script.py                       # defaults: misset 0-500 GeV, 50 GeV slices
python3 gen_slicing_script.py --width 25 --max 1000 # finer/narrower slicing
python3 gen_slicing_script.py --help                # all options
```

controllable: the slice variable (default `misset`), the MET range
(`--min`/`--max`), the slice width (`--width`), the run naming style
(`--name-format`, default `{var}_{lo}_{hi}` e.g. `misset_50_100`),
the number of events per slice, and the shower/detector/analysis
switches (defaults `OFF`/`OFF`/`MadAnalysis5`, ready for the LHE-level
stage; set `--detector Delphes` for the reconstruction-level stage).

The base script is the `w012j_MLM.mg5` of the
[phase_space_biasing](https://github.com/ammelsayed/phase_space_biasing)
repository with three modifications:

* `set ptj 30.0` and `set etaj 2.5` are removed (no jets in the final
  state);
* the multiplicity label `@0` is dropped from the generate line:
  `generate p p > lep nu QCD=0 QED=2`;
* instead of a single `launch`, the script contains one **named launch
  per slice**:

  ```bash
  launch w012j_MLM --name=misset_0_50
  ... (switches, run_card settings) ...
  set misset 0.0        # slice minimum
  set missetmax 50.0    # slice maximum
  done
  launch w012j_MLM --name=misset_50_100
  ...
  ```

  so the events of slice *i* are written to `Events/misset_{lo}_{hi}/`
  and each slice gets its own reported cross section σᵢ.

The `ptl1min 150.0` cut is untouched: the slicing variable is the MET,
defined in MG5aMC as the transverse momentum of the vector sum of the
invisible particles (the neutrino of the W decay).

## Run-card facts (verified in MG5aMC v3.7.2)

These were checked directly in the v3.7.2 templates and in the cuts
applied inside a generated process (`Template/LO/Cards/run_card.dat`,
`SubProcesses/*/cuts.f`) — they matter for the design of the study:

* The MET min cut is called `misset` — there is **no**
  `missetmin` parameter. The **maximum** cut is `missetmax`, and a
  negative value (`-1.0`, the default) disables it.
* **The last slice is always unbounded** (`missetmax = -1.0`, named e.g.
  `misset_450_inf`): bounding it would silently discard all events above
  the last edge.
* At LO the lepton and the neutrino are back-to-back, so
  `MET = ptl1 ≥ ptl1min = 150` GeV for every event of this sample:
  slices below 150 GeV contain no phase space (MG5aMC reports "Zero
  result" for them), which is why the default range starts at 150 GeV.
* The cut is inclusive on **both** edges: an event passes if
  `MET ≥ misset` **and** (`missetmax < 0` **or** `MET ≤ missetmax`).
  Adjacent slices therefore *overlap* at the exact boundary value
  (no gap): an event with `MET = 50.0` GeV passes both the
  `[0, 50]` and the `[50, 100]` slice. For a continuous observable this
  is a measure-zero coincidence — the border issues that actually
  matter are the per-slice normalization discontinuities (below) and,
  at reconstruction level, the migration of events across borders.

## What we will study

**Generation** (one pass serves both stages): an inclusive reference
(`ref`, 100k events, no MET cut) plus three slice matrices
(`w50`, `w100`, `w200`: 10k events per slice, slices of the MET from
150 GeV up to an unbounded last slice), all with Pythia8 + Delphes
enabled and the cards in `Cards/`, so each run yields both the LHE and
the Delphes ROOT file.

**Stage 1 — parton level (LHE files, analysed with pylhe).** The slices
of each width are combined into one weighted sample and compared to the
reference: cross-section additivity, closure of the combined spectra,
per-bin relative statistical uncertainty as a function of the slice
width, weight discontinuities and double-counting at the borders.

**Stage 2 — reconstruction level (Delphes, analysed with pyROOT +
libDelphes).** The same comparison on `MissingET.MET` and on the
leading-lepton pT from the `Electron`/`Muon` branches, plus the
migration of events across the parton-level borders caused by the
detector resolution.

**Results**: see [study/findings.md](study/findings.md) (full report)
and `study/plots/`. Headline: slicing preserves the total cross section
to <0.2% and buys a *flat* per-bin precision profile (~2.5% up to
~950 GeV with 50 GeV slices vs ~100% in the inclusive tail at 1 TeV);
the borders are benign at parton level (no gaps, no double-counting)
and are absorbed by the detector smearing at reco level; the price is a
large event-weight spread (N_eff drops below the inclusive reference)
and per-slice normalization uncertainties accumulating in the sum.

## References

* J. Alwall *et al.*, "The automated computation of tree-level and
  next-to-leading order differential cross sections, and their matching
  to parton shower simulations", JHEP **07** (2014) 079,
  [arXiv:1405.0301](https://arxiv.org/abs/1405.0301) — MG5aMC framework:
  run_card cuts, unweighting and per-run cross-section normalization.
* [LOEventGenerationBias wiki page](https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/LOEventGenerationBias)
  — MG5aMC's built-in event-generation biasing (available since
  v2.5.0), the modern in-run alternative to slicing for tail
  enhancement; studied in the companion repository.
* The MG5aMC v3.7.2 `Template/LO/Cards/run_card.dat` (definition of the
  `misset` / `missetmax` cuts) and the generated
  `SubProcesses/*/cuts.f` (their inclusive implementation), verified in
  this study and documented above.
* DELPHES 3 Collaboration (J. de Favereau *et al.*), "DELPHES 3, A
  modular framework for fast simulation of a generic collider
  experiment", JHEP **02** (2014) 057,
  [arXiv:1307.6346](https://arxiv.org/abs/1307.6346) — for the
  reconstruction-level stage.
* E. Conte, B. Fuks, G. Serret, "MadAnalysis 5, a user framework for
  qualitative analysis of collider events", Comput. Phys. Commun.
  **184** (2013) 222, [arXiv:1206.1599](https://arxiv.org/abs/1206.1599).
