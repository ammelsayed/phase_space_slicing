# phase_space_biasing

Study of the effects of **phase-space biasing** on V+jets production with
MadGraph5_aMC@NLO (MG5aMC): `pp → Z/W → leptons + 0,1,2 jets`, MLM merged,
13.6 TeV (`ebeam 6800`).

The idea: MG5aMC can bias the LO event generation so that more unweighted
events are produced in a region of phase space that is statistically
limited (e.g. the high-pT tail), while the physics content of the sample
is preserved by the corresponding event weights. This repository collects
the generation scripts and a set of biasing modules to quantify how much
statistics can be gained in the tails of the `ptl1`, `MET`, `HT` and
`ptj` distributions of V+jets samples.

## Repository layout

```
phase_space_biasing/
├── z012j_MLM.mg5                  # generation script, Z(→l+l-) + 0,1,2j (MLM)
├── w012j_MLM.mg5                  # generation script, W(→lν)   + 0,1,2j (MLM)
├── ptj_bias.f                     # verbatim copy of the shipped MG5aMC ptj_bias module
└── mg5_bias_modules/              # biasing modules (each one is a self-contained module dir)
    ├── ptj_bias/                  # the shipped module (reference/template)
    ├── ptl1_bias/                 # bias on the pt of the leading charged lepton
    ├── met_bias/                  # bias on the missing transverse energy
    ├── ht_bias/                   # bias on HT (scalar sum of visible pT)
    ├── ptl1_met_bias/             # combined bias: ptl1 × MET (no HT term)
    ├── ptl1_met_ht_bias/          # combined bias: ptl1 × MET × HT
    └── xsec_flat_bias/            # (experimental) bias on the cross section itself
```

Every module directory contains `<name>.f` (the Fortran implementation of
`bias_wgt`) and a `makefile`, following the layout of the shipped
`Template/LO/Source/BIAS/ptj_bias` module of MG5aMC.

## Generation scripts

`z012j_MLM.mg5` and `w012j_MLM.mg5` are self-contained MG5aMC scripts
(process definition + launch). Run them with:

```bash
mg5_aMC z012j_MLM.mg5
mg5_aMC w012j_MLM.mg5
```

Settings (as requested for this study):

| Setting | Value |
|---|---|
| Events | `set nevents 100_000` |
| Random seed | `set iseed 0` (chosen from the clock; no multi-seed scan) |
| Shower | `shower = OFF` (no Pythia8) |
| Detector | `detector = OFF` (no Delphes) |
| Analysis | `analysis = MadAnalysis5` |
| MadSpin / reweight | OFF |
| Merging | MLM settings (`ickkw`, `xqcut`, `maxjetflavor`) commented out — no matching |
| Cuts | `ptj 30`, `etaj 2.5`, `ptl1min 150`, `ptl2min 50` / `misset 50` (W), `etal 2.5` |
| PDF | `lhapdf`, `lhaid 331900` |

The process-definition block (`define`, `generate`, `output`, diagram
plotting) is uncommented, so each script regenerates the process from
scratch. Since Pythia8 and Delphes are OFF, the two
`./Cards_edited/pythia8_card_*.dat` / `delphes_card_*.tcl` lines are
commented out (they would otherwise be dangling answers to questions that
are no longer asked).

Prerequisites:

* A PDF set with id `331900` installed in LHAPDF.
* MadAnalysis5 for the analysis stage: install it (e.g. `install
  madanalysis5` from the MG5aMC prompt) and make sure `madanalysis5_path`
  is set in `mg5_configuration.txt`, otherwise `analysis = MadAnalysis5`
  cannot run.

## How MG5aMC phase-space biasing works

Documented in the [LOEventGenerationBias wiki
page](https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/LOEventGenerationBias),
available since MG5aMC v2.5.0. Summary:

* The run_card entry `bias_module` selects the module (`None` by default,
  the name of a shipped module, or **an absolute path to a custom module
  directory**).
* Its double-precision parameters are set through the run_card entry
  `bias_parameters` (a python-style dictionary), or from the launch
  prompt with `set bias_parameter <name> <value>`.
* During unweighting, events are kept with a probability proportional to
  `|original_weight × bias_weight|`, so the bias reshapes where in phase
  space the unweighted events land.
* Each module declares, inside the Fortran code, the common block

  ```fortran
  double precision stored_bias_weight
  logical impact_xsec, requires_full_event_info
  common/bias/stored_bias_weight,impact_xsec,requires_full_event_info
  ```

  * `impact_xsec = .False.` : the bias is **only a sampling device**. The
    total cross section and the physical distributions are unchanged;
    the unweighted events come out *weighted* and each event's bias
    weight is recorded in the LHE file as `<rwgt><wgt id='bias'>`,
    so any analysis can undo the bias a posteriori. **All modules in
    this repository use this mode.**
  * `impact_xsec = .True.` : the bias is a genuine modification of the
    matrix-element weight (it changes the cross section and the
    distributions themselves).
  * `requires_full_event_info = .False.` : the module only needs the
    event kinematics `p(0:3,nexternal)` (faster). Set it to `.True.`
    (see the `PY8_CKKW` example on the wiki) if you need the full LHE
    event record, colors, helicities, etc.
* A custom module is a directory named after the module containing
  `<name>.f` (with the `subroutine bias_wgt(p, original_weight,
  bias_weight)` entry point) and a `makefile`. Parameters are declared
  with a comment inside the subroutine,

  ```fortran
  C  parameters = {'ptl1_bias_target_ptl1': 300.0,
  C               'ptl1_bias_enhancement_power': 4.0}
  ```

  MG5aMC parses it, creates the corresponding run_card entries and
  writes their values into `Source/BIAS/bias.inc`, which the module
  includes after the variable declarations.
* Cluster-mode caveat (from the wiki): a bias module must not read
  external files at runtime — compile everything in.

## The biasing modules

All modules have the functional form
`bias_wgt = (observable / target)^power` and leave events without the
relevant object unbiased (weight 1). Defaults below are what the run_card
gets if you don't override `bias_parameters`.

| Module | Observable | Default parameters |
|---|---|---|
| `ptj_bias` (shipped, copied here) | pt of the hardest jet | `target 1000`, `power 4` |
| `ptl1_bias` | pt of the hardest charged lepton (e±, μ±, τ±) | `ptl1_bias_target_ptl1 300`, `ptl1_bias_enhancement_power 4` |
| `met_bias` | MET = \|Σ pT of invisible particles (neutrinos)\| | `met_bias_target_met 300`, `met_bias_enhancement_power 4` |
| `ht_bias` | HT = scalar Σ pT of all visible particles | `ht_bias_target_ht 1000`, `ht_bias_enhancement_power 4` |
| `ptl1_met_bias` | product of the ptl1 and MET biases above (no HT term) | all four ptl1/met parameters above |
| `ptl1_met_ht_bias` | product of the ptl1, MET and HT biases above | all six parameters above |
| `xsec_flat_bias` | (experimental) see below | `xsec_flat_bias_exponent 1.0`, `xsec_flat_bias_reference_weight 1.0` |

### xsec_flat_bias — biasing the cross section directly

Instead of biasing one kinematic distribution, this module tries to bias
the *cross section itself*: it uses the `original_weight` argument of
`bias_wgt`, which is proportional to the differential cross section at
that phase-space point, and inverts it:

```fortran
bias_wgt = xsec_flat_bias_reference_weight
           / |original_weight|^xsec_flat_bias_exponent
```

With exponent = 1 the combined weight `original_weight × bias_weight` is
flat, i.e. **every phase-space point is generated with the same rate** —
regions with a small cross section (which would normally be sampled
rarely) get as many events as the regions with a large one. `0 < exponent
< 1` gives a partial flattening, `0` recovers the unbiased generation.

It runs with `impact_xsec = .False.`, so the reported cross section and
the physical distributions are untouched; only the *population* of the
unweighted sample changes, and each event's bias weight is stored in the
LHE `<rwgt>` block for a-posteriori correction. Caveats (documented in
the module header): the flattening is in the generator's integration
measure (which already includes its own importance sampling), not
uniformly in a physical observable, and with exponent = 1 the bias
weights can be huge — if the unweighting efficiency degrades, lower the
exponent.

## Using the modules

Point `bias_module` at the absolute path of the module directory:

```bash
mg5_aMC> set bias_module /path/to/this/repo/mg5_bias_modules/ptl1_bias
mg5_aMC> set bias_parameters {'ptl1_bias_target_ptl1': 500.0, 'ptl1_bias_enhancement_power': 4.0}
```

or, per-parameter from the launch prompt:

```bash
mg5_aMC> set bias_parameter ptl1_bias_target_ptl1 500.0
```

To disable: `set bias_module None`.

## How to study phase-space biasing with this repository

1. **Reference sample**: generate `z012j_MLM` / `w012j_MLM` with
   `bias_module None`.
2. **Biased samples**: for each module, rerun the *same script* with the
   same `nevents` and cuts, changing only `bias_module` (and optionally
   its parameters).
3. **Compare**, per module:
   * the relative statistical uncertainty of the biased vs. reference
     distributions in the tails (`ptl1`, `MET`, `HT`, `ptj`), where the
     bias should buy you events;
   * the total cross section — with `impact_xsec = .False.` it must agree
     with the reference within uncertainties;
   * the central values — reweighting each biased event by
     `1/bias_weight` (or simply dividing the stored event weight by the
     `<wgt id='bias'>` entry) must reproduce the reference distributions,
     which is the basic closure test of the whole procedure.
4. **`xsec_flat_bias`**: additionally check how the unweighted events are
   spread over the phase space (population of normally rare regions) and
   the unweighting efficiency printed by MG5aMC, as a function of
   `xsec_flat_bias_exponent`.

## Results of the biasing study

The study was run on V→leptons samples only (the `add process … j` lines
were dropped, so there are no final-state jets): `pp → Z → τ⁺τ⁻` and
`pp → W → τν` at 13.6 TeV, 100 000 unweighted events each, all with the
same RNG seed (`iseed 0` → MG5aMC v3.7.2 used seed offset 21 for every
run) and `impact_xsec = .False.`. Since there are no jets, `ptj_bias`,
`ht_bias` and `ptl1_met_ht_bias` were **not** run, and the Z sample does
not use any MET-containing module. The tested matrix:

| Sample | Runs |
|---|---|
| Z | reference, `ptl1_bias`, `xsec_flat_bias` |
| W | reference, `ptl1_bias`, `met_bias`, `ptl1_met_bias` (default 4+4), `ptl1_met_bias` with `bias_parameters` powers 2+2, `xsec_flat_bias` |

Reproduce with `study/run_study.sh` (uses `/data/ammelsayed/runmg5.sh`);
`study/analyze_lhe.py` builds the plots in `study/plots/` and the numbers
below. MadAnalysis5 reports generated during the runs are collected in
`study/ma5/`.

Reported cross sections (with `impact_xsec = .False.` the bias must not
change the physics):

| Run | σ [pb] |
|---|---|
| Z reference | 0.812 ± 0.0007 |
| Z `ptl1_bias` | 0.8106 ± 0 |
| Z `xsec_flat_bias` | 0.8152 ± 0 |
| W reference | 2.301 ± 0.0027 |
| W `ptl1_bias` | 2.312 ± 0 |
| W `met_bias` | 2.312 ± 0 |
| W `ptl1_met_bias` (4+4) | **1.088 ± 0 (broken!)** |
| W `ptl1_met_bias` (2+2) | 2.312 ± 0 |
| W `xsec_flat_bias` | 2.295 ± 0 |

Event counts N and per-bin relative statistical uncertainty
δ = √(Σw²)/Σw above pT thresholds (format `N / δ / gain = δ_ref/δ`):

**Z sample, ptl1** (identical gains hold for ptl2 and HT, see below):

| threshold | reference | `ptl1_bias` | `xsec_flat_bias` |
|---|---|---|---|
| 300 GeV | 9 046 / 0.011 / – | 72 754 / 0.006 / **1.7×** | 8 986 / 0.013 / 0.8× |
| 500 GeV | 1 068 / 0.031 / – | 45 005 / 0.007 / **4.6×** | 1 100 / 0.036 / 0.9× |
| 800 GeV | 109 / 0.096 / – | 21 182 / 0.009 / **11.2×** | 107 / 0.116 / 0.8× |

**W sample, ptl1 (and MET, HT — degenerate, see below):**

| threshold | reference | `ptl1_bias` = `met_bias` = `ptl1_met_bias` (2+2) | `ptl1_met_bias` (4+4) | `xsec_flat_bias` |
|---|---|---|---|---|
| 300 GeV | 8 260 / 0.011 / – | 71 347 / 0.006 / **1.8×** | 99 903 / 0.050 / 0.2× | 9 089 / 0.013 / 0.9× |
| 500 GeV | 1 025 / 0.031 / – | 44 612 / 0.007 / **4.7×** | 98 959 / 0.019 / 1.7× | 1 159 / 0.035 / 0.9× |
| 800 GeV | 116 / 0.093 / – | 21 564 / 0.008 / **11.0×** | 93 280 / 0.010 / 9.7× | 132 / 0.103 / 0.9× |

Plots (top: physical weighted distributions — closure; middle: unweighted
population; bottom: δ_biased/δ_ref): [z_ptl1](study/plots/z_ptl1.png),
[z_ptl2](study/plots/z_ptl2.png), [z_ht](study/plots/z_ht.png),
[w_ptl1](study/plots/w_ptl1.png), [w_met](study/plots/w_met.png),
[w_ht](study/plots/w_ht.png).

Findings:

1. **Closure works.** With `impact_xsec = .False.` the weighted
   distributions of every healthy biased sample reproduce the reference
   over five decades (top panels): MG5aMC v3.7.2 stores the physical
   per-event weight in the LHE weight column (no `<rwgt>` block is
   written), so biased samples are directly usable in analyses.
2. **The tails gain a lot.** Above 300/500/800 GeV the same 100k-event
   budget yields respectively ~2×/~5×/~11× smaller statistical
   uncertainties, while the lowest bins lose at most a factor ~2–3.
3. **At LO without jets all key observables are degenerate.** The W is
   produced with zero transverse momentum, so `ptl1 = ptl2 = MET = HT`
   (Z: `ptl1 = ptl2 = HT/2`) event by event. Consequently `met_bias` is
   *numerically identical* to `ptl1_bias` for the W (the two runs produced
   bit-identical event files), and biasing the single variable ptl1
   improves **all** of ptl1, ptl2, MET and HT simultaneously.
4. **Over-biasing breaks the generation.** The combined `ptl1_met_bias`
   with default parameters (powers 4×4 = effective power 8 for the W)
   converged to a cross section 53% too low (1.088 instead of 2.301 pb,
   reported as "± 0" and "stable") and its weighted distribution does not
   close: the integration/unweighting cannot handle such an aggressive
   bias. With powers 2+2 (effective power 4) the same module is perfectly
   healthy. Practical rule: keep the *effective total* power of the bias
   ≈ 4 or below, and always check the reported σ against an unbiased run.
5. **Biasing the cross section directly is a no-op for the event
   population.** `xsec_flat_bias` (bias ∝ 1/|original_weight|, exponent 1)
   is the *healthiest* run of the study: the flattened total weight makes
   MG5aMC's integration trivially stable — σ within 0.3–0.6% of the
   reference on both samples, ~100% unweighting efficiency (fastest runs,
   ~1–2 min) — but its unweighted population is statistically
   indistinguishable from the reference (gains 0.8–0.9× in every tail,
   e.g. 107 events above ptl1 = 800 GeV vs 109 for the reference). The
   1/|w| factor cancels the importance sampling that MG5aMC's VEGAS
   integration has already built in, so no region of phase space ends up
   oversampled. To actually move statistics between regions one must bias
   an observable that compensates the *physical* falloff (as `ptl1_bias`
   does), not the event weight itself.

**Verdict** — for making sure there is enough statistics in *almost all*
distributions of these V→leptons samples, **`ptl1_bias` with its default
parameters (target 300 GeV, power 4) is the best module** for both the Z
and the W sample: it is the only knob needed at LO (all key observables
are degenerate with ptl1), it improves every tail by 2–11× at fixed event
budget, and it preserves the total cross section and the physical
distributions. Combined modules (`ptl1_met_bias`) only become useful when
the observables decouple — i.e. once jets/recoil are present — and must
then be run with moderate exponents (effective total power ≲ 4).
`xsec_flat_bias` is technically the most robust module (perfect σ
stability, 100% unweighting efficiency) but does not redistribute
statistics at all, so it does not help the tails.

## References

* [LOEventGenerationBias](https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/LOEventGenerationBias)
  — official MG5aMC wiki documentation of the biasing feature, including
  the custom-module interface and the `PY8_CKKW` example used as the
  advanced reference.
* The shipped module `Template/LO/Source/BIAS/ptj_bias` of MG5aMC —
  copied verbatim in this repository (`ptj_bias.f` and
  `mg5_bias_modules/ptj_bias/`); every custom module here follows its
  structure.
* MG5aMC is available since v2.5.0 with this feature; the framework is
  described in J. Alwall *et al.*, "The automated computation of
  tree-level and next-to-leading order differential cross sections, and
  their matching to parton shower simulations", JHEP **07** (2014) 079,
  [arXiv:1405.0301](https://arxiv.org/abs/1405.0301).
* MadAnalysis5: E. Conte, B. Fuks, G. Serret, "MadAnalysis 5, a user
  framework for qualitative analysis of collider events", Comput. Phys.
  Commun. **184** (2013) 222, [arXiv:1206.1599](https://arxiv.org/abs/1206.1599).
