# Phase-space slicing: findings

**Sample**: `pp → W → ℓν` at 13.6 TeV (MG5aMC v3.7.2), `ptl1min = 150` GeV,
no other MET/jet cuts, Pythia8 shower + Delphes reconstruction (the cards
in `Cards/`), 40 cores, same RNG seed for every run.

**Variants**:

| Variant | Runs | Events | Total events |
|---|---|---|---|
| `ref` (inclusive, no MET cut) | 1 | 100 000 | 100 000 |
| `w200` (200 GeV slices, 150→950 + [950, ∞)) | 5 | 10 000 | 50 000 |
| `w100` (100 GeV slices) | 9 | 10 000 | 90 000 |
| `w50` (50 GeV slices) | 17 | 10 000 | 170 000 |

Every event keeps the LHE weight of its run, rescaled so that
`Σw = σ_i` per slice; the combined sample of a variant is simply the
concatenation of its slices. Parton level = `unweighted_events.lhe.gz`
(read with `pylhe`); reconstruction level = `tag_1_delphes_events.root`
(read with pyROOT + `libDelphes`, observables: `MissingET.MET[0]` and
the highest `Electron.PT` / `Muon.PT`).

## 1. Cross-section additivity: exact

| Variant | Σσᵢ [pb] | σ(ref) [pb] | pull |
|---|---|---|---|
| w200 | 2.3010 ± 0.0091 | 2.3010 ± 0.0027 | +0.00 |
| w100 | 2.2976 ± 0.0090 | 2.3010 ± 0.0027 | −0.36 |
| w50  | 2.3003 ± 0.0064 | 2.3010 ± 0.0027 | −0.10 |

Slicing conserves the total cross section: the per-slice integrations add
up to the inclusive value within ~0.2%. The uncertainty of the summed
normalization (±0.006–0.009 pb) is 2–3× larger than the reference's
(±0.0027 pb) because each slice's integration uses only 10k events.

Per-slice closure: each σᵢ was compared with the reference's content in
the same MET range; the worst |pull| over all slices is
0.98 (w200) / 1.43 (w100) / 1.83 (w50) — consistent with statistics
(17 slices × ~1σ trials).

## 2. Parton level (LHE): the precision profile is the whole story

Per-bin relative statistical uncertainty δ = √(Σw²)/Σw of the combined
samples vs the reference (10 GeV bins, `study/plots/relerr_profiles.png`,
`closure_w50.png`):

* **Reference (black)**: δ grows from ~0.7% at 160 GeV to ~100% at
  1 TeV — the tail is starved.
* **w50 (red)**: flat at ~2.5% from 200 GeV to 950 GeV, rising to ~6%
  in the unbounded last slice.
* **Crossover at ~215 GeV**: below it the inclusive run wins (its 100k
  events pile up in the bulk); above it slicing wins, by a factor
  ~40 in δ at 1 TeV.
* The δ profile of the sliced samples shows a **sawtooth with the slice
  period**: within each slice the population falls (the physical
  spectrum drops ~3× per 50 GeV), resetting at every border. Narrower
  slices → smaller teeth → more uniform precision.
* Closure: the combined weighted spectra reproduce the reference within
  its statistical errors over four decades (ratio panels of
  `closure_w*.png`), with no visible steps at the borders.

## 3. The border issues, quantified (parton level)

* **Exact-border double counting**: zero. The MG5aMC
  `misset`/`missetmax` cuts are inclusive on both edges, but in ~1.4M
  combined events not a single one sat exactly on a border value (the
  observable is continuous).
* **Weight discontinuity**: the per-event weight jumps at every border
  by σᵢ₊₁/σᵢ, i.e. **×24.7 for w200, ×7.0 for w100, ×3.1 for w50**
  (`weight_steps.png`). This is *the* border artifact an analysis must
  handle: events at the same MET carry different weights depending on
  which slice generated them. Any analysis that treats events as
  unweighted (or that bins cross-border) mixes these scales.
* **Effective sample size**: N_eff = (Σw)²/Σw² of the combined samples
  is 1.1e4 (w200) / 1.4e4 (w100) / 2.2e4 (w50) vs 1e5 for the
  reference. Despite having 1.7× the events, the w50 combination is
  **~4.5× worse** than the reference for any measurement that is
  sensitive to the weight spread (e.g. the total rate). Slicing buys
  *uniform* precision, not better total precision.

## 4. Reconstruction level (Delphes): migration absorbs the borders

(`study/plots/closure_*_{met,ptl1}.png`, `relerr_profiles_delphes_*.png`,
`parton_vs_reco_w50.png`, `migration_w50.png`)

* **Migration is large and asymmetric**: for w50, 24 171 events (14%)
  reconstruct above their slice's upper border and **75 630 events (44%)
  reconstruct below their slice's lower border**. The asymmetry is
  kinematic: within each slice the spectrum falls steeply, so most
  events sit near the slice's low edge, and the MET resolution smears
  them across it. At the borders of the *first* slice the effect is
  mildest; it grows with MET because the spectrum gets steeper.
* **The combined reco-MET spectrum is smooth**: the stacked per-slice
  contributions (`migration_w50.png`) show each slice leaking into its
  neighbours such that the sum follows the reference everywhere — the
  per-event weights absorb the migration automatically. No visible
  step at any border in the ratio panels.
* **The δ profile at reco level keeps the sawtooth** (now with period =
  slice width, `relerr_profiles_delphes_met.png`): bins just above a
  border are fed by two slices (the migrated low-edge events of the
  upper slice + the up-migrating tail of the lower slice) and have the
  best precision; bins in the upper part of each slice have the worst.
  Narrow slices (w50) minimize the depth of the teeth.
* **ptl1 at reco level** (`relerr_profiles_delphes_ptl1.png`) behaves
  like MET, as it must at LO (MET ≡ ptl1 event by event before
  smearing): the same flat-vs-rising structure with a crossover at
  ~300 GeV. The lepton pT resolution is much smaller than the MET
  one, so *less* migration washes the borders out and the ptl1
  sawtooth is correspondingly sharper than the MET one.

## 5. Conclusions

1. **Slicing works**: with 10k events per slice, the combined sample
   reproduces the inclusive spectrum at parton *and* reconstruction
   level, and the total cross section is preserved to better than 0.2%.
2. **The slice width trades bulk precision for tail precision**:
   at fixed 10k/slice, narrower slices give a more *uniform* δ profile
   (w50: flat ~2.5% at parton level) and smaller border weight steps
   (×3.1 vs ×24.7 for w200), at the price of more runs and more total
   events.
3. **The borders are benign at parton level** (no gaps, no
   double-counting, exact additivity) and **absorbed by the detector
   smearing at reco level** — provided analyses use the event weights
   (which carry the per-slice normalization) and do not treat the
   combined sample as unweighted.
4. **The price of slicing is the weight spread**: N_eff drops below the
   inclusive reference for total-rate-like measurements, and the
   per-slice normalization uncertainties accumulate in the sum. Slicing
   is the right tool when the *distributions in the tail* are the
   target; it is the wrong tool for measuring the total cross section.

## Reproducing

```bash
python3 gen_slicing_script.py --ref --nevents 100_000 --shower Pythia8 \
    --detector Delphes --analysis OFF \
    --pythia8-card ../Cards/pythia8_card_TuneA14.dat \
    --delphes-card ../Cards/delphes_card_ATLAS_FatJetPTMin50.tcl \
    --process ref --output runs/ref.mg5
# same with --max 1000 --width {50,100,200} --nevents 10_000 --process w{50,100,200}
sh /data/ammelsayed/runmg5.sh runs/ref.mg5   # etc. (one per variant)

# parton level (pylhe)
/data/ammelsayed/venv/bin/python study/analyze_slicing.py runs study/plots
# reconstruction level (pyROOT + libDelphes)
source /home/ammelsayed/root/bin/thisroot.sh
/data/ammelsayed/venv/bin/python study/analyze_delphes.py runs study/plots
```
