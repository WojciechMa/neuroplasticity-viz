# Intrinsic excitability and f-I curves: project plan

Source: [assignment text](../assignemnt-text.txt) (Assignment 14*, computational).

## What we are doing

Neurons adapt not only through their synapses but also through their intrinsic excitability. Changing sodium, potassium or leak conductances shifts the firing-rate-versus-current (f-I) curve. The assignment has three parts:

1. **Data.** Get Patch-seq current-clamp recordings from the Allen Institute, keeping only cells that have both electrophysiology and a transcriptomic label.
2. **Empirical f-I curves.** Classify each cell as excitatory (E) or inhibitory (I) from its transcriptomic type and build an f-I curve for each group.
3. **Model perturbation.** Fit a simple neuron model to a representative E cell and a representative I cell, vary the leak conductance $g_L$, and compare the resulting f-I curves.

Part 3 is a what-if in a model. It does not show that the recorded cells underwent activity-dependent plasticity, and the write-up must say so.

## Data

### What the linked page gives

`data_from_website/cell_types_specimen_details.csv` is the "Download cell feature data" file from the Cell Feature Search page. It covers 2,333 cells (1,920 mouse, 413 human) from the older Cell Types survey:

- one row of summary features per cell (f-I slope, rheobase, input resistance, time constant, resting potential), plus Cre line and dendrite type;
- no transcriptomic labels and no voltage traces;
- none of its cells are in the Patch-seq tables.

So it cannot answer the assignment. At most it serves as a sanity check for our feature values.

### Patch-seq data we use: mouse visual cortex

Every cell below has a recording and a transcriptomic type, which is the assignment's filter.

| Group | Study | Transcriptomic type | Recordings |
|---|---|---|---|
| I | Gouwens et al. 2020 | Allen `20200711_patchseq_metadata_mouse.csv`, column `T-type Label`: 4,269 GABAergic cells; 165 without a t-type are dropped | [DANDI 000020](https://dandiarchive.org/dandiset/000020), files `sub-*_ses-<ephys_session_id>_icephys.nwb` |
| E | Sorensen et al. | `derived_data/inferred_met_types.csv` in [AllenInstitute/exc_vis_manuscript](https://github.com/AllenInstitute/exc_vis_manuscript): 1,528 glutamatergic cells | `exc_patchseq_file_manifest.csv` in the same repository: a DANDI download link per cell |

- 75 L2/3 IT cells appear in both tables with the same recording. Only the excitatory study gives them a t-type, so they count as E.
- 50 cells per group are drawn at random (`random_state=0`) from all labelled cells. Their 100 NWB files come to about 3 GB.
- Human Patch-seq exists, but the Allen human metadata covers only excitatory L2/3 cells (318). Human interneurons come from a separate study, so human data is not used.

The download URLs are in [config.py](../src/fi_analysis/config.py); [01_fetch_data.py](../scripts/01_fetch_data.py) downloads the metadata into `data/metadata/`.

### Inside an NWB file

- Voltage, sampled at 50 kHz, is in `acquisition/` (`CurrentClampSeries`). The command current is in `stimulus/presentation/`; the two are matched by sweep number.
- `stimulus_description` names the protocol. A survey of all 100 cells found two protocols:

| Protocol | Cells (E / I) | 1-s current steps | Current range |
|---|---|---|---|
| Older, full | 8 / 21 (+1 I with both) | `C1LSCOARSE*`, `C1LSFINEST*`, `subthreshold_*`, `Rheo_Supra*` | 20-pA steps from −110 pA to well above rheobase |
| Newer, shortened | 42 / 28 | `X1PS_SubThresh`, `X3LP_Rheo`, `X4PS_SupraThresh`, `Y3LP_Rheo`, `Y4PS_SupraThresh` | above threshold only up to about rheobase + 60 pA |

- Every cell reaches at least rheobase + 60 pA, so this is the shared range used for the statistics.
- Left out: the 0.5-s steps (`X2LP`, `Y2LP`, `Rheobase_finder`), ramps and short pulses. Some `X1PS` sweeps hold only the 10-ms test pulse, so the loader keeps only steps of 1 s ± 50 ms.
- Aborted sweeps are padded with exact zeros (not NaN). Before this was caught, they produced fake spikes at 0 mV and negative input resistances. The loader drops sweeps with more than 100 exact zeros before the step ends.

## Method

### E/I classification

The first word of a t-type label names its subclass:

- Inhibitory: `Lamp5`, `Sncg`, `Vip`, `Sst`, `Pvalb`, `Serpinf1`, `Meis2`
- Excitatory: `L2/3 IT`, `L4 IT`, `L5 ET`, `L5 NP`, `L6 CT`, `L6b`, and so on
- Anything else (unlabelled) is dropped.

### f-I curve of one cell

- Stimulus: the 1-s current steps listed above.
- Spike: an upward crossing of $dV/dt \ge 20$ mV/ms whose peak reaches above −20 mV. This puts the threshold a few mV higher than the Allen definition (5 % of the maximum $dV/dt$), which matters for $V_{th}$ in the model.
- Firing rate: number of spikes during the step ÷ step duration.
- Per-cell features:
  - rheobase (smallest step that evokes a spike);
  - f-I slope (linear fit from rheobase upwards);
  - maximum rate;
  - input resistance $R_{in}$ and membrane time constant $\tau_m$, from the hyperpolarising steps;
  - resting potential;
  - spike threshold and post-spike trough, from the rheobase sweep.

### Quality control

- Sweep: the 100-ms pre-step baseline must lie within 5 mV of the cell's median baseline and have SD ≤ 1 mV.
- Cell: it needs hyperpolarising steps, at least one spike, at least two points in the shared range, and a measurable $R_{in}$ and $\tau_m$. All 100 cells pass.

### Group f-I curves

- Average within each cell first, then across cells: mean ± SEM per current bin.
- Two x-axes: absolute current, and current above rheobase. The second separates threshold shifts from gain changes.
- Compare the features between E and I with a two-sided Mann-Whitney U test, for all cells and again for the newer protocol only. Slope and rate use the shared range.

### Model: leaky integrate-and-fire (LIF)

$$C\frac{dV}{dt} = -g_L(V-E_L) + I$$

When $V$ reaches $V_{th}$, the model fires a spike, resets to $V_r$ and stays there for $t_{ref}$. For a constant current the firing rate has a closed form:

$$I_{rh} = g_L(V_{th}-E_L),\qquad f(I)=\left[t_{ref}+\frac{C}{g_L}\ln\frac{I-g_L(V_r-E_L)}{I-g_L(V_{th}-E_L)}\right]^{-1}\quad (I>I_{rh})$$

Parameters from the recorded cell:

| Model | From recording |
|---|---|
| $g_L$ | $1/R_{in}$ |
| $C$ | $\tau_m/R_{in}$ |
| $E_L$ | resting potential |
| $V_{th}$ | $E_L + I_{rh}/g_L$, with $I_{rh}$ midway between the last silent and the first spiking step |
| $V_r$, $t_{ref}$ | least-squares fit to the recorded f-I points; $V_r$ at most 20 mV below rest |

Why not the measured values:

- The measured spike threshold (about −34 mV, from $dV/dt$) lies far above the effective LIF threshold (about −55 mV). Using it would make the model rheobase 2–5 times too high.
- A free fit drives $V_r$ to about −200 mV, because only a huge reset can flatten the LIF curve.
- The LIF rate rises almost vertically at rheobase, whereas the 1-s mean rate of a real (adapting) cell rises slowly. The fitted $t_{ref}$ (about 55–60 ms) is therefore an effective stand-in for adaptation, not a real refractory period.

**Perturbation:** scale $g_L$ by 0.5, 0.75, 1, 1.5 and 2, keeping $C$ fixed, so $\tau_m = C/g_L$ changes too.

**What to compare:** how rheobase and f-I slope move in the E and the I model, and whether a leak change alone can move the E curve towards the I curve.

**Expectation to check, not assume:** in a LIF model, more leak mainly shifts the curve to the right (raises rheobase). At high current the slope approaches $1/(C(V_{th}-V_r))$ whatever $g_L$ is, so the effect is subtractive rather than divisive (Holt & Koch 1997).

**Limits of the model:** no spike-frequency adaptation, no spike shape, and no separate Na/K conductances.

### Model extension: adding sodium and potassium channels

Na and K channels can be added to the integrate-and-fire model. Once the channels produce the spike, the threshold-and-reset rule is no longer needed, and the model becomes a single-compartment Hodgkin–Huxley-type model:

$$C\frac{dV}{dt} = -g_L(V-E_L) - \bar g_{Na}\, m^3 h\,(V-E_{Na}) - \bar g_K\, n^4 (V-E_K) - \bar g_M\, p\,(V-E_K) + I$$

- **Starting point:** the minimal cortical models of Pospischil et al. 2008. Use "regular spiking" for E (it includes a slow M-type K current $\bar g_M$, which gives adaptation) and "fast spiking" for I (no M current).
- **From the recorded cell:** $g_L$, $C$ and $E_L$, as for the LIF model.
- **From the template:** $\bar g_{Na}$ and $\bar g_K$ as densities (S/cm²) multiplied by the membrane area implied by $C$ at 1 µF/cm². Kinetics, $E_{Na} = 50$ mV and $E_K = -100$ mV follow ModelDB 123623.
- **Fitted:** the kinetic offset $V_T$ (−75 to −35 mV in 1-mV steps) and $\bar g_M$ (0 to 4·10⁻⁴ S/cm²), by grid search on the recorded f-I points.
- **Integration:** exponential Euler, dt = 0.05 ms (same spike counts as 0.005 ms), with 200 ms at zero current before each 1-s step.
- **Cost:**
  - The f-I curve has no closed form, so each point is a simulation. The simulations are vectorised over currents and parameters, so a full fit takes about 15 s.
  - The model has more parameters than one cell's f-I curve can pin down, so the channel kinetics come from the literature.
- **Gain:**
  - The same model can also vary $\bar g_{Na}$ and $\bar g_K$, which the assignment's introduction mentions.
  - It shows how leak interacts with adaptation and spike generation, which a LIF model cannot.
- **Middle ground:** the adaptive exponential integrate-and-fire model (AdEx; Brette & Gerstner 2005). Its exponential term stands in for Na activation, and its adaptation variable for a slow K current.

Order: first the LIF model, which is analytic and isolates the pure leak effect; then the Hodgkin–Huxley-type model as the extended version, with $g_L$, $\bar g_{Na}$ and $\bar g_K$ each scaled by 0.5–2. The two models are then compared on how leak acts. Results are in [report.md](report.md).

## Code layout

```
assignment/
  assignemnt-text.txt      original assignment
  data_from_website/       Cell Feature Search download (summary features only, see Data)
  docs/plan.md             this file
  docs/report.md           results and caveats
  pyproject.toml           package and dependencies
  src/fi_analysis/
    config.py              paths, download URLs, analysis constants
    data_access.py         metadata download and join, NWB download
    cell_class.py          t-type label -> E / I
    ephys.py               NWB sweeps -> spikes, rates, passive properties
    fi_curves.py           per-cell f-I and features, group curves, representative cell
    model.py               LIF (analytic f-I, fit, simulation) and HH-type model (simulation, fit)
  scripts/
    01_fetch_data.py       select cells, download NWB files
    02_build_fi_curves.py  per-cell features, QC, group f-I curves, E vs I statistics
    03_leak_perturbation.py  vary g_L in LIF fitted to representative cells
    04_hh_model.py         vary g_L, g_Na, g_K in HH-type model; LIF vs HH comparison
  tests/                   unit tests for the pure functions
  data/                    downloaded and processed data (git-ignored, created by the scripts)
  results/                 tables and figures (git-ignored, created by the scripts)
```

Setup (from `assignment/`):

```
python -m venv .venv
.venv\Scripts\activate
pip install -e .[dev]
pytest
```

The scripts are run in order (`python scripts/01_fetch_data.py`, and so on). All four have been run on the full sample.

## Decisions

| Topic | Decision | Reason |
|---|---|---|
| Data source | Allen Patch-seq: DANDI recordings plus Allen metadata | The Cell Feature Search file has no transcriptomic labels |
| Cells | Mouse visual cortex; 50 random cells per group, each with a recording and a t-type | The assignment's filter. Only mouse has both E and I cells |
| Groups | E vs I | As in the assignment |
| f-I curve | Mean rate over the 1-s step; plotted against pA and against pA above rheobase | Separates threshold shifts from gain changes |
| Comparison range | Statistics over rheobase to rheobase + 60 pA; plots over the full range, with that range marked | Every cell covers it, whatever its protocol |
| Representative cell | One recorded cell per group: the one closest to the group median rheobase and shared-range slope | A real cell, not an average |
| Model | LIF first, then a Hodgkin–Huxley-type model with Na/K/M channels | See the model sections |
| Perturbation | $g_L$ (and, in the HH-type model, $\bar g_{Na}$, $\bar g_K$) × 0.5, 0.75, 1, 1.5, 2; everything else fixed | Spans halving to doubling |
| Deliverable | Figures and [report.md](report.md), later a GitHub Pages page | No fixed format required |

## References

- Gouwens N. W. et al. (2020). Integrated morphoelectric and transcriptomic classification of cortical GABAergic cells. *Cell*.
- Sorensen S. A., Gouwens N. W. et al. Mouse excitatory visual cortex Patch-seq and whole-neuron morphology. Preprint: bioRxiv doi:10.1101/2023.11.25.568393. Code and data tables: [AllenInstitute/exc_vis_manuscript](https://github.com/AllenInstitute/exc_vis_manuscript). TODO: replace with the published version.
- Pospischil M. et al. (2008). Minimal Hodgkin–Huxley type models for different classes of cortical and thalamic neurons. *Biological Cybernetics*.
- Brette R. & Gerstner W. (2005). Adaptive exponential integrate-and-fire model as an effective description of neuronal activity. *Journal of Neurophysiology*.
- Holt G. R. & Koch C. (1997). Shunting inhibition does not have a divisive effect on firing rates. *Neural Computation*.
- Desai N. S., Rutherford L. C. & Turrigiano G. G. (1999). Plasticity in the intrinsic excitability of cortical pyramidal neurons. *Nature Neuroscience*.
