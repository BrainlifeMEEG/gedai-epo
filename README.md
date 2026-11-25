# GEDAI -- unsupervised EEG denoising (epoched data)

<!-- brainlife-app:begin (generated from the App manifest; edit the manifest, not this block) -->

> **Denoises the EEG channels of epoched data with GEDAI 0.4.2, an unsupervised leadfield-filtering method with automatic SENSAI thresholds, and outputs denoised epochs, the removed noise and a QC report.**

Contrasts each epoch's covariance with a template leadfield reference covariance and projects out non-brain components chosen by the SENSAI criterion; non-EEG channels pass through unchanged and the output is average-referenced. Use it for automated artifact removal in epoched EEG without manual ICA review. EEG channels should have standard 10-05 names or digitized positions and be high-pass filtered. GEDAI is licensed for non-commercial use only.

| | |
|---|---|
| 📥 **Inputs** | `neuro/meeg/mne/epochs` |
| 📤 **Outputs** | `neuro/meeg/mne/epochs` ×2 |
| 🐳 **Container** | `brainlifemeeg/mne-fsaverage:1.12.1` |
| ⚙️ **Resources** | 8 CPUs · 24 GB memory · 2 h walltime |
| 🔗 **Run it** | [doi.org/10.25663/brainlife.app.960](https://doi.org/10.25663/brainlife.app.960) |

#### 📥 Inputs

| Input | Datatype | Required | Description |
|---|---|---|---|
| `epo` | `neuro/meeg/mne/epochs` | yes | MNE epochs .fif with EEG channels |

#### 📤 Outputs

| Output | Datatype | Description |
|---|---|---|
| `out_dir` | `neuro/meeg/mne/epochs` | Denoised epochs |
| `out_noise` | `neuro/meeg/mne/epochs` | The signal GEDAI removed (EEG channels only) |

<details>
<summary><b>⚙️ Parameters</b> (12)</summary>

| Parameter | Type | Default | Description |
|---|---|---|---|
| `mode` | enum | `broadband` | GEDAI variant (epochs). Options: `broadband`, `multiband`. |
| `l_freq` | number |  | High-pass (Hz) before fitting. Empty/null = skip (a warning is issued if the epochs look unfiltered). |
| `noise_multiplier` | string | `3.0` | SENSAI noise weighting -- higher keeps more brain signal. A number > 0, or a preset: standard / auto (3.0), aggressive (1.5), conservative (6.0). |
| `sensai_method` | enum | `optimize` | SENSAI threshold search strategy. Options: `optimize`, `gridsearch`. |
| `broadband_pass` | boolean | `true` | multiband / adaptive: run an initial broadband GEDAI pre-pass. |
| `wavelet_low_cutoff` | string | `auto` | multiband / adaptive: wavelet bands below this (Hz) are zeroed. 'auto' uses the data high-pass. |
| `add_reference_channel` | string |  | Name of a flat reference channel to add before fitting (preserves rank for non-average references). Kept in the output. |
| `montage` | string |  | Name of an MNE built-in montage to attach before fitting (e.g. standard_1005, GSN-HydroCel-257). Needed when channel names are non-standard and no positions are stored. |
| `montage_handling` | enum | `auto` | How GEDAI's reference covariance is built. Options: `auto`, `leadfield`, `forward`, `channel_positions`, `drop_unmatched`. |
| `interpolate_bads` | boolean | `false` | Interpolate channels marked bad before fitting. |
| `output` | enum | `recombine` | What to write as the denoised output. Options: `recombine`, `eeg_only`. |
| `n_jobs` | number | `8` | Parallel jobs; match the allocated CPUs (currently 8). |

</details>

📚 **References used by this App:**
- Ros T, Ferat V, Huang Y, Colangelo C, Kia SM, Wolfers T, Vulliemoz S, Michela A. Return of the GEDAI: Unsupervised EEG Denoising based on Leadfield Filtering. bioRxiv (2025). [doi:10.1101/2025.10.04.680449](https://doi.org/10.1101/2025.10.04.680449)
- Hayashi S, Caron BA, Heinsfeld AS, et al. brainlife.io: a decentralized and open-source cloud platform to support neuroscience research. Nat Methods 21:809-813 (2024). [doi:10.1038/s41592-024-02237-2](https://doi.org/10.1038/s41592-024-02237-2)
- Gramfort A, et al. MEG and EEG data analysis with MNE-Python. Front Neurosci 7:267 (2013). [doi:10.3389/fnins.2013.00267](https://doi.org/10.3389/fnins.2013.00267)

<!-- brainlife-app:end -->

> ⚠️ **Licence / patent notice.** This app wraps
> [GEDAI](https://github.com/neurotuning/gedai), which is licensed under the
> **PolyForm Noncommercial License 1.0.0** and whose core algorithm is
> **patent pending**. The app as a whole may be used for **non-commercial
> purposes only** (research, education, evaluation). Commercial or proprietary
> use requires a separate licence (including a patent grant) from the GEDAI
> authors (tomas.ros@unige.ch). See `license.txt`.

## Description

GEDAI ("Generalized Eigenvalue De-Artifacting Instrument") is an **unsupervised**
EEG denoising method based on **leadfield filtering**. It contrasts the
covariance of each epoch with a physics-based reference covariance derived from a
template head model / leadfield (bundled with GEDAI, a ~343-electrode 10-05
montage). A generalized eigendecomposition keeps the components aligned with the
brain subspace and projects out the rest. The rejection threshold is selected
automatically by GEDAI's **SENSAI** criterion -- no manual component review, no
training data, no per-subject MRI or forward model.

This app runs GEDAI on the **EEG channels** of an MNE epochs `.fif` file. Non-EEG
channels (MEG, EOG, ECG, ...) are passed through untouched. For continuous data,
use the **`gedai`** app.

Reference: Ros T., Ferat V., Huang Y., Colangelo C., Kia S.M., Wolfers T.,
Vulliemoz S., Michela A. *Return of the GEDAI: Unsupervised EEG Denoising based
on Leadfield Filtering.* bioRxiv 2025. doi:10.1101/2025.10.04.680449

## Inputs

- **epo** (`neuro/meeg/mne/epochs`): an MNE epochs `.fif` file. Must contain EEG
  channels.

For best results the EEG channels should carry **standard 10-05 names**
(`Fp1`, `Fz`, `Cz`, ...) or have **digitised positions**, and the epochs should
be **high-pass filtered**.

## Outputs

- **out_dir/meg-epo.fif** (`neuro/meeg/mne/epochs`): the denoised epochs.
- **out_noise/meg-epo.fif** (`neuro/meeg/mne/epochs`): the *removed* signal (EEG
  channels only) -- `input EEG = denoised EEG + noise`.
- **out_report/report.html**: QC report (SENSAI fit, evoked GFP before/after,
  PSD before/after, ENOVA topomap).
- **out_figs/*.png**, **product.json**.

> GEDAI applies an **average reference** internally, so the denoised EEG output
> is average-referenced.

## Configuration parameters

| key | type | default | description |
|---|---|---|---|
| `mode` | `broadband` \| `multiband` | `broadband` | GEDAI variant. `broadband` = single pass (robust for any epoch length). `multiband` = per wavelet band, and needs epochs long enough for several wavelet levels (roughly a few seconds; short epochs will fail with a clear message). `adaptive` is **not** available for epochs -- use the `gedai` app. |
| `l_freq` | float \| null | `null` | High-pass (Hz) applied before fitting. `null` = skip (a warning is issued if the epochs look unfiltered). |
| `noise_multiplier` | float \| str | `3.0` | SENSAI noise weighting. Higher = keep more brain signal. A number > 0, or `standard` (3.0) / `aggressive` (1.5) / `conservative` (6.0). |
| `sensai_method` | `optimize` \| `gridsearch` | `optimize` | Threshold search strategy. |
| `broadband_pass` | bool | `true` | `multiband`: initial broadband pre-pass. |
| `wavelet_low_cutoff` | `auto` \| float \| null | `auto` | `multiband`: wavelet bands below this (Hz) are zeroed. `auto` uses the data high-pass. |
| `add_reference_channel` | str \| null | `null` | Add this flat reference channel before fitting (preserves rank for non-average references). Kept in the output. |
| `montage` | str \| null | `null` | MNE built-in montage to attach before fitting (e.g. `standard_1005`, `GSN-HydroCel-257`). |
| `montage_handling` | `auto` \| `leadfield` \| `forward` \| `channel_positions` \| `drop_unmatched` | `auto` | How GEDAI's reference covariance is built. `auto`: standard 10-05 names → bundled leadfield; else if every channel is positioned → an **fsaverage template forward solution** (+~35 s); else the 10-05-matched subset (leadfield) or positioned subset (forward). `forward` forces the fsaverage forward. `channel_positions` uses a fast distance covariance (**lower fidelity, tends to over-clean**). `leadfield` / `drop_unmatched` keep only 10-05-matched channels. |
| `interpolate_bads` | bool | `false` | Interpolate bad channels before fitting. |
| `output` | `recombine` \| `eeg_only` | `recombine` | `recombine`: full input copy with EEG replaced by the denoised signal. `eeg_only`: only the denoised EEG channels. |
| `n_jobs` | int | `8` | Parallel jobs. |

## Limitations

- **EEG only.** MEG channels are passed through unchanged.
- Needs standard 10-05 channel names **or** channel positions **or** a `montage`.
- `mode: adaptive` is not available -- use `gedai` on continuous data.
- The denoised EEG output is average-referenced.
- `mode: broadband` does not produce the SENSAI subspace metrics.

## Local testing

```bash
# edit config.json to point "epo" at a real epochs .fif, then:
pip install 'gedai==0.4.2'
# non-standard montages also need fsaverage: python -c "import mne; mne.datasets.fetch_fsaverage()"
python main.py
```

## Authors

- Maximilien Chaumon (https://github.com/dnacombo)
- GEDAI algorithm: T. Ros, V. Ferat et al. (https://github.com/neurotuning/gedai)

## Citations

If you use this app, please cite **GEDAI** and **brainlife.io**:

Ros T., Ferat V., Huang Y., Colangelo C., Kia S.M., Wolfers T., Vulliemoz S.,
Michela A. *Return of the GEDAI: Unsupervised EEG Denoising based on Leadfield
Filtering.* bioRxiv 2025. https://doi.org/10.1101/2025.10.04.680449

Hayashi, S., Caron, B.A., Heinsfeld, A.S. et al. brainlife.io: a decentralized
and open-source cloud platform to support neuroscience research. Nat Methods 21,
809-813 (2024). https://doi.org/10.1038/s41592-024-02237-2

Gramfort A. et al. MEG and EEG data analysis with MNE-Python. Front Neurosci 7,
267 (2013). https://doi.org/10.3389/fnins.2013.00267

## Funding Acknowledgement

brainlife.io is publicly funded. We kindly ask that you acknowledge the funding
below in your code and publications.

[![NSF-BCS-1734853](https://img.shields.io/badge/NSF_BCS-1734853-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1734853)
[![NSF-BCS-1636893](https://img.shields.io/badge/NSF_BCS-1636893-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1636893)
[![NSF-ACI-1916518](https://img.shields.io/badge/NSF_ACI-1916518-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1916518)
[![NSF-IIS-1912270](https://img.shields.io/badge/NSF_IIS-1912270-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1912270)
[![NIH-NIBIB-R01EB029272](https://img.shields.io/badge/NIH_NIBIB-R01EB029272-green.svg)](https://grantome.com/grant/NIH/R01-EB029272-01)
