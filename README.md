# GEDAI -- unsupervised EEG denoising (epoched data)

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

- **epo** (`neuro/mne`, `epo`): an MNE epochs `.fif` file. Must contain EEG
  channels.

For best results the EEG channels should carry **standard 10-05 names**
(`Fp1`, `Fz`, `Cz`, ...) or have **digitised positions**, and the epochs should
be **high-pass filtered**.

## Outputs

- **out_dir/epo.fif** (`neuro/mne`, `epo`): the denoised epochs.
- **out_dir/noise_epo.fif** (`neuro/mne`, `epo`): the *removed* signal (EEG
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
| `montage` | str \| null | `null` | MNE built-in montage to attach before fitting (e.g. `standard_1005`). |
| `montage_handling` | `auto` \| `leadfield` \| `channel_positions` \| `drop_unmatched` | `auto` | How to get the reference covariance when some channel names are not in the 10-05 leadfield. `auto`: leadfield → channel-position covariance → drop unmatched. |
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
