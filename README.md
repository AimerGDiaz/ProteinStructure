Aimer G. Diaz

# CaMV P6 structural analysis

This repository contains the structural models, prediction inputs,
analysis code, processed results, and computational figures supporting
the analysis of the *Cauliflower mosaic virus* (CaMV) P6 protein in the
manuscript:

> Condensation on dsRNA switches CaMV P6 from Target of Rapamycin
> activation for virus translation to viroplasm services.

The analysis compares P6 homodimers predicted in complex with a 30-bp
dsRNA duplex or with *Arabidopsis thaliana* TARGET OF RAPAMYCIN (TOR).
It also includes the P6IDR construct in which the conditionally folding
region at residues 224-247 was modified, modelled as a dimer with the
same dsRNA duplex. Models were generated with Protenix using the
prediction inputs and seeds provided here. AlphaFold 3 models used as
independent structural comparisons are deposited under
`Models/AlphaFold3/`.

## Repository contents

| Directory  | Contents                                                                                                                             |
|------------|--------------------------------------------------------------------------------------------------------------------------------------|
| `Inputs/`  | Protenix inputs containing the sequences and complex stoichiometries.                                                                |
| `Models/`  | Ranked Protenix coordinate models and confidence summaries, representative manuscript models, and the AlphaFold 3 comparison models. |
| `Raw/`     | Compressed Protenix full-data JSON outputs used for metric extraction.                                                               |
| `Code/`    | Scripts for extracting PAE, PDE, pLDDT, and contact-probability information and regenerating the main interface plots.               |
| `Data/`    | Analysis-ready R objects used by the plotting workflow.                                                                              |
| `Results/` | Per-residue interface summaries and supplementary result tables.                                                                     |
| `Figures/` | Computational plots generated from the deposited analyses.                                                                           |

## Deposited complexes

| Complex     | Composition                                                     | Protenix seed | Samples | Manuscript model                            |
|-------------|-----------------------------------------------------------------|--------------:|--------:|---------------------------------------------|
| P6-dsRNA    | Two 520-aa P6 chains and two complementary 30-nt RNA strands    |         95890 |       5 | `Models/Model_S1_P6-dsRNA.cif.gz`           |
| P6-TOR      | Two 520-aa P6 chains and one 2,481-aa *A. thaliana* TOR chain   |         88467 |       5 | `Models/Model_S2_P6-TOR.cif.gz`             |
| P6IDR-dsRNA | Two 509-aa P6IDR chains and two complementary 30-nt RNA strands |         95890 |       5 | `Models/Representative_P6-IDR-dsRNA.cif.gz` |

The supplementary wild-type model files are copies of sample 0 from the
corresponding five-model prediction sets.

In the P6IDR modelled construct, wild-type residues 224-247
(`WLTLGTKRPSSDPAPKEISFAPEI`) are replaced by the 13-residue sequence
`QTGVAYIPGAKCG`, producing a 509-aa chain. Relative to the deposited
wild-type input, the modelled construct also contains V142L.

## Requirements

- Python 3.10 or later for extracting matrices from Protenix JSON
  output.
- R 4.4 or later for the analysis and plotting scripts.
- R packages: `rmarkdown`, `ggplot2`, `dplyr`, `readr`, `tidyr`,
  `scales`, and `svglite`.

Create the R environment with:

``` r
install.packages(c("rmarkdown", "ggplot2", "dplyr", "readr", "tidyr", "scales", "svglite"))
```

## Reproducing metric extraction

The compressed full-data JSON files contain the model-wise PAE, PDE,
pLDDT, and contact-probability outputs. Extract the pairwise matrices
with:

``` bash
python Code/extract_protenix_metrics.py \
  --input-dir Raw/P6-dsRNA \
  --prefix P6_dimer_dsRNA \
  --output-dir Data/Matrices

python Code/extract_protenix_metrics.py \
  --input-dir Raw/P6-TOR \
  --prefix P6_dimer_TOR \
  --output-dir Data/Matrices

python Code/extract_protenix_metrics.py \
  --input-dir Raw/P6-IDR-dsRNA \
  --prefix P6_dimer_ivIDR_dsRNA \
  --output-dir Data/Matrices
```

Each command writes compressed CSV matrices for `token_pair_pae`,
`token_pair_pde`, and `contact_probs`, together with a table of
atom-level pLDDT values. Processing the P6-TOR full-data files requires
substantial RAM.

## Regenerating interface plots

The analysis-ready objects in `Data/` contain the per-residue values
used for the manuscript plots. From the repository root, run:

``` bash
Rscript Code/plot_p6_interfaces.R
```

The script writes regenerated SVG files to `Figures/Reproduced/`. The
publication versions are retained in `Figures/` for direct comparison.

## Interpretation

PAE and PDE report model confidence rather than experimental binding
affinity. Contact probability identifies model-supported contacts but
does not establish that an interaction occurs in vivo. The deposited
models should therefore be interpreted together with the experimental
condensation and translation assays reported in the manuscript.

## Citations

Please cite the accompanying manuscript

## License

Code is released under the GNU General Public License v3.0. The
manuscript, figures, model outputs, and third-party software remain
subject to their respective licences and citation requirements.
