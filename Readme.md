Condensation on dsRNA switches CaMV P6 from Target of Rapamycin
activation for virus translation to viroplasm service
================
Sanjana Holla\* 1, Aimer Gutiérrez-Díaz\* 1, Claudia Cortez 1, Emilie
Vantard 1, Eduardo Méndez-López 1, Athanasios Xhurxhi 2, Subhankar Sahu
2, Christophe Ritzenthaler 2 & Anders Hafrén+ 1

**affiliations:** 1. Department of Plant Biology, Uppsala BioCenter,
Swedish University of Agricultural Sciences and Linnean Center for Plant
Biology, Box 7080, 75007 Uppsala, Sweden. 2. Institut de Biologie
Moléculaire des Plantes, UPR2357 du Centre National de la Recherche
Scientifique, Université de Strasbourg, Strasbourg F-67084, France

- equal contribution
- **correspondence:** <anders.hafren@slu.se>

# Methods : Protein Complex Structure

P6 dimer complexes with either a 30bp dsRNA duplex or TOR were predicted
using Alphafold3 \[[1](#ref-abramson2024accurate)\] and Protenix
\[[2](#ref-bytedance2025protenix)\] with default parameters, retaining
per model ipTM, ptm, pLDDT, and GPDE confidence metrics to benchmark
confidence. The dsRNA sequence used was AAGCUUCAAAUUAAGUCAGCUCCUUAAAUG
and its reverse complement. From Protenix JSON output, Predicted Aligned
Error (PAE) and Predicted Distance Error (PDE) from all 5 models, and
the single-residue contact-probability (CP) matrix was extracted and
processed into per-residue level interface statistics by extracting the
minimal PAE or maximal contact of the Protein-Protein interaction
surface of P6-to-itself, P6 dimers, P6 to 30bp dsRNA duplex, or P6 bound
to TOR.

## CaMV P6 structural analysis

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

### Repository contents

| Directory  | Contents                                                                                                                             |
|------------|--------------------------------------------------------------------------------------------------------------------------------------|
| `Inputs/`  | Protenix inputs containing the sequences and complex stoichiometries.                                                                |
| `Models/`  | Ranked Protenix coordinate models and confidence summaries, representative manuscript models, and the AlphaFold 3 comparison models. |
| `Raw/`     | Compressed Protenix full-data JSON outputs used for metric extraction.                                                               |
| `Code/`    | Scripts for extracting PAE, PDE, pLDDT, and contact-probability information and regenerating the main interface plots.               |
| `Data/`    | Analysis-ready R objects used by the plotting workflow.                                                                              |
| `Results/` | Per-residue interface summaries and supplementary result tables.                                                                     |
| `Figures/` | Computational plots generated from the deposited analyses.                                                                           |

### Deposited complexes

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

### Requirements

- Python 3.10 or later for extracting matrices from Protenix JSON
  output.
- R 4.4 or later for the analysis and plotting scripts.
- R packages: `rmarkdown`, `ggplot2`, `dplyr`, `readr`, `tidyr`,
  `scales`, and `svglite`.

Create the R environment with:

``` r
install.packages(c("rmarkdown", "ggplot2", "dplyr", "readr", "tidyr", "scales", "svglite"))
```

### Reproducing metric extraction

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

### Interpretation

PAE and PDE report model confidence rather than experimental binding
affinity. Contact probability identifies model-supported contacts but
does not establish that an interaction occurs in vivo. The deposited
models should therefore be interpreted together with the experimental
condensation and translation assays reported in the manuscript.

# Methods : Disorder Prediction

Intrinsic disorder and conditional folding propensity were estimated
using
[AlphaFold-Disorder](https://github.com/BioComputingUP/AlphaFold-disorder)
\[[3](#ref-piovesan2022intrinsic)\].

# License

Code is released under the GNU General Public License v3.0. The
manuscript, figures, model outputs, and third-party software remain
subject to their respective licences and citation requirements.

# References

<div id="refs" class="references csl-bib-body">

<div id="ref-abramson2024accurate" class="csl-entry">

1\. Abramson J, Adler J, Dunger J, Evans R, Green T, Pritzel A, et al.
Accurate structure prediction of biomolecular interactions with
AlphaFold 3. Nature. 2024;630:493–500.

</div>

<div id="ref-bytedance2025protenix" class="csl-entry">

2\. Team BAA, Chen X, Zhang Y, Lu C, Ma W, Guan J, et al.
Protenix-advancing structure prediction through a comprehensive
AlphaFold3 reproduction. BioRxiv. 2025;2025–01.

</div>

<div id="ref-piovesan2022intrinsic" class="csl-entry">

3\. Piovesan D, Monzon AM, Tosatto SC. Intrinsic protein disorder and
conditional folding in AlphaFoldDB. Protein Science. 2022;31:e4466.

</div>

</div>
