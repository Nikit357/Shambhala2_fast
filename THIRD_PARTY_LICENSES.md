# Third-party components

The MIT license in `LICENSE` covers the original code in this repository. The following
vendored files originate elsewhere and remain under their upstream terms.

## Octave sources

| File | Origin | Upstream terms |
|---|---|---|
| `octave/CuBlock.m` | CuBlock — Junet V, Farrés J, Mas JM, Daura X. *CuBlock: a cross-platform normalization method for gene-expression microarrays*, Bioinformatics (2021). Vendored verbatim. | *(to be filled from the upstream distribution)* |
| `octave/CuBlock_fixed.m` | Modification of `CuBlock.m` above: k-means replaced by precomputed cluster labels. | *(to be filled — follows CuBlock)* |
| `octave/Shambhala2_piped.m` | Derived from `Shambhala2.m` — Borisov N, et al., doi:10.1093/bioinformatics/btab105, Zenodo record 6415067. Three changes: stdin input, stdout output, arguments injected via `--eval`. | *(to be filled)* |
| `octave/Shambhala2_piped_fixed.m`, `octave/Shambhala2_piped_preqn.m` | Further variants of `Shambhala2_piped.m` above. | *(to be filled — follows Shambhala2)* |
| `octave/readExpressionData.m` | From the Shambhala2 distribution, vendored verbatim. | *(to be filled — follows Shambhala2)* |
| `octave/quantilenorm.m`, `octave/kmeans.m` | Dependency-free Octave reimplementations written for this project, replacing functions unavailable in Conda Octave. | MIT (this repository) |

None of the Octave files carries a license header. The upstream terms must be established from
the original distributions and recorded above before this repository is published.

## Data

`Calibration_datasets/*.csv.gz` are previously published reference datasets, redistributed here
for reproducibility:

| File | Source |
|---|---|
| `P0_standard.csv.gz`, `Q0_standard.csv.gz` | Shambhala2 calibration and definitive references (Zenodo 6415067). Sample identifiers are public GEO `GSM…` and GTEx accessions. |
| `ANTE.csv.gz` | Normal-tissue reference. *(provenance to be filled)* |
| `OncoboxCancer.csv.gz` | Cancer-sample reference. *(provenance to be filled)* |

`tests/fixtures/*.csv` are derived from the Shambhala2 Zenodo record, transposed to the
rows-as-samples convention this pipeline uses.

## Python dependencies

Installed from PyPI at their own licenses; see `requirements.txt`. None is vendored into this
repository.
