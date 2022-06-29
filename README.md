Repository to perform RNAseq analysis in sea urchin embryos between ethanol treated and control embryos.

# Installation

To install the pipeline, you will need to have the following tools installed in
your environemnt:
1. Snakemake
2. Python 3
3. Python libraries `pandas` and `numpy`
4. R
5. R packages `sva`, `gsva`, `deseq2`, and `ggplot2`

# Data availability
The data is available on GEO at GSE207100

# Running the pipeline:
To run the pipeline, adust the necessary file paths in `config.yaml` to match
your file system, then simply navigate to this repository and issue the command:

```bash
snakemake -j<number_of_jobs_to_run>
```
