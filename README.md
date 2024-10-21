# Measurement and classification of bold-shy behaviours in medaka fish

This repository describes the analysis done in this work:

> Measurement and classification of bold-shy behaviours in medaka fish
>
> Saul Pierotti, Ian Brettell, Tomas Fitzgerald, Cathrin Herder, Narendar Aadepu, Christian Pylatiuk, Joachim Wittbrodt, Ewan Birney, Felix Loosli
> **bioR$\chi$iv** 2024.10.18.618696; doi: [https://doi.org/10.1101/2024.10.18.618696](https://doi.org/10.1101/2024.10.18.618696)

We studied differences in behaviour among medaka (_O. latipes_) strains.
We collected videos of medaka fish pairs, tracked their movements with the [idtracker.ai](https://idtracker.ai/latest/) package, and used a Hidden Markov Model (HMM) to classify behavioural modes.
We detected significant differences in behaviour among medaka strains, and also differences in the behaviour of the tank partner (always from the same strain, iCab) depending on the medaka strain they are paired with.

This is a [Nextflow](https://www.nextflow.io/) pipeline that can reproduce the whole tracking and HMM training and optimisation. The original videos used in the analysis are available on the [EBI Bioimage Archive](https://doi.org/10.6019/S-BIAD1421). Be sure to correctly install Nextflow and [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) before running the pipeline. Set the correct paths in the [params.yaml](params.yaml) and [samplesheets](samplesheets) files for your system when reproducing this analysis.

The command to run the pipeline is the following and must be run within the cloned repository:

```
nextflow run main.nf -params-file params.yaml
```

The figures and tables that we present in the manuscript can be reproduced from the respective [Jupyter](https://jupyter.org/) notebooks present in the [figures_and_tables](figures_and_tables) folder. Also here, make sure to set the correct paths for your system. Pre-computed partial and final results are available under the [results](results) folder.
