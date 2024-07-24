# Measurement and classification of bold-shy behaviours in medaka fish.

This repository describes the analysis done in this work: [DOI of preprint].
We studied differences in behaviour among medaka (*O. latipes*) strains.
We collected videos of medaka fish pairs, tracked their movements with the [idtracker.ai](https://idtracker.ai/latest/) package, and used a Hidden Markov Model (HMM) to classify behavioural modes.
We detected significant differences in behaviour among medaka strains, and also differences in the behaviour of the tank partner (always from the same strain, iCab) depending on the medaka strain they are paired with.

This is a [Nextflow](https://www.nextflow.io/) pipeline that can reproduce the whole tracking and HMM training and optimisation. The original videos used in the analysis are available at [DOI of data]. Be sure to set the correct paths in the [params.yaml](params.yaml) and [samplesheet](samplesheet) files for your system when reproducing this analysis.

The figures and tables that we present in the manuscript can be reproduced from the respective [Jupyter](https://jupyter.org/) notebooks present in the [figures_and_tables](figures_and_tables) folder. Also here, make sure to set the correct paths for your system. Pre-computed partial and final results are available under [results](results) the folder.