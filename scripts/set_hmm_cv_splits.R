#!/usr/bin/env Rscript

library("data.table")

set.seed(1)
n_folds <- 2

df <- fread("../samplesheets/trajectories_for_hmm.csv")[
  , .(id, cv_fold = sample(1:n_folds, .N, replace = TRUE))
]

fwrite(df, "../samplesheets/hmm_cv_splits.csv")