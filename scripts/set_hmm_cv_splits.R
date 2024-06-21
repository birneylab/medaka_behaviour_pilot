#!/usr/bin/env Rscript

library("data.table")

set.seed(1)

df <- fread("../samplesheets/trajectories_for_hmm.csv")
n_samples <- nrow(df)
n_samples_a <- n_samples %/% 2
n_samples_b <- n_samples - n_samples_a
cv_fold <- c(rep("A", n_samples_a), rep("B", n_samples_b)) |> sample()
df <- df[, .(id, cv_fold = ..cv_fold)]

fwrite(df, "../samplesheets/hmm_cv_splits.csv")