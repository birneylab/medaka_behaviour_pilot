#!/usr/bin/env Rscript

library("data.table")

set.seed(1)

df_raw <- fread("../samplesheets/trajectories_for_hmm.csv")
df_ref <- df_raw[, .(id = sprintf("%s_ref", id))]
df_test <- df_raw[, .(id = sprintf("%s_test", id))]
df <- rbind(df_ref, df_test)
n_samples <- nrow(df)
n_samples_a <- n_samples %/% 2
n_samples_b <- n_samples - n_samples_a
cv_fold <- c(rep("A", n_samples_a), rep("B", n_samples_b)) |> sample()
df[, cv_fold := ..cv_fold]

fwrite(df, "../samplesheets/hmm_cv_splits.csv")