#!/usr/bin/Rscript
source("main.R")
library(testthat)

csv <- paste0("data/verse_counts.tsv")

test_that("test data is loading correctly and reduced", {
  test_df <- load_n_trim(csv)
  expect_equal(names(test_df), c("vP0_1", "vP0_2", "vAd_1", "vAd_2"))
  expect_equal(dim(test_df), c(55416, 4))
  expect_equal(class(test_df), "data.frame")
})

test_that("deseq2 function is returning correct results", {
  load("data/mock_counts_df.RData") # loads the counts_df object into env
  coldata <- data.frame(condition = rep(c("day4", "day7"), each=2),
                        type="paired-end")
  row.names(coldata) <- c("vP4_1", "vP4_2", "vP7_1", "vP7_2")
  expect_warning(deseq <- run_deseq(counts_df, 
                                    coldata, 
                                    10, "condition_day7_vs_day4"))
  expect_equal(dim(deseq), c(19127, 6))
  expect_equal(class(deseq)[1], "DESeqResults")
  expect_equal(c("pvalue", "padj") %in% names(deseq), c(TRUE, TRUE))
})

test_that("edger function is returning correct results", {
  load("data/mock_counts_df.RData")
  group <- factor(rep(c(1,2), each=2))
  edger_res <- run_edger(counts_df, group)
  expect_equal(names(edger_res), c("logFC", "logCPM", "PValue"))
  expect_equal(dim(edger_res), c(15026, 3))
})

test_that("test limma + voom work + work together", {
  load("data/mock_counts_df.RData")
  group <- factor(rep(c(1,2), each=2))
  # design
  design <- data.frame(day4=1, day4vsday7=c(0, 0, 1, 1))
  row.names(design) <- c("vP4_1", "vP4_2", "vP7_1", "vP7_2")
  # limma + voom
  expect_warning(voom_res <- run_limma(counts_df, design, group)) # voom yelling, ignore
  expect_equal(dim(voom_res), c(15026, 6))
  expect_equal(names(voom_res), c("logFC", "AveExpr", "t", "P.Value", 
                                  "adj.P.Val", "B"))
})

test_that("ggplot and data formatting works", {
  deseq <- data.frame(pvalue = runif(1000, 1e-100, 1e-5), 
                      log2FoldChange = runif(1000, -9, 9),
                      padj = runif(1000, 1e-300, 1e-15))
  edger <- data.frame(PValue = runif(1000, 1e-100, 1e-5),
                      logFC = runif(1000, -9, 9),
                      padj = runif(1000, 1e-100, 1e-14))
  limma <- data.frame(P.Value = runif(1000, 1e-100, 1e-5),
                      logFC = runif(1000, -9, 9),
                      adj.P.Val = runif(1000, 1e-6, 1e-4))
  # combine_pval
  res <- combine_pval(deseq, edger, limma) 
  expect_equal(dim(res), c(3000, 2))
  expect_true(mean(res$pval) < 5e-5 && mean(res$pval) > 5e-7)
  expect_equal(table(res$package)[[1]], 1000)
  # create_facets
  facets <- create_facets(deseq, edger, limma)
  expect_equal(dim(facets), c(3000, 3))
  expect_equal(table(facets$package)[[2]], 1000)
  # theme_plot
  tibble(logFC = runif(3000, -9, 9),
         padj = runif(3000, 1e-300, 1e-5),
         package = rep(c("DESeq2", "edgeR", "Limma"), each=1000)) %>%
  theme_plot() -> volc_plot
  expect_equal(dim(volc_plot$data), c(3000, 3))
  expect_equal(class(volc_plot$layers[[1]]$geom)[1], "GeomPoint")
  expect_equal(class(volc_plot$facet$params), "list")
})

