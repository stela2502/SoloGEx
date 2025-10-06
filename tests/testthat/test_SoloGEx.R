library(testthat)
library(SoloGEx)

context("SoloGEx S4 basic workflow tests")

test_that("Constructor works and checks dimensions", {
  expr <- matrix(1:10, nrow = 5)
  colData <- data.frame(SampleID = paste0("S", 1:3))
  
  # Should fail (columns != rows)
  expect_error(SoloGEx(expr, colData))
  
  # Fix dimensions
  colData <- as.data.frame(matrix(1:10, ncol = 5, nrow = 2))
  edat = SoloGEx( expr, colData )
  expect_s4_class(edat, "SoloGEx")
  expect_equal(ncol(edat@expr), nrow(edat@colData))
})

test_that("compute_group_stats computes correct dimensions", {
  expr <- matrix(rnorm(10*4), nrow = 10)
  colData <- data.frame(SampleID = paste0("S", 1:4),
                        BiologicalGroup = rep(c("A","B"), each = 2))
  edat <- SoloGEx(expr, colData)
  edat <- compute_group_stats(edat, "BiologicalGroup")
  
  expect_equal(nrow(edat), nrow(expr))
  expect_equal(ncol(edat), 2*3) # 2 groups * (mean, sd, cv)
})

test_that("analyze_singletons works and produces expected columns", {
  # Reference
  ref_expr <- matrix(rnorm(10*4), nrow = 10)
  rownames(ref_expr) = paste(sep=".", "Gene", 1:10);
  ref_colData <- data.frame(SampleID = paste0("R", 1:4),
                            BiologicalGroup = rep(c("A","B"), each = 2))
  ref_edat <- SoloGEx(ref_expr, ref_colData)
  
  # Singleton
  singleton_expr <- matrix(rnorm(10*2), nrow = 10)
  rownames(singleton_expr) = paste(sep=".", "Gene", 1:10);
  singleton_colData <- data.frame(SampleID = paste0("S", 1:2))
  singleton_edat <- SoloGEx(singleton_expr, singleton_colData)
  
  singleton_edat <- analyze_singletons(singleton_edat, ref_edat, group_col="BiologicalGroup")
  
  expect_true(all(c("z_scores", "closest_group", "overall_dev") %in% names(singleton_edat@singleton_analysis)))
})

test_that("export_combined produces correct dimensions", {
  expr <- matrix(rnorm(10*2), nrow = 10)
  colData <- data.frame(SampleID = paste0("S", 1:2))
  edat <- SoloGEx(expr, colData)
  
  # Add dummy group_stats and singleton_analysis
  edat@group_stats <- matrix(rnorm(10*3), nrow=10)
  edat@singleton_analysis <- data.frame(Z_A=rnorm(10), closest_group=rep("A",10), overall_dev=rnorm(10))
  
  combined <- export_combined(edat)
  expect_true(is.data.frame(combined))
  expect_equal(nrow(combined), 10)
  expect_true(all(c("closest_group", "overall_dev") %in% colnames(combined)))
})
