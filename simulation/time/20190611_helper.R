.assignCells <- function(cellWeights) {
  if (is.null(dim(cellWeights))) {
    if (any(cellWeights == 0)) {
      stop("Some cells have no positive cell weights.")
    } else {
      return(matrix(1, nrow = length(cellWeights), ncol = 1))
    }
  } else {
    if (any(rowSums(cellWeights) == 0)) {
      stop("Some cells have no positive cell weights.")
    } else {
      # normalize weights
      normWeights <- sweep(cellWeights, 1,
        FUN = "/",
        STATS = apply(cellWeights, 1, sum)
      )
      # sample weights
      wSamp <- apply(normWeights, 1, function(prob) {
        rmultinom(n = 1, prob = prob, size = 1)
      })
      # If there is only one lineage, wSamp is a vector so we need to adjust for that
      if (is.null(dim(wSamp))) {
        wSamp <- matrix(wSamp, ncol = 1)
      } else {
        wSamp <- t(wSamp)
      }
      return(wSamp)
    }
  }
}


getWeightsBifurcation <- function(data, crv) {
  miles <- data$milestone_ids
  startMile <- data$prior_information$start_milestones
  endMiles <- data$prior_information$end_milestones
  trueT <- data$prior_information$timecourse_continuous

  cellsInStartMile <- data$prior_information$groups_id$cell_id[
    data$prior_information$groups_id$group_id == startMile]
  cellsInEndMile1 <- data$prior_information$groups_id$cell_id[
    data$prior_information$groups_id$group_id == endMiles[1]]
  cellsInEndMile2 <- data$prior_information$groups_id$cell_id[
    data$prior_information$groups_id$group_id == endMiles[2]]
  endMileOfLongestLineage <- endMiles[which.max(c(max(trueT[
    cellsInEndMile1]), max(trueT[cellsInEndMile2])))]

  # fix weights for start and end milestones
  weights <- matrix(0, nrow = length(data$cell_ids), ncol = 2)
  rownames(weights) <- data$cell_ids
  weights[cellsInStartMile, ] <- c(1 / 2, 1 / 2)
  if (endMiles[1] == endMileOfLongestLineage) {
    weights[cellsInEndMile1, 1] <- 1
    weights[cellsInEndMile2, 2] <- 1
  } else if (endMiles[2] == endMileOfLongestLineage) {
    weights[cellsInEndMile2, 1] <- 1
    weights[cellsInEndMile1, 2] <- 1
  }

  # for the intermediate milestone, use tradeSeq assignments
  sWeights <- slingCurveWeights(crv)
  set.seed(8)
  tradeSeqAssignments <- .assignCells(sWeights)
  cellsToFill <- rowSums(weights) == 0
  weights[cellsToFill, ] <- tradeSeqAssignments[cellsToFill, ]
  return(weights)
}


#### Monocle
# downloaded from https://github.com/cole-trapnell-lab/monocle-release on June 11, 2019.
# source(
#   "~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/methodsPaper/monocle-release/R/differential_expression.R"
#   )
# source(
#   "~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/methodsPaper/monocle-release/R/utils.R"
#   )
# source(
#   "~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/methodsPaper/monocle-release/R/expr_models.R"
#   )

BEAM_kvdb <- function(cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch",
                      reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)",
                      branch_states = NULL, branch_point = 1,
                      relative_expr = TRUE, branch_labels = NULL,
                      verbose = FALSE, cores = 1, ...) {
  branchTest_res <- branchTest_kvdb(cds,
    fullModelFormulaStr = fullModelFormulaStr,
    reducedModelFormulaStr = reducedModelFormulaStr, branch_states = branch_states,
    branch_point = branch_point, relative_expr = relative_expr,
    cores = cores, branch_labels = branch_labels, verbose = verbose,
    ...
  )
  cmbn_df <- branchTest_res[, 1:4]
  if (verbose) {
    message("pass branchTest")
  }
  fd <- fData(cds)[row.names(cmbn_df), ]
  cmbn_df <- cbind(cmbn_df, fd)
  if (verbose) {
    message("return results")
  }
  return(cmbn_df)
}


branchTest_kvdb <- function(cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch",
                            reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)",
                            branch_states = NULL,
                            branch_point = 1, relative_expr = TRUE, cores = 1,
                            branch_labels = NULL, verbose = FALSE, ...) {
  #### KVDB: Here, the Branch variable is created.
  # if ("Branch" %in% all.vars(terms(as.formula(fullModelFormulaStr)))) {
  #    cds_subset <- buildBranchCellDataSet(cds = cds, branch_states = branch_states,
  #        branch_point = branch_point, branch_labels = branch_labels,
  #        ...)
  # }
  # else cds_subset <- cds
  cds_subset <- cds
  branchTest_res <- differentialGeneTest(cds_subset,
    fullModelFormulaStr = fullModelFormulaStr,
    reducedModelFormulaStr = reducedModelFormulaStr, cores = cores,
    relative_expr = relative_expr, verbose = verbose
  )
  return(branchTest_res)
}