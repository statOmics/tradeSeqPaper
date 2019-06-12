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


getWeightsBifurcation <- function(data, crv){
  miles <- data$milestone_ids
  startMile <- data$prior_information$start_milestones
  endMiles <- data$prior_information$end_milestones
  trueT <- data$prior_information$timecourse_continuous

  cellsInStartMile <- data$prior_information$groups_id$cell_id[data$prior_information$groups_id$group_id==startMile]
  cellsInEndMile1 <- data$prior_information$groups_id$cell_id[data$prior_information$groups_id$group_id==endMiles[1]]
  cellsInEndMile2 <- data$prior_information$groups_id$cell_id[data$prior_information$groups_id$group_id==endMiles[2]]
  endMileOfLongestLineage <- endMiles[which.max(c(max(trueT[cellsInEndMile1]), max(trueT[cellsInEndMile2])))]

  # fix weights for start and end milestones
  weights <- matrix(0, nrow=length(data$cell_ids), ncol=2)
  rownames(weights) <- data$cell_ids
  weights[cellsInStartMile,] <- c(1/2, 1/2)
  if(endMiles[1] == endMileOfLongestLineage){
    weights[cellsInEndMile1 , 1] <- 1
    weights[cellsInEndMile2 , 2] <- 1
  } else if(endMiles[2] == endMileOfLongestLineage){
    weights[cellsInEndMile2 , 1] <- 1
    weights[cellsInEndMile1 , 2] <- 1
  }

  # for the intermediate milestone, use tradeSeq assignments
  sWeights <- slingCurveWeights(crv)
  set.seed(8)
  tradeSeqAssignments <- .assignCells(sWeights)
  cellsToFill <- rowSums(weights)==0
  weights[cellsToFill,] <- tradeSeqAssignments[cellsToFill,]
  return(weights)
}


#### Monocle
# downloaded from https://github.com/cole-trapnell-lab/monocle-release on June 11, 2019.
source("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/methodsPaper/monocle-release/R/differential_expression.R")
source("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/methodsPaper/monocle-release/R/utils.R")
source("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/methodsPaper/monocle-release/R/expr_models.R")



BEAM_kvdb <- function (cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch",
    reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", branch_states = NULL,
    branch_point = 1, relative_expr = TRUE, branch_labels = NULL,
    verbose = FALSE, cores = 1, ...)
{
    branchTest_res <- branchTest_kvdb(cds, fullModelFormulaStr = fullModelFormulaStr,
        reducedModelFormulaStr = reducedModelFormulaStr, branch_states = branch_states,
        branch_point = branch_point, relative_expr = relative_expr,
        cores = cores, branch_labels = branch_labels, verbose = verbose,
        ...)
    cmbn_df <- branchTest_res[, 1:4]
    if (verbose)
        message("pass branchTest")
    fd <- fData(cds)[row.names(cmbn_df), ]
    cmbn_df <- cbind(cmbn_df, fd)
    if (verbose)
        message("return results")
    return(cmbn_df)
}


branchTest_kvdb <- function (cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch",
    reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", branch_states = NULL,
    branch_point = 1, relative_expr = TRUE, cores = 1, branch_labels = NULL,
    verbose = FALSE, ...)
{
  #### KVDB: Here, the Branch variable is created.
    #if ("Branch" %in% all.vars(terms(as.formula(fullModelFormulaStr)))) {
    #    cds_subset <- buildBranchCellDataSet(cds = cds, branch_states = branch_states,
    #        branch_point = branch_point, branch_labels = branch_labels,
    #        ...)
    #}
    #else cds_subset <- cds
    cds_subset <- cds
    branchTest_res <- differentialGeneTest_kvdb(cds_subset, fullModelFormulaStr = fullModelFormulaStr,
        reducedModelFormulaStr = reducedModelFormulaStr, cores = cores,
        relative_expr = relative_expr, verbose = verbose)
    return(branchTest_res)
}


buildBranchCellDataSet_kvdb <- function (cds, progenitor_method = c("sequential_split", "duplicate"),
    branch_states = NULL, branch_point = 1, branch_labels = NULL,
    stretch = TRUE)
{
    if (is.null(pData(cds)$State) | is.null(pData(cds)$Pseudotime))
        stop("Please first order the cells in pseudotime using orderCells()")
    if (is.null(branch_point) & is.null(branch_states))
        stop("Please either specify the branch_point or branch_states to select subset of cells")
    if (!is.null(branch_labels) & !is.null(branch_states)) {
        if (length(branch_labels) != length(branch_states))
            stop("length of branch_labels doesn't match with that of branch_states")
        branch_map <- setNames(branch_labels, as.character(branch_states))
    }
    if (cds@dim_reduce_type == "DDRTree") {
        pr_graph_cell_proj_mst <- minSpanningTree(cds)
    }
    else {
        pr_graph_cell_proj_mst <- cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree
    }
    root_cell <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
    root_state <- pData(cds)[root_cell, ]$State
    pr_graph_root <- subset(pData(cds), State == root_state)
    if (cds@dim_reduce_type == "DDRTree") {
        closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
        root_cell_point_in_Y <- closest_vertex[row.names(pr_graph_root),
            ]
    }
    else {
        root_cell_point_in_Y <- row.names(pr_graph_root)
    }
    root_cell <- names(which(degree(pr_graph_cell_proj_mst, v = root_cell_point_in_Y,
        mode = "all") == 1, useNames = T))[1]
    paths_to_root <- list()
    if (is.null(branch_states) == FALSE) {
        for (leaf_state in branch_states) {
            curr_cell <- subset(pData(cds), State == leaf_state)
            if (cds@dim_reduce_type == "DDRTree") {
                closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
                curr_cell_point_in_Y <- closest_vertex[row.names(curr_cell),
                  ]
            }
            else {
                curr_cell_point_in_Y <- row.names(curr_cell)
            }
            curr_cell <- names(which(degree(pr_graph_cell_proj_mst,
                v = curr_cell_point_in_Y, mode = "all") == 1,
                useNames = T))[1]
            path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst,
                curr_cell, root_cell)
            path_to_ancestor <- names(unlist(path_to_ancestor$vpath))
            if (cds@dim_reduce_type == "DDRTree") {
                closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
                ancestor_cells_for_branch <- row.names(closest_vertex)[which(V(pr_graph_cell_proj_mst)[closest_vertex]$name %in%
                  path_to_ancestor)]
            }
            else if (cds@dim_reduce_type == "ICA") {
                ancestor_cells_for_branch <- path_to_ancestor
            }
            ancestor_cells_for_branch <- intersect(ancestor_cells_for_branch,
                colnames(cds))
            paths_to_root[[as.character(leaf_state)]] <- ancestor_cells_for_branch
        }
    }
    else {
        if (cds@dim_reduce_type == "DDRTree")
            pr_graph_cell_proj_mst <- minSpanningTree(cds)
        else pr_graph_cell_proj_mst <- cds@auxOrderingData$ICA$cell_ordering_tree
        mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
        branch_cell <- mst_branch_nodes[branch_point]
        mst_no_branch_point <- pr_graph_cell_proj_mst - V(pr_graph_cell_proj_mst)[branch_cell]
        path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst,
            branch_cell, root_cell)
        path_to_ancestor <- names(unlist(path_to_ancestor$vpath))
        for (backbone_nei in V(pr_graph_cell_proj_mst)[suppressWarnings(nei(branch_cell))]$name) {
            descendents <- bfs(mst_no_branch_point, V(mst_no_branch_point)[backbone_nei],
                unreachable = FALSE)
            descendents <- descendents$order[!is.na(descendents$order)]
            descendents <- V(mst_no_branch_point)[descendents]$name
            if (root_cell %in% descendents == FALSE) {
                path_to_root <- unique(c(path_to_ancestor, branch_cell,
                  descendents))
                if (cds@dim_reduce_type == "DDRTree") {
                  closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
                  path_to_root <- row.names(closest_vertex)[which(V(pr_graph_cell_proj_mst)[closest_vertex]$name %in%
                    path_to_root)]
                }
                else {
                  path_to_root <- path_to_root
                }
                closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
                path_to_root <- intersect(path_to_root, colnames(cds))
                paths_to_root[[backbone_nei]] <- path_to_root
            }
        }
    }
    all_cells_in_subset <- c()
    if (is.null(branch_labels) == FALSE) {
        if (length(branch_labels) != 2)
            stop("Error: branch_labels must have exactly two entries")
        names(paths_to_root) <- branch_labels
    }
    for (path_to_ancestor in paths_to_root) {
        if (length(path_to_ancestor) == 0) {
            stop("Error: common ancestors between selected State values on path to root State")
        }
        all_cells_in_subset <- c(all_cells_in_subset, path_to_ancestor)
    }
    all_cells_in_subset <- unique(all_cells_in_subset)
    common_ancestor_cells <- intersect(paths_to_root[[1]], paths_to_root[[2]])
    cds <- cds[, row.names(pData(cds[, all_cells_in_subset]))]
    Pseudotime <- pData(cds)$Pseudotime
    pData <- pData(cds)
    if (stretch) {
        max_pseudotime <- -1
        for (path_to_ancestor in paths_to_root) {
            max_pseudotime_on_path <- max(pData[path_to_ancestor,
                ]$Pseudotime)
            if (max_pseudotime < max_pseudotime_on_path) {
                max_pseudotime <- max_pseudotime_on_path
            }
        }
        branch_pseudotime <- max(pData[common_ancestor_cells,
            ]$Pseudotime)
        for (path_to_ancestor in paths_to_root) {
            max_pseudotime_on_path <- max(pData[path_to_ancestor,
                ]$Pseudotime)
            path_scaling_factor <- (max_pseudotime - branch_pseudotime)/(max_pseudotime_on_path -
                branch_pseudotime)
            if (is.finite(path_scaling_factor)) {
                branch_cells <- setdiff(path_to_ancestor, common_ancestor_cells)
                pData[branch_cells, ]$Pseudotime <- ((pData[branch_cells,
                  ]$Pseudotime - branch_pseudotime) * path_scaling_factor +
                  branch_pseudotime)
            }
        }
        pData$Pseudotime <- 100 * pData$Pseudotime/max_pseudotime
    }
    pData$original_cell_id <- row.names(pData)
    pData$original_cell_id <- row.names(pData)
    if (length(paths_to_root) != 2)
        stop("more than 2 branch states are used!")
    pData[common_ancestor_cells, "Branch"] <- names(paths_to_root)[1]
    progenitor_pseudotime_order <- order(pData[common_ancestor_cells,
        "Pseudotime"])
    if (progenitor_method == "duplicate") {
        ancestor_exprs <- exprs(cds)[, common_ancestor_cells]
        expr_blocks <- list()
        for (i in 1:length(paths_to_root)) {
            if (nrow(ancestor_exprs) == 1)
                exprs_data <- t(as.matrix(ancestor_exprs))
            else exprs_data <- ancestor_exprs
            colnames(exprs_data) <- paste("duplicate", i, 1:length(common_ancestor_cells),
                sep = "_")
            expr_lineage_data <- exprs(cds)[, setdiff(paths_to_root[[i]],
                common_ancestor_cells)]
            exprs_data <- cbind(exprs_data, expr_lineage_data)
            expr_blocks[[i]] <- exprs_data
        }
        ancestor_pData_block <- pData[common_ancestor_cells,
            ]
        pData_blocks <- list()
        weight_vec <- c()
        for (i in 1:length(paths_to_root)) {
            weight_vec <- c(weight_vec, rep(1, length(common_ancestor_cells)))
            weight_vec_block <- rep(1, length(common_ancestor_cells))
            new_pData_block <- ancestor_pData_block
            row.names(new_pData_block) <- paste("duplicate",
                i, 1:length(common_ancestor_cells), sep = "_")
            pData_lineage_cells <- pData[setdiff(paths_to_root[[i]],
                common_ancestor_cells), ]
            weight_vec_block <- c(weight_vec_block, rep(1, nrow(pData_lineage_cells)))
            weight_vec <- c(weight_vec, weight_vec_block)
            new_pData_block <- rbind(new_pData_block, pData_lineage_cells)
            new_pData_block$Branch <- names(paths_to_root)[i]
            pData_blocks[[i]] <- new_pData_block
        }
        pData <- do.call(rbind, pData_blocks)
        exprs_data <- do.call(cbind, expr_blocks)
    }
    else if (progenitor_method == "sequential_split") {
        pData$Branch <- names(paths_to_root)[1]
        branchA <- progenitor_pseudotime_order[seq(1, length(common_ancestor_cells),
            by = 2)]
        pData[common_ancestor_cells[branchA], "Branch"] <- names(paths_to_root)[1]
        branchB <- progenitor_pseudotime_order[seq(2, length(common_ancestor_cells),
            by = 2)]
        pData[common_ancestor_cells[branchB], "Branch"] <- names(paths_to_root)[2]
        zero_pseudotime_root_cell <- common_ancestor_cells[progenitor_pseudotime_order[1]]
        exprs_data <- cbind(exprs(cds), duplicate_root = exprs(cds)[,
            zero_pseudotime_root_cell])
        pData <- rbind(pData, pData[zero_pseudotime_root_cell,
            ])
        row.names(pData)[nrow(pData)] <- "duplicate_root"
        pData[nrow(pData), "Branch"] <- names(paths_to_root)[2]
        weight_vec <- rep(1, nrow(pData))
        for (i in 1:length(paths_to_root)) {
            path_to_ancestor <- paths_to_root[[i]]
            branch_cells <- setdiff(path_to_ancestor, common_ancestor_cells)
            pData[branch_cells, ]$Branch <- names(paths_to_root)[i]
        }
    }
    pData$Branch <- as.factor(pData$Branch)
    pData$State <- factor(pData$State)
    Size_Factor <- pData$Size_Factor
    fData <- fData(cds)
    colnames(exprs_data) <- row.names(pData)
    cds_subset <- newCellDataSet(as.matrix(exprs_data), phenoData = new("AnnotatedDataFrame",
        data = pData), featureData = new("AnnotatedDataFrame",
        data = fData), expressionFamily = cds@expressionFamily,
        lowerDetectionLimit = cds@lowerDetectionLimit)
    pData(cds_subset)$State <- as.factor(pData(cds_subset)$State)
    pData(cds_subset)$Size_Factor <- Size_Factor
    cds_subset@dispFitInfo <- cds@dispFitInfo
    return(cds_subset)
}


















differentialGeneTest_kvdb <- function (cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)",
    reducedModelFormulaStr = "~1", relative_expr = TRUE, cores = 1,
    verbose = FALSE)
{
    status <- NA
    if (class(cds)[1] != "CellDataSet") {
        stop("Error cds is not of type 'CellDataSet'")
    }
    all_vars <- c(all.vars(formula(fullModelFormulaStr)), all.vars(formula(reducedModelFormulaStr)))
    pd <- pData(cds)
    for (i in all_vars) {
        x <- pd[, i]
        if (any((c(Inf, NaN, NA) %in% x))) {
            stop("Error: Inf, NaN, or NA values were located in pData of cds in columns mentioned in model terms")
        }
    }
    if (relative_expr && cds@expressionFamily@vfamily %in% c("negbinomial",
        "negbinomial.size")) {
        if (is.null(sizeFactors(cds)) || sum(is.na(sizeFactors(cds)))) {
            stop("Error: to call this function with relative_expr==TRUE, you must first call estimateSizeFactors() on the CellDataSet.")
        }
    }
    if (cores > 1) {
        diff_test_res <- mcesApply(cds, 1, diff_test_helper_kvdb,
            c("BiocGenerics", "VGAM", "Matrix"), cores = cores,
            fullModelFormulaStr = fullModelFormulaStr, reducedModelFormulaStr = reducedModelFormulaStr,
            expressionFamily = cds@expressionFamily, relative_expr = relative_expr,
            disp_func = cds@dispFitInfo[["blind"]]$disp_func,
            verbose = verbose)
        diff_test_res
    }
    else {
        diff_test_res <- smartEsApply(cds, 1, diff_test_helper_kvdb,
            convert_to_dense = TRUE, fullModelFormulaStr = fullModelFormulaStr,
            reducedModelFormulaStr = reducedModelFormulaStr,
            expressionFamily = cds@expressionFamily, relative_expr = relative_expr,
            disp_func = cds@dispFitInfo[["blind"]]$disp_func,
            verbose = verbose)
        diff_test_res
    }
    diff_test_res <- do.call(rbind.data.frame, diff_test_res)
    diff_test_res$qval <- 1
    diff_test_res$qval[which(diff_test_res$status == "OK")] <- p.adjust(subset(diff_test_res,
        status == "OK")[, "pval"], method = "BH")
    diff_test_res <- merge(diff_test_res, fData(cds), by = "row.names")
    row.names(diff_test_res) <- diff_test_res[, 1]
    diff_test_res[, 1] <- NULL
    diff_test_res[row.names(cds), ]
}


diff_test_helper_kvdb <- function(x,
                             fullModelFormulaStr,
                             reducedModelFormulaStr,
                             expressionFamily,
                             relative_expr,
                             weights,
                             disp_func=NULL,
                             verbose=FALSE
                             ){

  reducedModelFormulaStr <- paste("f_expression", reducedModelFormulaStr, sep="")
  fullModelFormulaStr <- paste("f_expression", fullModelFormulaStr, sep="")

  x_orig <- x
  disp_guess <- 0

  if (expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")){
    if (relative_expr == TRUE)
    {
      x <- x / Size_Factor
    }
    f_expression <- round(x)
    if (is.null(disp_func) == FALSE){
      disp_guess <- calculate_NB_dispersion_hint(disp_func, round(x_orig))
      if (is.null(disp_guess) == FALSE && disp_guess > 0 && is.na(disp_guess) == FALSE  ) {
        # FIXME: In theory, we could lose some user-provided parameters here
        # e.g. if users supply zero=NULL or something.
        if (expressionFamily@vfamily == "negbinomial")
          expressionFamily <- negbinomial(isize=1/disp_guess)
        else
          expressionFamily <- negbinomial.size(size=1/disp_guess)
      }
    }
  }else if (expressionFamily@vfamily %in% c("uninormal")){
    f_expression <- x
  }else if (expressionFamily@vfamily %in% c("binomialff")){
    f_expression <- x
    #f_expression[f_expression > 1] <- 1
  }else{
    f_expression <- log10(x)
  }

  test_res <- tryCatch({
    if (expressionFamily@vfamily %in% c("binomialff")){
      if (verbose){
        full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily)
        reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily)
      }else{
        full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily))
        reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily))
      }
    }else{
      if (verbose){
        full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily)
        reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily)
      }else{
        full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily))
        reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily))
      }
    }

    #print(full_model_fit)
    #print(coef(reduced_model_fit))
    compareModels(list(full_model_fit), list(reduced_model_fit))
  },
  #warning = function(w) { FM_fit },
  error = function(e) {
    if(verbose)
      print (e);
      data.frame(status = "FAIL", family=expressionFamily@vfamily, pval=1.0, qval=1.0)
    #data.frame(status = "FAIL", pval=1.0)
  }
  )
  test_res
}
