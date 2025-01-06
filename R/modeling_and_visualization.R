#' Plot Overlap Group Lasso Feature Tree
#'
#' Visualizes the feature hierarchy from a group lasso model using a dendrogram.
#' Highlights features based on their beta coefficients, with different colors
#' representing the internal and leaf nodes.
#' @param modelFit Fitted group lasso model.
#' @param trainingData Data frame used for training the model.
#' @param treeData Data containing hierarchical structure of features.
#' @param nodesToRemove Vector of nodes to remove from the hierarchy.
#' @return A ggplot object representing the feature tree.
#' @importFrom dplyr mutate filter full_join
#' @importFrom tidyr pivot_wider
#' @importFrom ggtree ggtree geom_tiplab geom_nodepoint
#' @export
plotOverlapTree <- function(modelFit, trainingData, treeData, nodesToRemove) {
  # Function to sort each sublist
  sortSublist <- function(x) {
    return(sort(x))
  }

  # Predict the variables from the model
  predictedVars <- predict(modelFit, type = "vars", latent = T)

  # Extract coefficients and filter the relevant features
  resultData <- coef(modelFit, latent = T) %>%
    as.matrix() %>%
    as.data.frame() %>%
    dplyr::mutate(feature = rownames(.)) %>%
    dplyr::filter(feature %like% "_logit") %>%
    dplyr::filter(!feature %like% "Intercept") %>%
    dplyr::rename("coefficientValue" = V1) %>%
    dplyr::mutate(node = as.numeric(sub("grp([0-9]+).*", "\\1", feature))) %>%
    dplyr::filter(!is.na(node)) %>% # Exclude rows where `node` is NA
    dplyr::group_by(node) %>%
    dplyr::summarise(betaValues = sum(coefficientValue, na.rm = TRUE), .groups = "drop")

  # colnames(coefficients)[1] <- 'coefficientValue'

  # Initialize lists for groups, clusters, and beta values
  groupList <- list()
  clustersList <- list()
  betaValues <- list()

  # Modify tree data
  treeClusterData <- treeData$data
  treeClusterData$clustersList <- lapply(treeClusterData$clusters, function(x) gsub("_", " ", x))
  treeClusterData$clustersList <- lapply(treeClusterData$clustersList, function(x) gsub(" logit", "", x))

  # Function to remove specified nodes from each list
  removeNodes <- function(x, nodes) {
    return(setdiff(x, nodes))
  }

  # Apply the node removal and sorting to clusters
  treeClusterData$clustersList <- lapply(treeClusterData$clustersList, removeNodes, nodes = nodesToRemove)
  treeClusterData$clustersList <- lapply(treeClusterData$clustersList, sortSublist)

  # If we have an empty model
  if (nrow(resultData) < 1) {
    message("No non-zero proportions.. Will generate null tree")
    resultData <- treeClusterData
    resultData$betaInternal <- rep(0, nrow(resultData))
    resultData$betaLeaf <- rep(0, nrow(resultData))
  } else {
    # Merge results with tree data and compute additional metrics
    resultData <- resultData %>%
      dplyr::full_join(treeClusterData, by = "node") %>%
      dplyr::mutate(betaInternal = if_else(isTip == FALSE, betaValues, 0)) %>%
      dplyr::mutate(betaLeaf = if_else(isTip == TRUE, betaValues, 0)) %>%
      mutate_at(vars(betaValues, betaLeaf, betaInternal), ~ replace_na(., 0))
  }

  resultData$label <- str_replace(resultData$label, paste0("_", "logit"), "")
  resultData$label <- str_replace(resultData$label, "_", " ")

  # Create ggplot object for dendrogram visualization
  plot <- ggtree(resultData, layout = "dendrogram", hang = 0) +
    geom_tiplab(as_ylab = TRUE, geom = "text", size = 12, color = "black") +
    geom_nodepoint(aes(subset = betaInternal != 0, color = betaInternal), size = 10) +
    geom_nodepoint(aes(color = betaInternal), size = 0) +
    scale_color_gradient2(
      mid = "grey", high = muted("purple"),
      midpoint = 0, name = "Beta Internal",
      guide = guide_colorbar(order = 2)
    ) +
    new_scale_color() +
    geom_tippoint(aes(subset = betaLeaf != 0, color = betaLeaf), size = 10) +
    scale_color_gradient2(
      low = muted("green"), high = muted("orange"), mid = "grey",
      midpoint = 0, name = "Beta Leaf",
      guide = guide_colorbar(order = 3)
    ) +
    geom_tippoint(aes(color = betaLeaf), size = 0)

  return(plot)
}


#' Visualize Group Lasso Model Tree
#'
#' Creates a dendrogram visualization for the group lasso model. Optionally,
#' displays a heatmap overlay of the feature coefficients.
#' @param fit Fitted group lasso model.
#' @param tree Tree structure data used to define feature hierarchies.
#' @param type Character; specifies type of grouping, either "cluster" or "marker", default is "cluster".
#' @param heatmap Logical; if TRUE, overlays a heatmap of feature coefficients, default is TRUE.
#' @param trainingData Data frame of training data for use in the visualization.
#' @param nodesToRemove Vector of nodes to remove from the tree.
#' @param title Character; title of the plot, default is "Feature Tree".
#' @return A ggplot object representing the dendrogram with optional heatmap overlay.
#' @importFrom dplyr filter mutate select
#' @importFrom tidyr pivot_wider
#' @importFrom ggtree ggtree geom_tiplab
#' @export
visualiseModelTree <- function(fit, tree, type = "cluster", heatmap = TRUE,
                               trainingData, nodesToRemove = NULL,
                               title = "Feature Tree") {
  coefs <- coef(fit) %>%
    as.matrix() %>%
    as.data.frame()
  coefs$feature <- rownames(coefs)

  # Plot the base tree structure
  p <- plotOverlapTree(
    modelFit = fit, trainingData = trainingData,
    treeData = tree, nodesToRemove = nodesToRemove
  )

  if (!heatmap) {
    p <- p + labs(title = "Proportions") +
      theme(
        legend.title = element_text(color = "black", size = 12),
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
        legend.text = element_text(color = "black", size = 12)
      )
    return(p)
  }

  # Filter the non-zero coefficients and exclude specific patterns
  coefs <- coefs %>%
    dplyr::filter(!feature %like% "logit") %>%
    dplyr::filter(!feature %like% "(Intercept)") %>%
    dplyr::filter(!feature %like% "Age") %>%
    dplyr::filter(!feature %like% "GenderM")

  # Split feature strings into clusters and markers based on the specified type
  splits <- strsplit(coefs$feature, "_", fixed = TRUE)

  if (type == "cluster") {
    coefs$cellType <- sapply(splits, function(x) paste(x[1], x[2], sep = "_"))
    coefs$marker <- sapply(splits, function(x) x[3])
  } else {
    coefs$cellType <- gsub("_.*", "", coefs$feature)
    coefs$marker <- sub(".*_", "", coefs$feature)
  }

  # Reshape data to wide format for heatmap plotting
  coefsWide <- coefs %>%
    dplyr::select(cellType, marker, V1) %>%
    pivot_wider(names_from = cellType, values_from = V1) %>%
    as.data.frame()

  rownames(coefsWide) <- coefsWide$marker
  coefsWide <- coefsWide %>% dplyr::select(-marker)

  colnames(coefsWide) <- str_replace(colnames(coefsWide), "logit", "")
  colnames(coefsWide) <- str_replace(colnames(coefsWide), "_", " ")

  # Determine color scale limits and colors for the heatmap
  valRange <- range(coefsWide)
  if (valRange[1] < 1e-16) valRange[1] <- -valRange[2]

  colorBreaks <- c(valRange[1], 0, valRange[2])
  colors <- c("#4575B4", "white", "#D73027")

  # Order the rows and apply the heatmap to the plot
  coefsWide <- coefsWide[order(rownames(coefsWide), decreasing = FALSE), ]
  gheatmap(p,
    data = t(coefsWide), offset = 0, hjust = 0.5, colnames_offset_y = -0.4,
    colnames_offset_x = 0, color = "black", font.size = 4, width = 2.5
  ) +
    labs(title = title) +
    theme(
      legend.title = element_text(color = "black", size = 12),
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      legend.text = element_text(color = "black", size = 12),
      legend.key.size = unit(2, "cm"),
      legend.key.height = unit(1, "cm"),
      legend.key.width = unit(1, "cm")
    ) +
    scale_fill_gradientn(
      limits = valRange,
      colours = colors[c(1, seq_along(colors), length(colors))],
      values = c(0, scales::rescale(colorBreaks, from = valRange), 1),
      name = "Beta Means",
      guide = guide_colorbar(order = 1)
    )
}


#' Plot Expression of Significant Features
#'
#' Creates a boxplot to visualize the expression levels of significant features
#' identified by the model, grouped by outcome. Optionally overlays statistical
#' significance labels on the plot.
#' @param stats Data frame containing statistical results for each feature.
#' @param sigFeatures Data frame of significant features and their expression levels.
#' @param outcome Character; name of the outcome variable to group by.
#' @param title Character; title of the plot, default is an empty string.
#' @param test Logical; if TRUE, applies significance testing and displays p-values, default is FALSE.
#' @return A ggplot object representing the boxplot of significant features.
#' @importFrom ggpubr ggboxplot stat_pvalue_manual
#' @importFrom dplyr group_by select mutate
#' @importFrom reshape2 melt
#' @export
plotSigFeatures <- function(stats, sigFeatures, outcome, title = "", type = "density") {
  pltData <- sigFeatures %>%
    as.data.frame() %>%
    dplyr::group_by(!!dplyr::sym(outcome)) %>%
    dplyr::select(stats$feature) %>%
    as.data.frame() %>%
    melt() %>%
    dplyr::mutate(feature = variable)

  if (type == "density") {
    p <- pltData %>%
      ggplot(aes(value, color = !!dplyr::sym(outcome))) +
      geom_density() +
      facet_wrap(~feature, scales = "free") +
      theme_minimal() +
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 16),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = "top",
        strip.text.x = element_text(size = 16, colour = "black", face = "bold")
      ) +
      scale_color_jco()
  } else {
    p <- pltData %>%
      ggboxplot(
        x = outcome, y = "value",
        color = outcome, palette = "jco",
        facet.by = "feature", scale = "free_y"
      ) +
      labs(title = title) +
      xlab("Gensini") +
      ylab("Mean Expression") +
      theme_minimal() +
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 16),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = "None",
        strip.text.x = element_text(size = 16, colour = "black", face = "bold")
      ) +
      stat_pvalue_manual(
        stats,
        label = "p.adj", tip.length = 0,
        size = 6, color = "red"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  }


  return(p)
}


#' Plot AUC for Model Evaluation
#'
#' Plots the ROC curve and calculates AUC for model predictions. Also returns
#' a data frame with prediction results.
#' @param fit Fitted model object.
#' @param xTest Matrix of test data.
#' @param yTest Vector of true labels.
#' @param title Character; plot title.
#' @return A list with a ggplot object for the ROC curve and a data frame of predictions.
#' @importFrom ROCR prediction performance
#' @importFrom ggplot2 ggplot aes geom_line geom_abline labs theme_bw
#' @export
plotAUC <- function(fit, xTest, yTest, title = "") {
  # Get predictions from the model
  preds <- predict(fit, X = xTest, type = "response")
  predsClass <- predict(fit, X = xTest, type = "class")
  predsData <- data.frame(Truth = yTest, Predicted = preds, PredictedClass = predsClass)

  # Calculate ROC curve
  predObj <- ROCR::prediction(preds, yTest)
  perf <- ROCR::performance(predObj, "tpr", "fpr")
  aucValue <- ROCR::performance(predObj, measure = "auc")@y.values[[1]]

  # Prepare data for plotting
  rocData <- data.frame(
    FPR = perf@x.values[[1]],
    TPR = perf@y.values[[1]]
  )

  # Generate ROC curve plot
  aucPlot <- ggplot(rocData, aes(x = FPR, y = TPR)) +
    geom_line(color = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "False Positive Rate", y = "True Positive Rate") +
    ggtitle(paste(title, round(aucValue, 2))) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16)
    )

  return(list(plot = aucPlot, preds = predsData))
}

#' Plot Elbow Graph for BICs
#'
#' Generates an elbow plot of BIC values to aid in selecting the optimal alpha value.
#' @param bics Data frame with BIC values for each alpha.
#' @return A ggplot object showing the elbow plot.
#' @importFrom ggplot2 ggplot aes geom_vline labs theme element_text
#' @importFrom ggpubr ggline
#' @export
plotElbow <- function(bics) {
  # Prepare data for plotting
  plotData <- data.frame(
    Alpha = seq(0.01, 0.1, 0.01),
    BIC = bics$bics[1:10]
  )

  # Generate elbow plot
  elbowPlot <- ggpubr::ggline(
    plotData,
    x = "Alpha", y = "BIC", group = 1, color = "steelblue"
  ) +
    geom_vline(xintercept = which.min(plotData$BIC), linetype = 2, color = "red") +
    labs(
      x = "Alpha", y = "BIC",
      title = "Optimal Alpha using the Elbow Method"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold")
    )

  return(elbowPlot)
}

#' Select Optimal Alpha for Group Lasso Model
#'
#' Searches for the best alpha value for the group lasso model by minimizing BIC.
#' @param xTrain Matrix of training features.
#' @param yTrain Vector of training labels.
#' @param groups List of groups for overlapping group lasso.
#' @param penalty Character; penalty type, default is "cMCP".
#' @param weights Vector of weights for the lasso penalty.
#' @param alphaSearch Numeric vector of alpha values to search.
#' @param seed Integer; random seed for reproducibility.
#' @return A list with the best fitted model, best BIC, best alpha, and BIC values for each alpha.
#' @export
selectAlpha <- function(xTrain, yTrain, groups, penalty = "cMCP", weights,
                        alphaSearch = seq(0.01, 0.1, by = 0.01), seed = 1994) {
  bestBIC <- Inf
  bics <- numeric(length(alphaSearch))

  max_lambda <- max(abs(t(xTrain) %*% yTrain)) / nrow(xTrain) # Based on the largest gradient
  lambda_grid <- seq(0.001 * max_lambda, max_lambda, length.out = 1000) # Minimum value set to 0.05

  for (i in seq_along(alphaSearch)) {
    alpha <- alphaSearch[i]
    message(paste("Testing alpha:", alpha))

    cvFit <- cv.grpregOverlap(xTrain, yTrain, groups,
      nfolds = 3, seed = seed,
      lambda = lambda_grid,
      penalty = penalty, family = "binomial", alpha = alpha
    )

    fit <- grpregOverlap(xTrain, yTrain, groups,
      family = "binomial", alpha = alpha,
      returnX.latent = TRUE, returnOverlap = FALSE, lambda = cvFit$lambda.min,
      penalty = penalty
    )
    currentBIC <- BIC(fit) %>% round(2)
    bics[i] <- currentBIC

    if (currentBIC < bestBIC) {
      bestAlpha <- alpha
      bestBIC <- currentBIC
      bestFit <- fit
      message(paste("Current best BIC:", currentBIC))
    }
  }

  return(list(bestFit = bestFit, bestBIC = bestBIC, bestAlpha = bestAlpha, bics = bics))
}


selectLambdaParallel <- function(xTrain, yTrain, groups, penalty = "cMCP",
                                 alpha = 1, lambdaSearch = seq(0.01, 1, by = 0.01),
                                 seed = 1994, BPPARAM = MulticoreParam(progressbar = TRUE, force.GC = F, workers = 10)) {
  set.seed(seed)
  bestBIC <- Inf
  bestLambda <- NA
  bestFit <- NULL

  # Define function to process each lambda
  processLambda <- function(lambda) {
    fit <- tryCatch(
      {
        grpregOverlap(xTrain, yTrain, groups,
          family = "binomial", alpha = alpha,
          lambda = lambda, penalty = penalty, seed = seed,
          returnX.latent = TRUE, returnOverlap = FALSE
        )
      },
      error = function(e) NULL
    )

    # Compute BIC if model fitting was successful
    if (!is.null(fit)) {
      currentBIC <- tryCatch(
        {
          round(BIC(fit), 2)
        },
        error = function(e) Inf
      )
      return(list(lambda = lambda, BIC = currentBIC, fit = fit))
    } else {
      return(list(lambda = lambda, BIC = Inf, fit = NULL))
    }
  }

  # Run the processLambda function in parallel over lambdaSearch with progress bar
  results <- bplapply(lambdaSearch, processLambda, BPPARAM = BPPARAM)

  # Combine results and find the best lambda
  bics <- sapply(results, function(res) res$BIC)
  bestIndex <- which.min(bics)

  if (length(bestIndex) > 0 && bics[bestIndex] != Inf) {
    bestBIC <- bics[bestIndex]
    bestLambda <- results[[bestIndex]]$lambda
    bestFit <- results[[bestIndex]]$fit
  }

  # Print optimal lambda
  message(paste("Optimal lambda:", bestLambda, "with BIC:", bestBIC))

  # Return results
  return(list(bestFit = bestFit, bestBIC = bestBIC, bestLambda = bestLambda, bics = bics))
}

#' Select Optimal Lambda for Group Lasso Model
#'
#' Searches for the best lambda value for the group lasso model by minimizing BIC.
#' @param xTrain Matrix of training features.
#' @param yTrain Vector of training labels.
#' @param groups List of groups for overlapping group lasso.
#' @param penalty Character; penalty type, default is "cMCP".
#' @param weights Vector of weights for the lasso penalty.
#' @param alpha Numeric; fixed alpha value for the group lasso.
#' @param lambdaSearch Numeric vector of lambda values to search.
#' @param seed Integer; random seed for reproducibility.
#' @return A list with the best fitted model, best BIC, best lambda, and BIC values for each lambda.
#' @export
selectLambda <- function(xTrain, yTrain, groups, penalty = "cMCP",
                         alphaRange = seq(0.1, 1, length = 10),
                         lambdaSearch = seq(0.01, 1, length = 10), seed = 1994) {
    set.seed(seed)
    
    # Initialize variables to store the best results
    bestBic <- Inf
    bestLambda <- NA
    bestAlpha <- NA
    bestFit <- NULL
    allBics <- matrix(NA, nrow = length(alphaRange), ncol = length(lambdaSearch))
    rownames(allBics) <- paste0("Alpha_", alphaRange)
    colnames(allBics) <- paste0("Lambda_", round(lambdaSearch, 4))
    
    # Loop over alpha values
    for (alphaIdx in seq_along(alphaRange)) {
        alpha <- alphaRange[alphaIdx]
        cat("Testing alpha =", alpha, "\n")
        
        # Fit model for the given alpha and lambda range
        fit <- tryCatch(
            {
                grpregOverlap(
                    xTrain, yTrain, groups,
                    family = "binomial", alpha = alpha,
                    lambda = lambdaSearch, penalty = penalty, seed = seed,
                    returnX.latent = TRUE, returnOverlap = FALSE
                )
            },
            error = function(e) {
                message("Error fitting model for alpha =", alpha, ": ", e$message)
                return(NULL)
            }
        )
        
        # If the model fitting fails, skip to the next alpha
        if (is.null(fit)) next
        
        # Extract deviance and calculate BIC for each lambda
        deviances <- fit$deviance
        n <- nrow(xTrain)
        validLambda <- length(deviances) # Actual number of valid lambdas
        bic <- deviances + (fit$df * log(n))
        allBics[alphaIdx, seq_len(validLambda)] <- bic # Store valid BICs
        
        # Find the best lambda for this alpha
        minBicIdx <- computeElbow(bic)
        minBic <- bic[minBicIdx]
        
        if (minBic < bestBic) {
            bestBic <- minBic
            bestLambda <- fit$lambda[minBicIdx] # Use actual lambda values
            bestAlpha <- alpha
            bestFit <- fit
            cat("New best BIC:", bestBic, "at alpha =", bestAlpha, "lambda =", bestLambda, "\n")
        }
    }
    
    # Return the best results and all BICs
    return(list(
        bestFit = bestFit,
        bestBic = bestBic,
        bestLambda = bestLambda,
        bestAlpha = bestAlpha,
        allBics = allBics
    ))
}

#' Fit Overlapping Group Lasso Model
#'
#' Fits the group lasso model to data with the selected alpha. If alpha is not provided,
#' searches for optimal alpha using BIC.
#' @param xTrain Matrix of training features.
#' @param yTrain Vector of training labels.
#' @param groups List of groups for overlapping group lasso.
#' @param alpha Numeric; alpha value for regularization.
#' @param alphaSearch Numeric vector of alpha values to search, if alpha is NULL.
#' @param penalty Character; penalty type, default is "cMCP".
#' @param modelSummary Logical; should the model summary be printed
#' @param seed Integer; random seed for reproducibility.
#' @return List containing the fitted model.
#' @export
fitModel <- function(xTrain, yTrain, groups, alpha = NULL, lambda = NULL,
                     alphaSearch = seq(0.01, 0.1, 0.01),
                     penalty = "cMCP", modelSummary = F,
                     seed = 1994) {
  # Determine optimal alpha if not provided
  message("Computing optimal alpha using the BIC method")
  # Fit the model with selected alpha
  max_lambda <- max(abs(t(xTrain) %*% yTrain)) / nrow(xTrain) # Based on the largest gradient
  lambda_grid <- seq(0.05, 1, length.out = 10000) # Minimum value set to 0.05

  if (is.null(alpha) && is.null(lambda)) {
    resBIC <- selectLambda(xTrain, yTrain,
      groups = groups, penalty = penalty,
      lambdaSearch = lambda_grid, seed = seed
    )
    # print(length(resBIC$bics))
    # lambda <- lambda_grid[computeElbow(resBIC$bics)]
    lambda <- resBIC$bestLambda
    alpha <- resBIC$bestAlpha
  }

  message(paste("Optimal alpha value:", alpha))
  message(paste("Optimal lambda value:", lambda))
  message("Fitting the final model")

  # # Fit the model with selected alpha
  # cvFit <- cv.grpregOverlap(xTrain, yTrain, groups, penalty = penalty,
  #                           family = 'binomial', alpha = alpha, nfolds = 3,
  #                           seed = seed, lamba = lambda_grid)
  # return(cvFit)
  if (modelSummary) {
    par(mfrow = c(2, 2))
    grpreg:::plot.cv.grpreg(cvFit, type = "cve")
    grpreg:::plot.cv.grpreg(cvFit, type = "pred")
    grpreg:::plot.cv.grpreg(cvFit, type = "snr")
    grpreg:::plot.cv.grpreg(cvFit, type = "rsq")
    print(summary(cvFit))
  }

  # print(cvFit$lambda.min)

  fit <- grpregOverlap(xTrain, yTrain, groups,
    family = "binomial", alpha = alpha,
    returnX.latent = TRUE, returnOverlap = FALSE, lambda = lambda,
    penalty = penalty
  )
  message("Model fitting complete")

  return(list(fit = fit))
}

#' Plot Heatmap of Marker Means Coefficients from Group Lasso Model
#'
#' Creates a heatmap of coefficients from the group lasso model, grouped by cell type clusters.
#' @param fit Fitted group lasso model.
#' @param type Character; type of clustering, default is "cluster".
#' @param order Optional vector to order clusters in the heatmap.
#' @param markerOrder Optional vector to order markers in the heatmap.
#' @return ggplot object representing the heatmap.
#' @importFrom dplyr filter select
#' @importFrom tidyr pivot_wider
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn theme element_text
#' @export
plotHeatmap <- function(fit, type = "cluster", order = NULL, markerOrder = NULL) {
  # Extract and filter coefficients
  coefs <- coef(fit) %>%
    as.matrix() %>%
    as.data.frame() %>%
    dplyr::mutate(feature = rownames(.)) %>%
    dplyr::filter(!grepl("logit|Intercept|Age|GenderM", feature))

  # Split feature into Cell_Type and Marker
  splits <- strsplit(coefs$feature, "_", fixed = TRUE)
  if (type == "cluster") {
    coefs$Cell_Type <- sapply(splits, function(x) paste(x[1], x[2], sep = "_"))
    coefs$Marker <- sapply(splits, function(x) x[3])
  } else {
    coefs$Cell_Type <- gsub("_.*", "", coefs$feature)
    coefs$Marker <- sub(".*_", "", coefs$feature)
  }

  # Reshape data for heatmap
  coefsWide <- coefs %>%
    dplyr::select(Cell_Type, Marker, V1) %>%
    tidyr::pivot_wider(names_from = Cell_Type, values_from = V1) %>%
    as.data.frame()
  rownames(coefsWide) <- coefsWide$Marker
  coefsWide <- coefsWide %>% dplyr::select(-Marker)

  colnames(coefsWide) <- gsub("logit", "", colnames(coefsWide))
  colnames(coefsWide) <- gsub("_", " ", colnames(coefsWide))

  # Set color range and breaks
  valRange <- range(coefsWide, na.rm = TRUE)
  colourBreaks <- c(valRange[1], 0, valRange[2])
  if (min(valRange) < 0) {
    colours <- c("#4575B4", "white", "#D73027")
  } else {
    colours <- c("white", "#4575B4", "#D73027")
  }

  # Reorder columns if specified
  if (!is.null(order)) {
    coefsWide <- coefsWide[, order, drop = FALSE]
  }

  # Melt data for ggplot
  coefsWide$Gene <- rownames(coefsWide)
  coefsMelted <- reshape2::melt(coefsWide)
  coefsMelted$Gene <- ordered(coefsMelted$Gene, levels = rev(sort(unique(coefsMelted$Gene))))
  coefsMelted <- coefsMelted[order(coefsMelted$Gene, decreasing = T), ]
  colnames(coefsMelted)[2:3] <- c("Cluster", "Value")

  # Plot heatmap
  heatmapPlot <- ggplot(coefsMelted, aes(x = Cluster, y = Gene, fill = Value)) +
    geom_tile(color = "black", lwd = 0.05, linetype = 1) +
    ggtitle("Marker Means") +
    theme(
      legend.title = element_text(color = "black", size = 12),
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      legend.text = element_text(color = "black", size = 12),
      axis.text.x = element_text(
        angle = 90, vjust = 0.5,
        hjust = 1, size = 12, colour = "black"
      ),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.key.size = unit(2, "cm"),
      legend.key.height = unit(1, "cm"),
      legend.key.width = unit(1, "cm")
    ) +
    scale_fill_gradientn(
      limits = valRange,
      colours = colours,
      values = scales::rescale(colourBreaks, from = valRange),
      name = "Beta Means"
    )

  return(heatmapPlot)
}

#' Train and Predict Cell Type Classifier Using caret
#'
#' Trains a cell type classifier using LDA (or other specified models) and performs
#' predictions on test data.
#' @param trainData Data frame; training data with cell type labels.
#' @param testData Data frame; test data to predict cell types.
#' @param model Character; machine learning model to use, default is "lda".
#' @return Vector of predicted cell types for the test data.
#' @importFrom caret train trainControl
#' @importFrom doParallel registerDoParallel
#' @export
cellTypeClassifier <- function(trainData, testData, model = "lda") {
  # Set up cross-validation control
  fitControl <- trainControl(method = "cv", number = 3)

  # Fit the model
  message("Fitting cell type model")
  set.seed(1994)
  fit <- caret::train(
    cellTypes ~ .,
    data = trainData, method = model,
    trControl = fitControl, preprocess = c("range")
  )
  print(fit)

  # Predict on test data
  message("Predicting on test data")
  testCellTypes <- predict(fit, testData)

  return(testCellTypes)
}

#' Compute Elbow Point for BIC Plot
#'
#' Identifies the optimal elbow point in a BIC plot, commonly used for selecting
#' the alpha parameter in group lasso models.
#' @param vals Numeric vector of BIC values across different alpha values.
#' @return Integer; index of the elbow point, indicating the optimal alpha.
#' @export
computeElbow <- function(vals) {
  # Calculate differences between consecutive BIC values
  diffs <- diff(vals)
  # diffs <- diffs[-1]
  # Identify the index of the largest change, interpreted as the elbow
  optKb <- which.max(abs(diffs)) + 1

  return(optKb)
}
