#' Generate an Extended Q-matrix with Rule Interactions and Collinearity Diagnostics
#'
#' This function generates interaction terms between rules within the same component,
#' extends the Q-matrix, and evaluates the resulting matrix for collinearity issues
#' using eigenvalues, condition indices, and variance inflation factors (VIF).
#' If severe collinearity is detected, it attempts to iteratively remove problematic
#' interaction terms while keeping the original rules untouched.
#'
#' @param Q A binary matrix of items by rules (original Q-matrix). Each row represents
#'   an item and each column represents a rule. Values should be 0 or 1.
#' @param M_list A list where each element contains the indices of rules that belong
#'   to the same component/dimension. For example, list(c(1,2,3), c(4,5)) indicates
#'   that rules 1,2,3 belong to component 1 and rules 4,5 belong to component 2.
#' @param max_condition_index Numeric. Maximum acceptable condition index. Default is 30.
#'   Values above this threshold indicate severe collinearity.
#' @param min_eigenvalue Numeric. Minimum acceptable eigenvalue. Default is 0.1.
#'   Values below this threshold may indicate linear dependence.
#' @param plot_diagnostics Logical. Whether to generate diagnostic plots. Default is TRUE.
#' @param verbose Logical. Whether to print detailed diagnostic information. Default is TRUE.
#' @param save_to_global Logical. Whether to save results to global environment. Default is TRUE.
#'
#' @return A list containing:
#'   \item{Q_extended}{The extended Q-matrix with interaction terms}
#'   \item{M_list_extended}{Updated component list including interaction terms}
#'   \item{diagnostics}{List with collinearity diagnostics}
#'   \item{removed_interactions}{Vector of removed interaction names (if any)}
#'   \item{plots}{List of diagnostic plots (if plot_diagnostics = TRUE)}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input parameters
#'   \item Generates interaction terms for rules within the same component
#'   \item Performs collinearity diagnostics using multiple methods
#'   \item Attempts to resolve severe collinearity by removing problematic interactions
#'   \item Generates diagnostic plots and summaries
#' }
#'
#' Collinearity is assessed using:
#' \itemize{
#'   \item Condition indices (based on eigenvalues of correlation matrix)
#'   \item Variance Inflation Factors (VIF)
#'   \item Matrix rank assessment
#'   \item Eigenvalue analysis
#' }
#'
#' @examples
#' \dontrun{
#' # Create a sample Q-matrix (5 items, 4 rules)
#' Q <- matrix(c(1,1,0,0,0,
#'               1,0,1,0,0,
#'               0,1,1,0,0,
#'               0,0,0,1,1), nrow=5, ncol=4, byrow=FALSE)
#'
#' # Define components (rules 1-2 in component 1, rules 3-4 in component 2)
#' M_list <- list(c(1,2), c(3,4))
#'
#' # Generate extended Q-matrix with interactions
#' result <- generate_Q_with_interactions(Q, M_list)
#'
#' # Access results
#' extended_Q <- result$Q_extended
#' diagnostics <- result$diagnostics
#' plots <- result$plots
#' }
#'
#' @references
#' Belsley, D. A., Kuh, E., & Welsch, R. E. (1980). Regression diagnostics:
#'   Identifying influential data and sources of collinearity. John Wiley & Sons.
#' O'Brien, R. M. (2007). A caution regarding rules of thumb for variance
#'   inflation factors. Quality & Quantity, 41(5), 673-690.
#' Hair, J. F., Anderson, R. E., Tatham, R. L., & Black, W. C. (1995).
#'   Multivariate data analysis. Prentice Hall.
#'
#' @importFrom graphics par plot hist barplot
#' @importFrom grDevices dev.new
#' @export
generate_Q_with_interactions <- function(Q,
                                         M_list,
                                         max_condition_index = 30,
                                         min_eigenvalue = 0.1,
                                         plot_diagnostics = TRUE,
                                         verbose = TRUE,
                                         save_to_global = TRUE) {

  # Package-level environment for results (avoids global env modification)
  .GMLTM_env <- new.env(parent = emptyenv())

  # ============================================================================
  # 1. INPUT VALIDATION
  # ============================================================================

  if (!is.matrix(Q)) {
    stop("Q must be a matrix.")
  }

  if (!is.numeric(Q) || !all(Q %in% c(0, 1))) {
    stop("Q must be a binary matrix containing only 0s and 1s.")
  }

  if (!is.list(M_list) || length(M_list) == 0) {
    stop("M_list must be a non-empty list.")
  }

  if (any(unlist(M_list) > ncol(Q)) || any(unlist(M_list) < 1)) {
    stop("All indices in M_list must be between 1 and ncol(Q).")
  }

  if (max_condition_index <= 0 || min_eigenvalue <= 0) {
    stop("max_condition_index and min_eigenvalue must be positive.")
  }

  # ============================================================================
  # 2. INITIALIZE VARIABLES
  # ============================================================================

  K <- ncol(Q)  # Number of original rules
  I <- nrow(Q)  # Number of items

  # Assign column names if not present
  if (is.null(colnames(Q))) {
    colnames(Q) <- paste0("Rule_", 1:K)
  }

  # Initialize extended matrices and lists
  Q_extended <- Q
  new_rule_names <- colnames(Q)
  M_list_extended <- M_list
  rule_counter <- K
  interaction_log <- data.frame(
    component = integer(),
    rule_i = integer(),
    rule_j = integer(),
    interaction_name = character(),
    non_zero_items = integer(),
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat("=== Q-MATRIX INTERACTION GENERATION ===\n")
    cat("Original Q-matrix dimensions:", I, "items x", K, "rules\n")
    cat("Number of components:", length(M_list), "\n\n")
  }

  # ============================================================================
  # 3. GENERATE INTERACTION TERMS
  # ============================================================================

  interactions_generated <- 0

  for (m in seq_along(M_list)) {
    reglas_comp <- M_list[[m]]

    if (length(reglas_comp) < 2) {
      if (verbose) {
        cat("Component", m, ": Only", length(reglas_comp),
            "rule(s). No interactions generated.\n")
      }
      next
    }

    component_interactions <- 0

    for (i in 1:(length(reglas_comp) - 1)) {
      for (j in (i + 1):length(reglas_comp)) {
        regla_i <- reglas_comp[i]
        regla_j <- reglas_comp[j]

        # Create interaction term
        nueva_inter <- Q[, regla_i] * Q[, regla_j]
        non_zero_count <- sum(nueva_inter)

        # Only add interaction if it has non-zero entries
        if (non_zero_count > 0) {
          Q_extended <- cbind(Q_extended, nueva_inter)
          rule_counter <- rule_counter + 1

          nombre_inter <- paste0("Inter_", regla_i, "*", regla_j)
          new_rule_names <- c(new_rule_names, nombre_inter)
          M_list_extended[[m]] <- c(M_list_extended[[m]], rule_counter)

          # Log interaction
          interaction_log <- rbind(interaction_log, data.frame(
            component = m,
            rule_i = regla_i,
            rule_j = regla_j,
            interaction_name = nombre_inter,
            non_zero_items = non_zero_count,
            stringsAsFactors = FALSE
          ))

          interactions_generated <- interactions_generated + 1
          component_interactions <- component_interactions + 1
        }
      }
    }

    if (verbose) {
      cat("Component", m, ": Generated", component_interactions, "interactions\n")
    }
  }

  colnames(Q_extended) <- new_rule_names

  if (verbose) {
    cat("\nTotal interactions generated:", interactions_generated, "\n")
    cat("Extended Q-matrix dimensions:", I, "items x", ncol(Q_extended), "rules\n\n")
  }

  # ============================================================================
  # 4. COLLINEARITY DIAGNOSTICS FUNCTIONS
  # ============================================================================

  calculate_vif <- function(X) {
    if (ncol(X) < 2) return(rep(1, ncol(X)))

    vif_values <- numeric(ncol(X))
    for (i in 1:ncol(X)) {
      y <- X[, i]
      x <- X[, -i, drop = FALSE]

      # Handle perfect collinearity
      if (qr(x)$rank < ncol(x)) {
        vif_values[i] <- Inf
      } else {
        r_squared <- cor(y, x %*% solve(t(x) %*% x) %*% t(x) %*% y)^2
        vif_values[i] <- 1 / (1 - r_squared)
      }
    }
    names(vif_values) <- colnames(X)
    return(vif_values)
  }

  perform_collinearity_diagnostics <- function(X, label = "") {
    if (verbose && label != "") {
      cat("=== COLLINEARITY DIAGNOSTICS:", label, "===\n")
    }

    # Basic matrix properties
    rank_X <- qr(X)$rank
    full_rank <- rank_X == ncol(X)

    if (verbose) {
      cat("Matrix rank:", rank_X, "/", ncol(X),
          ifelse(full_rank, "(Full rank)", "(Rank deficient)"), "\n")
    }

    # Correlation matrix and eigenanalysis
    cor_X <- cor(X)
    eig_vals <- eigen(cor_X)$values
    cond_indices <- sqrt(max(eig_vals) / eig_vals)

    # VIF calculation (only if matrix has full rank)
    vif_values <- NULL
    if (full_rank && ncol(X) > 1) {
      tryCatch({
        vif_values <- calculate_vif(X)
      }, error = function(e) {
        if (verbose) cat("Warning: Could not calculate VIF values\n")
        vif_values <<- rep(NA, ncol(X))
      })
    }

    # Collinearity assessment
    severe_collinearity <- any(cond_indices > max_condition_index, na.rm = TRUE)
    moderate_collinearity <- any(cond_indices > 10 & cond_indices <= max_condition_index, na.rm = TRUE)
    low_eigenvalues <- sum(eig_vals < min_eigenvalue)
    high_vif <- if (!is.null(vif_values)) sum(vif_values > 10, na.rm = TRUE) else 0

    if (verbose) {
      cat("Condition indices range: [", round(min(cond_indices), 2), ", ",
          round(max(cond_indices), 2), "]\n")
      cat("Eigenvalues below", min_eigenvalue, ":", low_eigenvalues, "\n")
      if (!is.null(vif_values)) {
        cat("VIF values above 10:", high_vif, "\n")
      }

      # Warnings
      if (!full_rank) {
        cat("[WARNING] WARNING: Matrix is rank deficient!\n")
      }
      if (severe_collinearity) {
        cat("[WARNING] WARNING: Severe collinearity detected (CI >", max_condition_index, ")!\n")
      } else if (moderate_collinearity) {
        cat("[WARNING] WARNING: Moderate collinearity detected (CI > 10)!\n")
      }
      if (low_eigenvalues > 0) {
        cat("[WARNING] WARNING:", low_eigenvalues, "eigenvalue(s) below", min_eigenvalue, "\n")
      }
      cat("\n")
    }

    return(list(
      rank = rank_X,
      full_rank = full_rank,
      correlation_matrix = cor_X,
      eigenvalues = eig_vals,
      condition_indices = cond_indices,
      vif_values = vif_values,
      severe_collinearity = severe_collinearity,
      moderate_collinearity = moderate_collinearity,
      low_eigenvalues = low_eigenvalues,
      high_vif = high_vif
    ))
  }

  # ============================================================================
  # 5. DIAGNOSE ORIGINAL EXTENDED MATRIX
  # ============================================================================

  # Save original versions
  Q_extended_original <- Q_extended
  M_list_original <- M_list_extended

  if (save_to_global) {
    .GMLTM_env$Q_extended_original <- Q_extended_original
    .GMLTM_env$M_list_original <- M_list_original
  }

  # Perform diagnostics on original extended matrix
  diagnostics_original <- perform_collinearity_diagnostics(Q_extended, "ORIGINAL EXTENDED MATRIX")

  # ============================================================================
  # 6. ATTEMPT TO RESOLVE COLLINEARITY ISSUES
  # ============================================================================

  Q_improved <- Q_extended
  M_list_improved <- M_list_extended
  removed_interactions <- character(0)
  diagnostics_improved <- diagnostics_original

  if (diagnostics_original$severe_collinearity || !diagnostics_original$full_rank) {
    if (verbose) {
      cat("=== ATTEMPTING TO RESOLVE COLLINEARITY ===\n")
    }

    # Identify interaction columns (exclude original rules)
    interaction_indices <- (K + 1):ncol(Q_extended)
    interaction_indices <- interaction_indices[interaction_indices <= ncol(Q_extended)]

    if (length(interaction_indices) > 0) {
      to_remove <- c()
      max_iterations <- length(interaction_indices)
      iteration <- 0

      repeat {
        iteration <- iteration + 1
        if (iteration > max_iterations) break

        # Current matrix excluding already marked for removal
        current_indices <- setdiff(1:ncol(Q_extended), to_remove)
        Q_current <- Q_extended[, current_indices, drop = FALSE]

        # Check current diagnostics
        temp_diagnostics <- perform_collinearity_diagnostics(Q_current, "")

        # Stop if acceptable
        if (!temp_diagnostics$severe_collinearity && temp_diagnostics$full_rank) {
          break
        }

        # Find most problematic interaction to remove
        remaining_interactions <- intersect(interaction_indices, current_indices)
        if (length(remaining_interactions) == 0) break

        # Strategy: remove interaction with highest condition index contribution
        worst_col <- remaining_interactions[1]  # Simple strategy
        to_remove <- c(to_remove, worst_col)

        if (verbose) {
          cat("Iteration", iteration, ": Removing", colnames(Q_extended)[worst_col], "\n")
        }
      }

      # Apply removals if successful
      if (length(to_remove) > 0 && length(to_remove) < length(interaction_indices)) {
        Q_improved <- Q_extended[, -to_remove, drop = FALSE]
        removed_interactions <- colnames(Q_extended)[to_remove]

        # Update M_list
        M_list_improved <- lapply(M_list_extended, function(m) {
          m[!m %in% to_remove]
        })

        # Re-diagnose improved matrix
        diagnostics_improved <- perform_collinearity_diagnostics(Q_improved, "IMPROVED MATRIX")

        if (verbose) {
          cat("Successfully removed", length(removed_interactions), "interactions:\n")
          cat(paste(removed_interactions, collapse = ", "), "\n\n")
        }
      } else {
        if (verbose) {
          cat("[ERROR] Could not resolve collinearity issues without removing all interactions.\n\n")
        }
        if (save_to_global) {
          .GMLTM_env$Q_extended_failed <- Q_extended
          .GMLTM_env$M_list_failed <- M_list_extended
        }
      }
    }
  }

  # ============================================================================
  # 7. SAVE FINAL RESULTS
  # ============================================================================

  if (save_to_global) {
    .GMLTM_env$Q_extended <- Q_improved
    .GMLTM_env$M_list <- M_list_improved

    if (verbose) {
      cat("Results saved to .GMLTM_env (accessible via the returned list):\n")
      cat("- Q_extended: Final extended Q-matrix (", nrow(Q_improved), "x", ncol(Q_improved), ")\n")
      cat("- M_list: Updated component list\n")
      cat("- Q_extended_original: Original extended matrix\n")
      cat("- M_list_original: Original component list\n")
      if (!is.null(.GMLTM_env$Q_extended_failed)) {
        cat("- Q_extended_failed: Failed correction attempt\n")
      }
      cat("\n")
    }
  }

  # ============================================================================
  # 8. GENERATE DIAGNOSTIC PLOTS
  # ============================================================================

  plots <- list()

  if (plot_diagnostics) {
    if (verbose) cat("=== GENERATING DIAGNOSTIC PLOTS ===\n")

    # Generate plots immediately and store them as functions
    tryCatch({

      # Plot 1: Condition Indices Comparison
      if (verbose) cat("Generating condition indices plot...\n")

      if (length(removed_interactions) > 0) {
        # Show comparison plots in RStudio plots panel
        par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

        # Original
        barplot(diagnostics_original$condition_indices,
                main = "Condition Indices - Original",
                ylab = "Condition Index",
                las = 2, cex.names = 0.7, col = "lightblue")
        abline(h = max_condition_index, col = "red", lty = 2, lwd = 2)
        abline(h = 10, col = "orange", lty = 2, lwd = 2)
        legend("topright", legend = c("Severe (>30)", "Moderate (>10)"),
               col = c("red", "orange"), lty = 2, cex = 0.8)

        # Improved
        barplot(diagnostics_improved$condition_indices,
                main = "Condition Indices - Improved",
                ylab = "Condition Index",
                las = 2, cex.names = 0.7, col = "lightgreen")
        abline(h = max_condition_index, col = "red", lty = 2, lwd = 2)
        abline(h = 10, col = "orange", lty = 2, lwd = 2)
        legend("topright", legend = c("Severe (>30)", "Moderate (>10)"),
               col = c("red", "orange"), lty = 2, cex = 0.8)

      } else {
        # Show single plot in RStudio plots panel
        par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

        barplot(diagnostics_original$condition_indices,
                main = "Condition Indices",
                ylab = "Condition Index",
                las = 2, cex.names = 0.7, col = "lightblue")
        abline(h = max_condition_index, col = "red", lty = 2, lwd = 2)
        abline(h = 10, col = "orange", lty = 2, lwd = 2)
        legend("topright", legend = c("Severe (>30)", "Moderate (>10)"),
               col = c("red", "orange"), lty = 2, cex = 0.8)
      }

      # Store as function for later use
      plots$condition_indices <- function() {
        if (length(removed_interactions) > 0) {
          par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

          barplot(diagnostics_original$condition_indices,
                  main = "Condition Indices - Original",
                  ylab = "Condition Index",
                  las = 2, cex.names = 0.7, col = "lightblue")
          abline(h = max_condition_index, col = "red", lty = 2, lwd = 2)
          abline(h = 10, col = "orange", lty = 2, lwd = 2)

          barplot(diagnostics_improved$condition_indices,
                  main = "Condition Indices - Improved",
                  ylab = "Condition Index",
                  las = 2, cex.names = 0.7, col = "lightgreen")
          abline(h = max_condition_index, col = "red", lty = 2, lwd = 2)
          abline(h = 10, col = "orange", lty = 2, lwd = 2)
        } else {
          barplot(diagnostics_original$condition_indices,
                  main = "Condition Indices",
                  ylab = "Condition Index",
                  las = 2, cex.names = 0.7, col = "lightblue")
          abline(h = max_condition_index, col = "red", lty = 2, lwd = 2)
          abline(h = 10, col = "orange", lty = 2, lwd = 2)
        }
      }

      # Plot 2: Eigenvalue Distribution
      if (verbose) cat("Generating eigenvalues plot...\n")

      # Wait for user to see previous plot (if interactive)
      if (interactive() && length(removed_interactions) > 0) {
        readline(prompt = "Press Enter to see eigenvalues plot...")
      }

      if (length(removed_interactions) > 0) {
        par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

        plot(diagnostics_original$eigenvalues, type = "b",
             main = "Eigenvalues - Original",
             ylab = "Eigenvalue", xlab = "Component",
             col = "blue", pch = 19, lwd = 2)
        abline(h = min_eigenvalue, col = "red", lty = 2, lwd = 2)
        legend("topright", legend = paste("Threshold =", min_eigenvalue),
               col = "red", lty = 2, cex = 0.8)

        plot(diagnostics_improved$eigenvalues, type = "b",
             main = "Eigenvalues - Improved",
             ylab = "Eigenvalue", xlab = "Component",
             col = "darkgreen", pch = 19, lwd = 2)
        abline(h = min_eigenvalue, col = "red", lty = 2, lwd = 2)
        legend("topright", legend = paste("Threshold =", min_eigenvalue),
               col = "red", lty = 2, cex = 0.8)
      } else {
        par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

        plot(diagnostics_original$eigenvalues, type = "b",
             main = "Eigenvalues",
             ylab = "Eigenvalue", xlab = "Component",
             col = "blue", pch = 19, lwd = 2)
        abline(h = min_eigenvalue, col = "red", lty = 2, lwd = 2)
        legend("topright", legend = paste("Threshold =", min_eigenvalue),
               col = "red", lty = 2, cex = 0.8)
      }

      plots$eigenvalues <- function() {
        if (length(removed_interactions) > 0) {
          par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

          plot(diagnostics_original$eigenvalues, type = "b",
               main = "Eigenvalues - Original",
               ylab = "Eigenvalue", xlab = "Component",
               col = "blue", pch = 19, lwd = 2)
          abline(h = min_eigenvalue, col = "red", lty = 2, lwd = 2)

          plot(diagnostics_improved$eigenvalues, type = "b",
               main = "Eigenvalues - Improved",
               ylab = "Eigenvalue", xlab = "Component",
               col = "darkgreen", pch = 19, lwd = 2)
          abline(h = min_eigenvalue, col = "red", lty = 2, lwd = 2)
        } else {
          plot(diagnostics_original$eigenvalues, type = "b",
               main = "Eigenvalues",
               ylab = "Eigenvalue", xlab = "Component",
               col = "blue", pch = 19, lwd = 2)
          abline(h = min_eigenvalue, col = "red", lty = 2, lwd = 2)
        }
      }

      # Plot 3: Correlation Heatmap
      if (verbose) cat("Generating correlation heatmap...\n")

      # Wait for user to see previous plot (if interactive)
      if (interactive()) {
        readline(prompt = "Press Enter to see correlation heatmap...")
      }

      par(mfrow = c(1, 1), mar = c(5, 4, 4, 6))

      cor_matrix <- diagnostics_improved$correlation_matrix

      # Create color palette
      colors <- colorRampPalette(c("blue", "white", "red"))(100)

      # Create the heatmap
      image(1:ncol(cor_matrix), 1:nrow(cor_matrix),
            as.matrix(cor_matrix),
            col = colors,
            main = "Correlation Matrix Heatmap",
            xlab = "Rules", ylab = "Rules",
            axes = FALSE)

      # Add contour lines for high correlations
      contour(1:ncol(cor_matrix), 1:nrow(cor_matrix),
              as.matrix(cor_matrix),
              levels = c(-0.8, -0.5, 0.5, 0.8),
              add = TRUE, col = "black", lwd = 1)

      # Add axes
      axis(1, at = 1:ncol(cor_matrix), labels = colnames(cor_matrix),
           cex.axis = 0.7, las = 2)
      axis(2, at = 1:nrow(cor_matrix), labels = rownames(cor_matrix),
           cex.axis = 0.7, las = 2)

      # Add color scale legend
      legend("right", legend = c("1.0", "0.5", "0.0", "-0.5", "-1.0"),
             fill = colors[c(100, 75, 50, 25, 1)],
             title = "Correlation", xpd = TRUE, inset = c(-0.1, 0))

      plots$correlation_heatmap <- function() {
        cor_matrix <- diagnostics_improved$correlation_matrix
        colors <- colorRampPalette(c("blue", "white", "red"))(100)

        image(1:ncol(cor_matrix), 1:nrow(cor_matrix),
              as.matrix(cor_matrix),
              col = colors,
              main = "Correlation Matrix Heatmap",
              xlab = "Rules", ylab = "Rules",
              axes = FALSE)

        contour(1:ncol(cor_matrix), 1:nrow(cor_matrix),
                as.matrix(cor_matrix),
                levels = c(-0.8, -0.5, 0.5, 0.8),
                add = TRUE, col = "black", lwd = 1)

        axis(1, at = 1:ncol(cor_matrix), labels = colnames(cor_matrix),
             cex.axis = 0.7, las = 2)
        axis(2, at = 1:nrow(cor_matrix), labels = rownames(cor_matrix),
             cex.axis = 0.7, las = 2)
      }

      # Plot 4: VIF values (if available)
      if (!is.null(diagnostics_improved$vif_values)) {
        vif_vals <- diagnostics_improved$vif_values[is.finite(diagnostics_improved$vif_values)]
        if (length(vif_vals) > 0) {
          if (verbose) cat("Generating VIF plot...\n")

          # Wait for user to see previous plot (if interactive)
          if (interactive()) {
            readline(prompt = "Press Enter to see VIF plot...")
          }

          par(mfrow = c(1, 1), mar = c(8, 4, 4, 2))

          barplot(vif_vals, main = "Variance Inflation Factors",
                  ylab = "VIF", las = 2, cex.names = 0.7,
                  col = ifelse(vif_vals > 10, "red",
                               ifelse(vif_vals > 5, "orange", "lightblue")))
          abline(h = 10, col = "red", lty = 2, lwd = 2)
          abline(h = 5, col = "orange", lty = 2, lwd = 2)
          legend("topright", legend = c("High (>10)", "Moderate (>5)", "Acceptable (<=5)"),
                 fill = c("red", "orange", "lightblue"), cex = 0.8)

          plots$vif_values <- function() {
            vif_vals <- diagnostics_improved$vif_values[is.finite(diagnostics_improved$vif_values)]
            if (length(vif_vals) > 0) {
              par(mfrow = c(1, 1), mar = c(8, 4, 4, 2))
              barplot(vif_vals, main = "Variance Inflation Factors",
                      ylab = "VIF", las = 2, cex.names = 0.7,
                      col = ifelse(vif_vals > 10, "red",
                                   ifelse(vif_vals > 5, "orange", "lightblue")))
              abline(h = 10, col = "red", lty = 2, lwd = 2)
              abline(h = 5, col = "orange", lty = 2, lwd = 2)
            }
          }
        }
      }

      if (verbose) {
        cat("Diagnostic plots generated successfully!\n")
        cat("- Use result$plots$condition_indices() to recreate condition indices plot\n")
        cat("- Use result$plots$eigenvalues() to recreate eigenvalues plot\n")
        cat("- Use result$plots$correlation_heatmap() to recreate correlation heatmap\n")
        if (!is.null(plots$vif_values)) {
          cat("- Use result$plots$vif_values() to recreate VIF plot\n")
        }
        cat("\n")
      }

    }, error = function(e) {
      if (verbose) {
        cat("Warning: Could not generate some diagnostic plots:", e$message, "\n")
      }

      # Provide basic plotting functions even if immediate plotting fails
      plots$condition_indices <- function() {
        if (length(removed_interactions) > 0) {
          par(mfrow = c(1, 2))
          barplot(diagnostics_original$condition_indices, main = "Condition Indices - Original")
          barplot(diagnostics_improved$condition_indices, main = "Condition Indices - Improved")
        } else {
          barplot(diagnostics_original$condition_indices, main = "Condition Indices")
        }
      }

      plots$eigenvalues <- function() {
        if (length(removed_interactions) > 0) {
          par(mfrow = c(1, 2))
          plot(diagnostics_original$eigenvalues, type = "b", main = "Eigenvalues - Original")
          plot(diagnostics_improved$eigenvalues, type = "b", main = "Eigenvalues - Improved")
        } else {
          plot(diagnostics_original$eigenvalues, type = "b", main = "Eigenvalues")
        }
      }
    })
  }

  # ============================================================================
  # 9. PREPARE RETURN OBJECT
  # ============================================================================

  result <- list(
    Q_extended = Q_improved,
    M_list_extended = M_list_improved,
    Q_original = Q,
    M_list_original_input = M_list,
    diagnostics = list(
      original = diagnostics_original,
      improved = diagnostics_improved
    ),
    removed_interactions = removed_interactions,
    interaction_log = interaction_log,
    plots = plots,
    summary = list(
      original_dimensions = c(I, K),
      extended_dimensions = c(nrow(Q_improved), ncol(Q_improved)),
      interactions_generated = interactions_generated,
      interactions_removed = length(removed_interactions),
      final_collinearity_status = ifelse(
        diagnostics_improved$severe_collinearity, "Severe",
        ifelse(diagnostics_improved$moderate_collinearity, "Moderate", "Acceptable")
      )
    )
  )

  if (verbose) {
    cat("=== SUMMARY ===\n")
    cat("Original Q-matrix:", result$summary$original_dimensions[1], "x",
        result$summary$original_dimensions[2], "\n")
    cat("Extended Q-matrix:", result$summary$extended_dimensions[1], "x",
        result$summary$extended_dimensions[2], "\n")
    cat("Interactions generated:", result$summary$interactions_generated, "\n")
    cat("Interactions removed:", result$summary$interactions_removed, "\n")
    cat("Final collinearity status:", result$summary$final_collinearity_status, "\n")
    cat("=======================================\n")
  }

  return(result)
}
