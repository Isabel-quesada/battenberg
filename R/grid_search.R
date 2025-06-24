runASCAT_enhanced = function(lrr, baf, lrrsegmented, bafsegmented, chromosomes, dist_choice,
                               distancepng = NA, copynumberprofilespng = NA, nonroundedprofilepng = NA,
                               cnaStatusFile = "copynumber_solution_status.txt", gamma = 0.55, allow100percent,
                               reliabilityFile=NA, min.ploidy=1.6, max.ploidy=4.8, min.rho=0.1, max.rho=1.0,
                               min.goodness=63, uninformative_BAF_threshold = 0.51, chr.names, analysis="paired",
                               enable_robust_fallback = TRUE, verbose = TRUE, parallel = TRUE,
                               n_cores = NULL) {

  start_time <- Sys.time()

  # Setup data processing (identical to original)
  ch = chromosomes
  b = bafsegmented
  r = lrrsegmented[names(bafsegmented)]

  dist_min_psi = max(min.ploidy-0.6, 0)
  dist_max_psi = max.ploidy+0.6
  dist_min_rho = max(min.rho-0.03, 0.05)
  dist_max_rho = max.rho+0.03

  s = ASCAT::make_segments(r,b)
  dist_matrix_info <- create_distance_matrix(s, dist_choice, gamma,
                                           uninformative_BAF_threshold=uninformative_BAF_threshold,
                                           min_psi=dist_min_psi, max_psi=dist_max_psi,
                                           min_rho=dist_min_rho, max_rho=dist_max_rho)
  d = dist_matrix_info$distance_matrix
  minimise = dist_matrix_info$minimise

  TheoretMaxdist = sum(rep(0.25,dim(s)[1]) * s[,"length"],na.rm=T)

  if( !(minimise) ) {
    d = -d
  }

  if (verbose) cat("Phase 1: Original grid search optimization...\n")

  # PHASE 1: Run original algorithm exactly (fastest path)
  original_results <- run_original_battenberg_search(d, s, gamma, min.ploidy, max.ploidy,
                                                     min.rho, max.rho, min.goodness,
                                                     TheoretMaxdist, minimise, allow100percent, verbose)

  # Check if original method succeeded with good quality
  if (length(original_results) > 0) {
    best_original <- original_results[[which.max(sapply(original_results, function(x) x$goodness))]]

    # Accept original result if it's high quality OR if robust fallback is disabled
    if (!enable_robust_fallback || best_original$goodness >= (min.goodness + 5)) {
      optimization_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

      if (verbose) {
        cat("Original search succeeded in", round(optimization_time, 2), "seconds\n")
        cat("Best solution: rho =", round(best_original$rho, 3), ", psi =", round(best_original$psi, 3),
            ", ploidy =", round(best_original$ploidy, 3), ", goodness =", round(best_original$goodness, 2), "\n")
      }

      # Generate output using original result
      result <- finalize_battenberg_result(best_original, original_results, optimization_time,
                                          analysis, distancepng, copynumberprofilespng, nonroundedprofilepng,
                                          d, minimise, b, r, s, gamma, ch, lrr, bafsegmented,
                                          chr.names, reliabilityFile, cnaStatusFile, verbose)
      return(result)
    }

    if (verbose) {
      cat("Original search found solution (goodness=", round(best_original$goodness, 2),
          ") but quality is marginal. Trying enhanced search...\n")
    }
  } else {
    if (verbose) cat("Original search found no solutions. Trying enhanced search...\n")
  }

  # PHASE 2: Enhanced search only if needed
  if (!enable_robust_fallback) {
    if (verbose) cat("Robust fallback disabled. No solution found.\n")
    write.table("no copy number solutions found", file=cnaStatusFile, quote=F, col.names=F, row.names=F)
    return(list(psi = NA, rho = NA, ploidy = NA, convergence_info = list(converged = FALSE)))
  }

  if (verbose) cat("Phase 2: Enhanced targeted search...\n")

  # Setup parallel processing only if needed
  if (parallel) {
    hpc_info <- setup_parallel_environment(n_cores, verbose)
    n_cores <- hpc_info$n_cores
  }

  # Smart targeted search based on distance matrix analysis
  enhanced_results <- run_enhanced_targeted_search(d, s, gamma, min.ploidy, max.ploidy,
                                                   min.rho, max.rho, min.goodness,
                                                   TheoretMaxdist, minimise, allow100percent,
                                                   parallel, n_cores, verbose)

  # Combine results if we have both
  all_results <- c(original_results, enhanced_results)

  if (length(all_results) == 0) {
    if (verbose) cat("Phase 3: Fallback search with relaxed constraints...\n")

    # PHASE 3: Last resort with relaxed constraints
    fallback_results <- run_fallback_search(d, s, gamma, min.ploidy, max.ploidy,
                                           min.rho, max.rho, min.goodness,
                                           TheoretMaxdist, minimise, allow100percent, verbose)
    all_results <- fallback_results
  }

  # Final result processing
  if (length(all_results) == 0) {
    write.table("no copy number solutions found", file=cnaStatusFile, quote=F, col.names=F, row.names=F)
    if (verbose) cat("No suitable copy number solution found\n")
    return(list(psi = NA, rho = NA, ploidy = NA, convergence_info = list(converged = FALSE)))
  }

  # Select best result
  best_result <- all_results[[which.max(sapply(all_results, function(x) x$goodness))]]
  optimization_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  if (verbose) {
    cat("Enhanced search completed in", round(optimization_time, 2), "seconds\n")
    cat("Best solution: rho =", round(best_result$rho, 3), ", psi =", round(best_result$psi, 3),
        ", ploidy =", round(best_result$ploidy, 3), ", goodness =", round(best_result$goodness, 2), "\n")
  }

  result <- finalize_battenberg_result(best_result, all_results, optimization_time,
                                      analysis, distancepng, copynumberprofilespng, nonroundedprofilepng,
                                      d, minimise, b, r, s, gamma, ch, lrr, bafsegmented,
                                      chr.names, reliabilityFile, cnaStatusFile, verbose)
  return(result)
}

#' Run original Battenberg search algorithm exactly
run_original_battenberg_search <- function(d, s, gamma, min.ploidy, max.ploidy, min.rho, max.rho,
                                          min.goodness, TheoretMaxdist, minimise, allow100percent, verbose) {

  results <- list()

  # Original algorithm: 7x7 local search
  for (i in 4:(dim(d)[1]-3)) {
    for (j in 4:(dim(d)[2]-3)) {
      m = d[i,j]
      seld = d[(i-3):(i+3),(j-3):(j+3)]
      seld[4,4] = max(seld)

      if(min(seld) > m) {
        psi = as.numeric(rownames(d)[i])
        rho = as.numeric(colnames(d)[j])

        # Calculate solution
        solution <- calculate_battenberg_solution(psi, rho, s, gamma, min.ploidy, max.ploidy,
                                                 min.rho, max.rho, min.goodness, m,
                                                 TheoretMaxdist, minimise, allow100percent)

        if (!is.null(solution)) {
          results[[length(results) + 1]] <- solution
        }
      }
    }
  }

  # Handle 100% aberrant case if needed (original logic)
  if (allow100percent && length(results) == 0) {
    cold = which(as.numeric(colnames(d)) > 1)
    d[,cold] = 1E20

    for (i in 4:(dim(d)[1]-3)) {
      for (j in 4:(dim(d)[2]-3)) {
        m = d[i,j]
        seld = d[(i-3):(i+3),(j-3):(j+3)]
        seld[4,4] = max(seld)

        if(min(seld) > m) {
          psi = as.numeric(rownames(d)[i])
          rho = as.numeric(colnames(d)[j])

          solution <- calculate_battenberg_solution(psi, rho, s, gamma, min.ploidy, max.ploidy,
                                                   min.rho, max.rho, min.goodness, m,
                                                   TheoretMaxdist, minimise, allow100percent,
                                                   skip_zero_check = TRUE)

          if (!is.null(solution)) {
            results[[length(results) + 1]] <- solution
          }
        }
      }
    }
  }

  if (verbose && length(results) > 0) {
    cat("  Original search found", length(results), "solutions\n")
  }

  return(results)
}

#' Enhanced targeted search - focuses on most promising regions
run_enhanced_targeted_search <- function(d, s, gamma, min.ploidy, max.ploidy, min.rho, max.rho,
                                        min.goodness, TheoretMaxdist, minimise, allow100percent,
                                        parallel, n_cores, verbose) {

  # Identify promising regions from distance matrix
  finite_d <- d[is.finite(d)]
  if (length(finite_d) == 0) return(list())

  # Focus on top 10% of distance values
  threshold <- quantile(finite_d, 0.1, na.rm = TRUE)
  promising_indices <- which(d <= threshold & is.finite(d), arr.ind = TRUE)

  if (nrow(promising_indices) == 0) return(list())

  # Limit to reasonable number of points for speed
  max_points <- min(50, nrow(promising_indices))
  if (nrow(promising_indices) > max_points) {
    sample_indices <- sample(nrow(promising_indices), max_points)
    promising_indices <- promising_indices[sample_indices, ]
  }

  if (verbose) cat("  Testing", nrow(promising_indices), "promising regions\n")

  # Optimize each promising point
  optimize_point <- function(idx) {
    tryCatch({
      row_idx <- promising_indices[idx, 1]
      col_idx <- promising_indices[idx, 2]

      psi <- as.numeric(rownames(d)[row_idx])
      rho <- as.numeric(colnames(d)[col_idx])
      m <- d[row_idx, col_idx]

      # Enhanced local search around this point
      best_solution <- NULL
      best_distance <- Inf

      # Search 5x5 neighborhood
      for (di in -2:2) {
        for (dj in -2:2) {
          ni <- row_idx + di
          nj <- col_idx + dj

          if (ni >= 1 && ni <= nrow(d) && nj >= 1 && nj <= ncol(d) && is.finite(d[ni, nj])) {
            psi_test <- as.numeric(rownames(d)[ni])
            rho_test <- as.numeric(colnames(d)[nj])
            m_test <- d[ni, nj]

            solution <- calculate_battenberg_solution(psi_test, rho_test, s, gamma, min.ploidy, max.ploidy,
                                                     min.rho, max.rho, min.goodness, m_test,
                                                     TheoretMaxdist, minimise, allow100percent)

            if (!is.null(solution) && solution$distance < best_distance) {
              best_solution <- solution
              best_distance <- solution$distance
            }
          }
        }
      }

      return(best_solution)
    }, error = function(e) {
      return(NULL)
    })
  }

  # Run optimization
  if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    tryCatch({
      solutions <- future.apply::future_lapply(1:nrow(promising_indices), optimize_point, future.seed = TRUE)
      solutions <- solutions[!sapply(solutions, is.null)]
    }, error = function(e) {
      if (verbose) cat("  Parallel processing failed, using sequential\n")
      solutions <- lapply(1:nrow(promising_indices), optimize_point)
      solutions <- solutions[!sapply(solutions, is.null)]
    })
  } else {
    solutions <- lapply(1:nrow(promising_indices), optimize_point)
    solutions <- solutions[!sapply(solutions, is.null)]
  }

  if (verbose && length(solutions) > 0) {
    cat("  Enhanced search found", length(solutions), "solutions\n")
  }

  return(solutions)
}

#' Fallback search with relaxed constraints
run_fallback_search <- function(d, s, gamma, min.ploidy, max.ploidy, min.rho, max.rho,
                               min.goodness, TheoretMaxdist, minimise, allow100percent, verbose) {

  # More relaxed constraints
  relaxed_min.ploidy <- max(1.0, min.ploidy - 0.3)
  relaxed_max.ploidy <- min(6.0, max.ploidy + 0.5)
  relaxed_min.rho <- max(0.05, min.rho - 0.05)
  relaxed_max.rho <- min(1.0, max.rho + 0.05)
  relaxed_min.goodness <- max(30, min.goodness - 25)

  if (verbose) {
    cat("  Using relaxed constraints: ploidy=[", round(relaxed_min.ploidy, 2), ",", round(relaxed_max.ploidy, 2),
        "], rho=[", round(relaxed_min.rho, 2), ",", round(relaxed_max.rho, 2),
        "], goodness>=", round(relaxed_min.goodness, 1), "\n")
  }

  # Test a broader but still targeted set of points
  rho_values <- as.numeric(colnames(d))
  psi_values <- as.numeric(rownames(d))

  # Create strategic grid
  test_rhos <- seq(relaxed_min.rho, relaxed_max.rho, length.out = 20)
  test_psis <- seq(relaxed_min.ploidy, relaxed_max.ploidy, length.out = 20)

  results <- list()

  for (test_rho in test_rhos) {
    for (test_psi in test_psis) {
      rho_idx <- which.min(abs(rho_values - test_rho))
      psi_idx <- which.min(abs(psi_values - test_psi))

      if (is.finite(d[psi_idx, rho_idx])) {
        m <- d[psi_idx, rho_idx]

        solution <- calculate_battenberg_solution(test_psi, test_rho, s, gamma,
                                                 relaxed_min.ploidy, relaxed_max.ploidy,
                                                 relaxed_min.rho, relaxed_max.rho,
                                                 relaxed_min.goodness, m,
                                                 TheoretMaxdist, minimise, allow100percent)

        if (!is.null(solution)) {
          results[[length(results) + 1]] <- solution
        }
      }
    }
  }

  if (verbose && length(results) > 0) {
    cat("  Fallback search found", length(results), "solutions\n")
  }

  return(results)
}

#' Calculate Battenberg solution with error handling
calculate_battenberg_solution <- function(psi, rho, s, gamma, min.ploidy, max.ploidy, min.rho, max.rho,
                                         min.goodness, distance_value, TheoretMaxdist, minimise,
                                         allow100percent, skip_zero_check = FALSE) {

  tryCatch({
    # Validate inputs
    if (is.na(psi) || is.na(rho) || !is.finite(psi) || !is.finite(rho) || rho <= 0) {
      return(NULL)
    }

    # Calculate copy numbers
    nA <- (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
    nB <- (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho

    if (any(is.na(nA)) || any(is.na(nB)) || any(!is.finite(nA)) || any(!is.finite(nB))) {
      return(NULL)
    }

    # Calculate ploidy
    ploidy <- sum((nA+nB) * s[,"length"]) / sum(s[,"length"])

    if (is.na(ploidy) || !is.finite(ploidy)) {
      return(NULL)
    }

    # Calculate goodness of fit
    if(minimise) {
      goodnessOfFit <- (1 - distance_value/TheoretMaxdist) * 100
    } else {
      goodnessOfFit <- -distance_value/TheoretMaxdist * 100
    }

    if (is.na(goodnessOfFit) || !is.finite(goodnessOfFit)) {
      return(NULL)
    }

    # Check constraints
    if (ploidy < min.ploidy || ploidy > max.ploidy ||
        rho < min.rho || rho > max.rho ||
        goodnessOfFit < min.goodness) {
      return(NULL)
    }

    # Check zero conditions unless skipped
    if (!skip_zero_check && !allow100percent) {
      percentzero <- (sum((round(nA)==0)*s[,"length"])+sum((round(nB)==0)*s[,"length"]))/sum(s[,"length"])
      perczeroAbb <- (sum((round(nA)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1))+sum((round(nB)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1)))/sum(s[,"length"]*ifelse(s[,"b"]==0.5,0,1))

      if (is.na(perczeroAbb)) perczeroAbb <- 0

      if (!(percentzero > 0.01 || perczeroAbb > 0.1)) {
        return(NULL)
      }
    }

    return(list(
      psi = psi,
      rho = min(rho, 1.0),
      ploidy = ploidy,
      goodness = goodnessOfFit,
      distance = distance_value
    ))

  }, error = function(e) {
    return(NULL)
  })
}

#' Finalize result and generate outputs
finalize_battenberg_result <- function(best_result, all_results, optimization_time,
                                      analysis, distancepng, copynumberprofilespng, nonroundedprofilepng,
                                      d, minimise, b, r, s, gamma, ch, lrr, bafsegmented,
                                      chr.names, reliabilityFile, cnaStatusFile, verbose) {

  psi_opt1 <- best_result$psi
  rho_opt1 <- best_result$rho
  ploidy_opt1 <- best_result$ploidy
  goodnessOfFit_opt1 <- best_result$goodness

  # Create convergence information
  n_solutions <- length(all_results)
  convergence_info <- list(
    converged = n_solutions > 0,
    n_solutions_found = n_solutions,
    optimization_time = optimization_time,
    method_used = if(n_solutions == 1) "original" else "enhanced"
  )

  write.table(paste(n_solutions, "solutions found"),
              file=cnaStatusFile, quote=F, col.names=F, row.names=F)

  # Generate plots (same as original)
  generate_plots_battenberg(analysis, distancepng, copynumberprofilespng, nonroundedprofilepng,
                           d, psi_opt1, rho_opt1, ploidy_opt1, goodnessOfFit_opt1, minimise,
                           b, r, s, gamma, ch, lrr, bafsegmented, chr.names, reliabilityFile)

  return(list(
    psi = psi_opt1,
    rho = rho_opt1,
    ploidy = ploidy_opt1,
    convergence_info = convergence_info
  ))
}

#' Helper functions for parallel processing
setup_parallel_environment <- function(n_cores, verbose = TRUE) {
  hpc_detected <- FALSE
  hpc_scheduler <- "unknown"
  
  if (Sys.getenv("SLURM_JOB_ID", unset = "") != "") {
    hpc_detected <- TRUE
    hpc_scheduler <- "SLURM"
    if (is.null(n_cores)) {
      slurm_cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
      n_cores <- max(1, slurm_cpus)
    }
  } else if (Sys.getenv("PBS_JOBID", unset = "") != "") {
    hpc_detected <- TRUE
    hpc_scheduler <- "PBS/Torque"
    if (is.null(n_cores)) {
      pbs_ncpus <- as.numeric(Sys.getenv("PBS_NCPUS", unset = "1"))
      n_cores <- max(1, pbs_ncpus)
    }
  }
  
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  if (verbose && hpc_detected) {
    cat("HPC environment detected:", hpc_scheduler, "with", n_cores, "cores\n")
  }
  
  return(list(n_cores = n_cores, hpc_detected = hpc_detected, scheduler = hpc_scheduler))
}

#' Generate plots
generate_plots_battenberg <- function(analysis, distancepng, copynumberprofilespng, nonroundedprofilepng,
                                     d, psi_opt1, rho_opt1, ploidy_opt1, goodnessOfFit_opt1, minimise,
                                     b, r, s, gamma, ch, lrr, bafsegmented, chr.names, reliabilityFile) {
  
  if (analysis == "paired") {
    psi_opt1_plot <- psi_opt1
    rho_opt1_plot <- rho_opt1
    
    if (!is.na(distancepng)) {
      png(filename = distancepng, width = 1000, height = 1000, res = 1000/7, type = "cairo")
    }
    ASCAT::ascat.plotSunrise(-d, psi_opt1_plot, rho_opt1_plot, minimise)
    if (!is.na(distancepng)) { dev.off() }
  }
  
  nAfull <- (rho_opt1-1-(b-1)*2^(r/gamma)*((1-rho_opt1)*2+rho_opt1*psi_opt1))/rho_opt1
  nBfull <- (rho_opt1-1+b*2^(r/gamma)*((1-rho_opt1)*2+rho_opt1*psi_opt1))/rho_opt1
  nA <- pmax(round(nAfull),0)
  nB <- pmax(round(nBfull),0)
  
  if(!is.na(reliabilityFile)){
    rBacktransform <- gamma*log((rho_opt1*(nA+nB)+(1-rho_opt1)*2)/((1-rho_opt1)*2+rho_opt1*psi_opt1),2)
    bBacktransform <- (1-rho_opt1+rho_opt1*nB)/(2-2*rho_opt1+rho_opt1*(nA+nB))
    rConf <- ifelse(abs(rBacktransform)>0.15,pmin(100,pmax(0,100*(1-abs(rBacktransform-r)/abs(r)))),NA)
    bConf <- ifelse(bBacktransform!=0.5,pmin(100,pmax(0,ifelse(b==0.5,100,100*(1-abs(bBacktransform-b)/abs(b-0.5))))),NA)
    
    write.table(data.frame(segmentedBAF=b,backTransformedBAF=bBacktransform,confidenceBAF=bConf,
                          segmentedR=r,backTransformedR=rBacktransform,confidenceR=rConf,
                          nA=nA,nB=nB,nAfull=nAfull,nBfull=nBfull), 
                reliabilityFile, sep=",", row.names=F)
  }
  
  if (!is.na(copynumberprofilespng)) {
    png(filename = copynumberprofilespng, width = 2000, height = 500, res = 200, type = "cairo")
  }
  ASCAT::ascat.plotAscatProfile(n1all = nA, n2all = nB, heteroprobes = TRUE, 
                               ploidy = ploidy_opt1, rho = rho_opt1, goodnessOfFit = goodnessOfFit_opt1, 
                               nonaberrant = FALSE, ch = ch, lrr = lrr, bafsegmented = bafsegmented, 
                               chrs = chr.names)
  if (!is.na(copynumberprofilespng)) { dev.off() }
  
  if (!is.na(nonroundedprofilepng)) {
    png(filename = nonroundedprofilepng, width = 2000, height = 500, res = 200, type = "cairo")
  }
  ASCAT::ascat.plotNonRounded(ploidy = ploidy_opt1, rho = rho_opt1, goodnessOfFit = goodnessOfFit_opt1, 
                             nonaberrant = FALSE, nAfull = nAfull, nBfull = nBfull, 
                             bafsegmented = bafsegmented, ch = ch, lrr = lrr, chrs = chr.names)
  if (!is.na(nonroundedprofilepng)) { dev.off() }
}
