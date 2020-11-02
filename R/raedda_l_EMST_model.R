raedda_l_EMST_model <- function(X_train,
                                class_train,
                                n_relevant_variables,
                                alpha_train,
                                # The proportion of obs in the Xtrain to be trimmed
                                model_name,
                                swap_step,
                                ctrl_init,
                                ctrl_GA,
                                ...) {
  N_train <- nrow(X_train)
  ltrain <- mclust::unmap(class_train)
  D <- ncol(X_train)
  G <- ncol(ltrain)
  class_train <- as.factor(class_train)
  classLabel <- levels(class_train)

  # Core algorithm ----------------------------------------------------------


  if (alpha_train != 0) {
    nsamp <- ctrl_init$n_samp
    max_iter_init <- ctrl_init$max_iter

    N_train_trim <- N_train - ceiling(N_train * alpha_train)
    robust_result <-
      patterned_EMST(
        nsamp = nsamp,
        # I perform the estimation starting from nsamp J_g subsets of sample size (d+1), inspired to what done in \cite{Hubert2018}
        X_train = X_train,
        class_train = class_train,
        ltrain = ltrain,
        alpha_train = alpha_train,
        n_relevant_variables = n_relevant_variables,
        model_name = model_name,
        swap_step = swap_step,
        ctrl_GA = ctrl_GA,
        max_iter_init = max_iter_init
      )
    ll <- robust_result$ll
    fitm <- robust_result$fitm
    F_subset <- robust_result$relevant_variables
    E_subset <- robust_result$irrelevant_variables
    m_E_F <- robust_result$m_E_F
    V_E_F <- robust_result$V_E_F
    G_E_F <- robust_result$G_E_F

  } else{
    N_train_trim <- NULL
    # M step: computed for the whole set of columns D

    fitm <-
      mclust::mstep(modelName = model_name,
                    data = X_train,
                    z = ltrain)
    S_R_w <- fitm$parameters$variance$sigma
    fitm_one_group <-
      mclust::mstep(modelName = model_name,
                    data = X_train,
                    z = rep(1, nrow(X_train)))
    # fitm_one_group is used afterwards for updating the regression parameters
    S_R <- fitm_one_group$parameters$variance$sigma[, , 1]

    # S step
    if (model_name == "EEI") {
      F_position <-
        order(diag(S_R_w[, , 1]) / diag(S_R)) # (5.21 b) Ritter pag 209
      F_subset <- (1:D)[F_position <= n_relevant_variables]
    } else if (model_name == "VVI") {
      F_position <-
        order(log(apply(S_R_w, 3, diag) / diag(S_R)) %*% fitm$parameters$pro) # (5.21) Ritter pag 209
      F_subset <- (1:D)[F_position <= n_relevant_variables]
    } else {
      if (swap_step == "exhaustive") {
        H_F_value <-
          apply(
            possible_col_subset,
            1,
            H_F_subset,
            S_R_w = S_R_w,
            S_R = S_R,
            fitm = fitm
          )
        F_subset <-
          possible_col_subset[which.min(H_F_value), , drop = TRUE]
      } else if (swap_step == "ga") {
        swap_step_GA <- kofnGA::kofnGA(
          n = D,
          k = n_relevant_variables,
          popsize = ctrl_GA$popsize,
          keepbest = ctrl_GA$keepbest,
          ngen = ctrl_GA$ngen,
          tourneysize = ctrl_GA$tourneysize,
          mutprob = ctrl_GA$mutprob,
          mutfrac = ctrl_GA$mutfrac,
          initpop = ctrl_GA$initpop,
          cluster = ctrl_GA$cluster,
          sharedmemory = ctrl_GA$sharedmemory,
          OF = H_F_subset,
          S_R_w = S_R_w,
          S_R = S_R,
          fitm = fitm,
          verbose = 0
        )
        F_subset <- swap_step_GA$bestsol
      }
    }
    E_subset <- setdiff(1:D, F_subset)

    # T-step (used for computing log-lik afterwards)

    G_E_F <-
      S_R[E_subset, F_subset, drop = FALSE] %*% MASS::ginv(S_R[F_subset, F_subset, drop =
                                                                 FALSE])
    m_E_F <-
      as.vector(fitm_one_group$parameters$mean[E_subset, , drop = FALSE] - G_E_F %*%
                  fitm_one_group$parameters$mean[F_subset, , drop = FALSE])
    V_E_F <-
      S_R[E_subset, E_subset, drop = FALSE] - S_R[E_subset, F_subset, drop = FALSE] %*%
      MASS::ginv(S_R[F_subset, F_subset, drop = FALSE]) %*%
      S_R[F_subset, E_subset, drop = FALSE]

  }
  # Checking if errors in the procedure -------------------------------------

  fitetrain <-
    tryCatch(
      estep_EMST(
        data = X_train,
        fitm = fitm,
        F_subset = F_subset
      ),
      error = function(e) {
        list(z = NA)
      }
    )
  emptyz <- ifelse(all(!is.na(fitetrain$z)), yes = FALSE, no = TRUE)

  # Results Collection ------------------------------------------------------

  if (!emptyz) {
    res <- list()
    res$N_train <- N_train
    res$N_train_after_trimming <- N_train_trim
    res$alpha_train <- alpha_train
    res$d <- D
    res$relevant_variables <- F_subset
    res$irrelevant_variables <- E_subset
    res$G <- G
    res$model_name <- model_name
    res$parameters <- fitm$parameters
    #Name the different groups
    names(res$parameters$pro) <- classLabel
    if (D == 1) {
      names(res$parameters$mean) <- classLabel
    } else {
      colnames(res$parameters$mean) <- classLabel
    }

    ztrain <- fitetrain$z
    cltrain <-
      factor(sapply(map(ztrain), function(i)
        classLabel[i]), levels = classLabel) # I classify a posteriori also the trimmed units
    pos_trimmed_train <- NULL
    cltrain_after_trimming <- NULL

    grouping_component_ll <-
      tryCatch(
        sapply(1:fitm$G, function(g)
          log(fitm$parameters$pro[g]) + mvnfast::dmvn(
            X = as.matrix(X_train[, F_subset]),
            mu = fitm$parameters$mean[F_subset, g],
            sigma = fitm$parameters$variance$sigma[F_subset, F_subset , g],
            log = T
          )),
        error = function(e)
          - Inf,
        warning = function(w)
          - Inf
      )


    ind_D_Xtrain_cdens <-
      cbind(1:N_train, mclust::map(ltrain)) # matrix with obs and respective group from ltrain

    grouping_ll <-
      grouping_component_ll[ind_D_Xtrain_cdens] # I compute D_g conditioning ON the fact that I know the supposed true class

    eigen_V_E_F <-
      eigen(V_E_F, symmetric = TRUE, only.values = TRUE)$values # I force symmetric=TRUE to avoid numerical problems

    eigen_tol <-  1e-12 # FIX: maybe function arg?

    is_singular_V_E_F <- any(eigen_V_E_F < eigen_tol)

    if (is_singular_V_E_F) {
      # when D is large residual scatter matrix may be singular.
      # In such a case, the data still allows us to estimate a singular normal distribution on a subspace using the generalized inverse (pag 208 Ritter2015)
      no_grouping_ll <-
        dsmvnorm(
          x = as.matrix(X_train[, E_subset, drop = FALSE]) - as.matrix(X_train[, F_subset, drop = FALSE]) %*%
            t(G_E_F),
          mean = m_E_F,
          sigma = V_E_F,
          log = T,
          eigen_sigma = eigen_V_E_F,
          eigen_tol = eigen_tol
        )
    } else {
      no_grouping_ll <-
        mvnfast::dmvn(
          X = as.matrix(X_train[, E_subset, drop = FALSE]) - as.matrix(X_train[, F_subset, drop = FALSE]) %*%
            t(G_E_F),
          mu = m_E_F,
          sigma = V_E_F,
          log = T
        )
    }

    unit_likelihood <- grouping_ll + no_grouping_ll
    ll <- sum(unit_likelihood)

    if (alpha_train != 0) {
      pos_trimmed_train <-
        which(unit_likelihood <= (sort(unit_likelihood, decreasing = F)[[ceiling(N_train * alpha_train)]]))

      cltrain_after_trimming <- cltrain
      cltrain_after_trimming <-
        factor(cltrain, levels = c(classLabel, "0"))
      cltrain_after_trimming[pos_trimmed_train] <- "0"
        ll <- sum(unit_likelihood[-pos_trimmed_train])
    }

    res$train <- list()
    res$train$z <- ztrain
    res$train$cl <- cltrain
    res$train$cl_after_trimming <- cltrain_after_trimming
    res$train$obs_trimmed <- pos_trimmed_train
    res$train$alpha_train <- alpha_train
    bic_all <-
      2 * ll - (
        mclust::nMclustParams(
          # parameters for the Grouping part
          modelName = model_name,
          d = n_relevant_variables,
          G = G,
          noise = FALSE,
          equalPro = FALSE
        ) + mclust::nMclustParams(
          # parameters for the No Grouping part
          modelName = model_name,
          d = D,
          G = 1,
          noise = FALSE,
          equalPro = FALSE
        )
      ) * log(fitm$n) # if alpha_train !=0 this is trimmed BIC
    res$ll <- ll
    res$bic <- bic_all
  }
  else {
    res <- list()
    res$error <-
      "Either training groups too small or trimming level too large for this model type"
    res$ll <- NA
    res$bic <- NA
  }
  res
}
