patterned_EMST_new <- # in this version I do not consider the regression part
  # main function that will be repeated nsamp times, run only if alpha_train!=0
  function(nsamp,
           X_train,
           class_train,
           ltrain,
           alpha_train,
           n_relevant_variables,
           model_name,
           swap_step,
           ctrl_GA,
           max_iter_init) {
    robust_start <-
      function(X_train,
               class_train,
               ltrain,
               alpha_train,
               n_relevant_variables,
               model_name,
               swap_step,
               max_iter_init,
               ctrl_GA) {
        Ntrain <- nrow(X_train)
        D <- ncol(X_train)
        J_ind <-
          tryCatch(
            c(sapply(
              split(seq(class_train), class_train),
              sample,
              size = D + 1,
              replace = FALSE
            )),
            error = function(e)
              seq(class_train)
          )

        #Init with the random (p + 1)-subset J_g
        fitm <- tryCatch(
          mclust::mstep(
            modelName = model_name,
            data = X_train[J_ind, , drop = F],
            z = ltrain[J_ind, ]
          ),
          error = function(e) {
            list(parameters = NA)
          }
        )
        X_train_fit <-
          NULL # I initialize its value for preventing erros when the EM algorithm fails
        emptyz <- TRUE
        F_subset <- sample(1:D, size = n_relevant_variables, replace = F) # random select n_selected_variables columns for initialization
        if(swap_step=="exhaustive"){
          possible_col_subset <- gtools::combinations(n = D, # D= ncol(X_train)
                                                      r = n_relevant_variables
          ) # this is used in the swap step
        }

        # Initializing R^{(0)} using only the grouping info and the conditional group densities:

        # grouping_component_ll <-
        #   tryCatch(
        #     sapply(1:fitm$G, function(g)
        #       log(fitm$parameters$pro[g]) + mvnfast::dmvn(
        #         X = as.matrix(X_train[,F_subset]), # does not work with data.frame objects
        #         mu = fitm$parameters$mean[F_subset, g],
        #         sigma = fitm$parameters$variance$sigma[F_subset,F_subset , g],
        #         log = T
        #       )),
        #     error = function(e)
        #       - Inf,
        #     warning = function(w)
        #       - Inf
        #   )

        D_Xtrain_cond <- tryCatch(
          sapply(1:fitm$G, function(g)
            mvnfast::dmvn(
              X = as.matrix(X_train[, F_subset]),
              # does not work with data.frame objects
              mu = fitm$parameters$mean[F_subset, g],
              sigma = fitm$parameters$variance$sigma[F_subset, F_subset , g],
              log = T
            ), simplify = TRUE),
          error = function(e) {
              NA
          }
        )
        if(any(is.na(D_Xtrain_cond))){
          return(list(
            ll = NA,
            fitm = NA,
            relevant_variables = NA
          )) # in any of any singular sigma[F_subset, F_subset ,] return ll=NA and exit the function
        }
        ind_D_Xtrain_cdens <-
          cbind(1:Ntrain, mclust::map(ltrain)) # matrix with obs and respective group from ltrain
        D_Xtrain <-
          D_Xtrain_cond[ind_D_Xtrain_cdens] # I compute D_g conditioning ON the fact that I know the supposed true class
        pos_trimmed_train <-
          which(D_Xtrain <= (sort(D_Xtrain, decreasing = F)
                             [[ceiling(Ntrain * alpha_train)]]))

        # comp_density <- sapply(1:fitm$G, function(g)
        #   mvnfast::dmvn(
        #     X = as.matrix(X_train[, F_subset]),
        #     # does not work with data.frame objects
        #     mu = fitm$parameters$mean[F_subset, g],
        #     sigma = fitm$parameters$variance$sigma[F_subset, F_subset , g],
        #     log = T
        #   ), simplify = TRUE)
        #
        # grouping_ll <-
        #   grouping_component_ll[ind_D_X_train_cdens] # I compute D_g conditioning ON the fact that I know the supposed true class
        #
        # pos_trimmed_train <-
        #   which(grouping_ll <= (sort(grouping_ll, decreasing = F)[[ceiling(Ntrain * alpha_train)]]))

        if (length(pos_trimmed_train) != ceiling(Ntrain * alpha_train)) {
          #condition due to the fact that for some models (the simplest ones usually) it might happen that 2 obs have exactly the same comp prob, leading to errors
          pos_trimmed_train <-
            pos_trimmed_train[1:ceiling(Ntrain * alpha_train)]
        }
        pos_old <- pos_trimmed_train
        ltrain_fit <- ltrain[-pos_trimmed_train, , drop = F]
        X_train_fit <- X_train[-pos_trimmed_train, , drop = F]
        # the estimation will continue running until for two consecutive iterations exactly
        # the same obs are trimmed
        criterion <- TRUE
        iter_init <- 0

        while (criterion) {

          iter_init <- iter_init + 1

          # M step: computed for the whole set of columns D

          fitm <-
            mclust::mstep(modelName = model_name,
                          data = X_train_fit,
                          z = ltrain_fit)
          S_R_w <- fitm$parameters$variance$sigma
          fitm_one_group <-
            mclust::mstep(modelName = model_name,
                          data = X_train_fit,
                          z = rep(1, nrow(X_train_fit)))
          # fitm_one_group is used afterwards for updating the regression parameters
          S_R <- fitm_one_group$parameters$variance$sigma[, , 1]

          # S step
          if(model_name=="EEI") {
            F_position <- order(diag(S_R_w[,,1])/diag(S_R)) # (5.21 b) Ritter pag 209
            F_subset <- (1:D)[F_position<=n_relevant_variables]
          } else if(model_name == "VVI") {
            F_position <-
              order(log(apply(S_R_w, 3, diag)/ diag(S_R)) %*% fitm$parameters$pro) # (5.21) Ritter pag 209
            F_subset <- (1:D)[F_position<=n_relevant_variables]
          } else {
            if(swap_step=="exhaustive"){
              H_F_value <- apply(possible_col_subset, 1, H_F_subset, S_R_w=S_R_w, S_R=S_R, fitm = fitm)
              F_subset <-
                possible_col_subset[which.min(H_F_value), , drop = TRUE]
            } else if (swap_step=="ga"){
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

          # E_subset <- setdiff(1:D, F_subset)
          #
          # # T step
          #
          # G_E_F <-
          #   S_R[E_subset, F_subset, drop = FALSE] %*% MASS::ginv(S_R[F_subset, F_subset, drop =
          #                                                              FALSE])
          # m_E_F <-
          #   as.vector(fitm_one_group$parameters$mean[E_subset, , drop = FALSE] - G_E_F %*%
          #               fitm_one_group$parameters$mean[F_subset, , drop = FALSE])
          # V_E_F <-
          #   S_R[E_subset, E_subset, drop = FALSE] - S_R[E_subset, F_subset, drop = FALSE] %*%
          #   MASS::ginv(S_R[F_subset, F_subset, drop = FALSE]) %*%
          #   S_R[F_subset, E_subset, drop = FALSE]

          # compute the ll for the entire model

          D_Xtrain_cond <- tryCatch(
            sapply(1:fitm$G, function(g)
              mvnfast::dmvn(
                X = as.matrix(X_train[, F_subset]),
                # does not work with data.frame objects
                mu = fitm$parameters$mean[F_subset, g],
                sigma = fitm$parameters$variance$sigma[F_subset, F_subset , g],
                log = T
              ), simplify = TRUE),
            error = function(e) {
              NA
            }
          )
          if(any(is.na(D_Xtrain_cond))){
            return(list(
              ll = NA,
              fitm = NA,
              relevant_variables = NA
            )) # in any of any singular sigma[F_subset, F_subset ,] return ll=NA and exit the function
          }

          ind_D_Xtrain_cdens <-
            cbind(1:Ntrain, mclust::map(ltrain)) # matrix with obs and respective group from ltrain
          D_Xtrain <-
            D_Xtrain_cond[ind_D_Xtrain_cdens] # I compute D_g conditioning ON the fact that I know the supposed true class
          pos_trimmed_train <-
            which(D_Xtrain <= (sort(D_Xtrain, decreasing = F)
                               [[ceiling(Ntrain * alpha_train)]]))

          # eigen_V_E_F <-
          #   eigen(V_E_F, symmetric = TRUE, only.values = TRUE)$values # I force symmetric=TRUE to avoid numerical problems
          #
          # eigen_tol <-  1e-12 # FIX: maybe function arg?
          #
          # is_singular_V_E_F <- any(eigen_V_E_F < eigen_tol)
          #
          # if(is_singular_V_E_F){ # when D is large residual scatter matrix may be singular.
          #   # In such a case, the data still allows us to estimate a singular normal distribution on a subspace using the generalized inverse (pag 208 Ritter2015)
          #   no_grouping_ll <-
          #     dsmvnorm(
          #       x = as.matrix(X_train[, E_subset, drop = FALSE]) - as.matrix(X_train[, F_subset, drop = FALSE]) %*%
          #         t(G_E_F),
          #       mean = m_E_F,
          #       sigma = V_E_F,
          #       log = T,
          #       eigen_sigma = eigen_V_E_F,
          #       eigen_tol = eigen_tol
          #     )
          # } else {
          #   no_grouping_ll <-
          #     mvnfast::dmvn(
          #       X = as.matrix(X_train[, E_subset, drop = FALSE]) - as.matrix(X_train[, F_subset, drop = FALSE]) %*%
          #         t(G_E_F),
          #       mu = m_E_F,
          #       sigma = V_E_F,
          #       log = T
          #     )
          # }

          # unit_likelihood <- grouping_ll + no_grouping_ll
          # unit_likelihood <- grouping_ll
          # pos_trimmed_train <-
          #   which(unit_likelihood <= (sort(unit_likelihood, decreasing = F)[[ceiling(Ntrain * alpha_train)]]))

          if (length(pos_trimmed_train) != ceiling(Ntrain * alpha_train)) {
            #condition due to the fact that for some models (the simplest ones usually) it might happen that 2 obs have exactly the same comp prob, leading to errors
            pos_trimmed_train <-
              pos_trimmed_train[1:ceiling(Ntrain * alpha_train)]
          }

          ltrain_fit <- ltrain[-pos_trimmed_train, , drop = F]
          X_train_fit <- X_train[-pos_trimmed_train, , drop = F]

          # FIX this piece needs to be deleted: just for checking monot
          # lDensity_train <-
          #   sapply(1:fitm$G, function(g)
          #     mvnfast::dmvn(
          #       X = as.matrix(X_train_fit[,F_subset]),
          #       mu = fitm$parameters$mean[F_subset, g],
          #       sigma = fitm$parameters$variance$sigma[F_subset, F_subset, g],
          #       log = T
          #     ))
          #
          # E_subset <- setdiff(1:D, F_subset)
          #
          # # T step
          #
          # G_E_F <-
          #   S_R[E_subset, F_subset, drop = FALSE] %*% MASS::ginv(S_R[F_subset, F_subset, drop =
          #                                                              FALSE])
          # m_E_F <-
          #   as.vector(fitm_one_group$parameters$mean[E_subset, , drop = FALSE] - G_E_F %*%
          #               fitm_one_group$parameters$mean[F_subset, , drop = FALSE])
          # V_E_F <-
          #   S_R[E_subset, E_subset, drop = FALSE] - S_R[E_subset, F_subset, drop = FALSE] %*%
          #   MASS::ginv(S_R[F_subset, F_subset, drop = FALSE]) %*%
          #   S_R[F_subset, E_subset, drop = FALSE]
          #
          # eigen_V_E_F <-
          #   eigen(V_E_F, symmetric = TRUE, only.values = TRUE)$values # I force symmetric=TRUE to avoid numerical problems
          #
          # eigen_tol <-  1e-12 # FIX: maybe function arg?
          #
          # is_singular_V_E_F <- any(eigen_V_E_F < eigen_tol)
          #
          # if(is_singular_V_E_F){ # when D is large residual scatter matrix may be singular.
          #   # In such a case, the data still allows us to estimate a singular normal distribution on a subspace using the generalized inverse (pag 208 Ritter2015)
          #   no_grouping_ll <-
          #     dsmvnorm(
          #       x = as.matrix(X_train[, E_subset, drop = FALSE]) - as.matrix(X_train[, F_subset, drop = FALSE]) %*%
          #         t(G_E_F),
          #       mean = m_E_F,
          #       sigma = V_E_F,
          #       log = T,
          #       eigen_sigma = eigen_V_E_F,
          #       eigen_tol = eigen_tol
          #     )
          # } else {
          #   no_grouping_ll <-
          #     mvnfast::dmvn(
          #       X = as.matrix(X_train[, E_subset, drop = FALSE]) - as.matrix(X_train[, F_subset, drop = FALSE]) %*%
          #         t(G_E_F),
          #       mu = m_E_F,
          #       sigma = V_E_F,
          #       log = T
          #     )
          # }
          #
          # mTau_train <-
          #   matrix(log(fitm$parameters$pro),
          #          nrow(X_train_fit),
          #          fitm$G,
          #          byrow = TRUE)
          #
          # lDensity_train <-
          #   sapply(1:fitm$G, function(g)
          #     mvnfast::dmvn(
          #       X = as.matrix(X_train_fit[,F_subset, drop=FALSE]),
          #       mu = fitm$parameters$mean[F_subset, g],
          #       sigma = fitm$parameters$variance$sigma[F_subset, F_subset, g],
          #       log = T
          #     ))
          #
          # sum_train <- mTau_train + lDensity_train
          # mat_train <- ltrain_fit * sum_train
          #
          # # cat(sum(lDensity_train*ltrain_fit)+sum(no_grouping_ll[-pos_trimmed_train]), "\n")
          # cat(sum(lDensity_train*ltrain_fit), "\n")
          # cat(sum(mat_train)+sum(no_grouping_ll[-pos_trimmed_train]),"selected variables",F_subset, "\n")
          # FIX this piece needs to be deleted: just for checking monot

          if (all(pos_old == pos_trimmed_train)) {
            criterion <- FALSE
          } else {
            pos_old <- pos_trimmed_train
            criterion <- (criterion) & (iter_init < max_iter_init)
          }
          #
          #           else {
          #             criterion <- FALSE
          #           }
        }

        mTau_train <-
          matrix(log(fitm$parameters$pro),
                 nrow(X_train_fit),
                 fitm$G,
                 byrow = TRUE)

        lDensity_train <-
          tryCatch(
            sapply(1:fitm$G, function(g)
              mvnfast::dmvn(
                X = as.matrix(X_train_fit[, F_subset, drop = FALSE]),
                mu = fitm$parameters$mean[F_subset, g],
                sigma = fitm$parameters$variance$sigma[F_subset, F_subset, g],
                log = T
              )),
            error = function(e) {
              return(list(ll = NA))
            }
          )
        sum_train <- mTau_train + lDensity_train
        mat_train <- ltrain_fit * sum_train
        ll <- sum(mat_train)

        # ll <- sum(unit_likelihood[-pos_trimmed_train])
        # cat(ll, "\n")
        # } else {
        #   ll <- NA
        # }

        return(
          list(
            ll = ll,
            fitm = fitm,
            relevant_variables = F_subset
            # irrelevant_variables = E_subset,
            # m_E_F = m_E_F,
            # V_E_F = V_E_F,
            # G_E_F = G_E_F
          )
        )
      }
    n_init <-
      replicate(
        nsamp,
        expr = robust_start(
          X_train=X_train,
          class_train=class_train,
          ltrain=ltrain,
          alpha_train=alpha_train,
          n_relevant_variables=n_relevant_variables,
          model_name=model_name,
          swap_step=swap_step,
          max_iter_init=max_iter_init,
          ctrl_GA=ctrl_GA
        ),
        simplify = TRUE
      )

    if (!all(is.na(n_init[1, ]))) {
      ind <- which.max(n_init[1,])
    } else {
      ind <-
        1 #if no initial subsample worked, get the first one and then the alg will break down further in the code
    }
    return(n_init[, ind])

  }
