
# Utility control functions -----------------------------------------------

control_init_EMST <- function(n_samp = 50,
                         n_start_extra_classes = 30,
                         max_iter = 1e02) {
  # Initialization control paramaters
  list(
    n_samp = n_samp,
    n_start_extra_classes = n_start_extra_classes,
    max_iter = max_iter
  )
}

# Function ctrlGA sets parameters of the genetic algorithm
# used for swap step.
# Arguments correspond to those of function kofnGA in the kofnGA package.

control_GA <- function (popsize = 200,
                        keepbest = floor(popsize / 10),
                        ngen = 500,
                        tourneysize = max(ceiling(popsize / 10), 2),
                        mutprob = 0.01,
                        mutfrac = NULL,
                        initpop = NULL,
                        cluster = NULL,
                        sharedmemory = FALSE
)
{
  list(popsize = popsize,
       keepbest = keepbest,
       ngen = ngen,
       tourneysize = tourneysize,
       mutprob = mutprob,
       mutfrac = mutfrac,
       initpop = initpop,
       cluster = cluster,
       sharedmemory = sharedmemory)
}
