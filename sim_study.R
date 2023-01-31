# Code that executes the simulation study. 
rm(list = ls())

source("functions.R")
library(relevent)
library(survival)
library(future)

#' ============================================================================
#  1. User-defined settings ---------------------------------------------------
#' ============================================================================
# Set a seed 
set.seed(2386)

# Choose the number of replications 
R <- 100

# M 
M <- 5000

# N
N <- 20

# Choose a scenario 
scenario <- 1

# Halflife (tau) and duration weight (psi) parameter per scenario
if (scenario == 1) {
	# Scenario 1 
	psi <- c(rate = 0.6, duration = 0.6)
	tau <- c(rate = 300, duration = 300)
} else if (scenario == 2) {
	# Scenario 2 
	psi <- c(rate = 0.6, duration = 0.6)
	tau <- c(rate = 5000, duration = 5000)
} else if (scenario == 3) {
	# Scenario 3 
	psi <- c(rate = 0.6, duration = 0.6)
	tau <- c(rate = 7530, duration = 7530)
} else if (scenario == 4) {
	# Scenario 4 
	psi <- c(rate = 0.3, duration = 0.3)
	tau <- c(rate = 300, duration = 300)
} else if (scenario == 5) {
	# Scenario 5 
	psi <- c(rate = 0, duration = 0)
	tau <- c(rate = 300, duration = 300)
} else if (scenario == 6) {
	# Scenario 6
	psi <- c(rate = 0.3, duration = 0.6)
	tau <- c(rate = 300, duration = 5000)
} 

# Effects
effects_rate <- ~ 1 + outdegreeSender(scaling = "std") + 
	inertia(scaling = "std") + otp(scaling = "std")
effects_duration <- ~ 1 + indegreeReceiver(scaling = "std") + 
	inertia(scaling = "std") + otp(scaling = "std")

# Parameters
beta <- c(baseline = -8e0, outdegreeSender = 2e-01, inertia = 2e-01, 
	otp = 2e-01)
theta <- c(baseline = 1, indegreeReceiver = 2e-01, inertia = -2e-01, 
	otp = -2e-01)

#' ============================================================================
#  3. Generate the data ------------------------------------------------------
#' ============================================================================
# Generate using generateSequence
cl <- makeCluster(4)
clusterEvalQ(cl, c(library(remstats), library(progress)))
clusterExport(cl, c("generateSequence", "effects_rate", "effects_duration", 
	"beta", "theta", "tau", "psi", "N", "M"))

dat <- parLapply(cl, 1:R, function(r) {
	generateSequence(N = N, M = M, effects_rate = effects_rate, 
		effects_duration = effects_duration, 
		beta = beta, theta = theta, tau = tau, psi = psi, 
		verbose = FALSE)
})

# Save generated data
save.image(file = paste0("Data/sim_dat", scenario, ".RData"))

#' ============================================================================
#  4. Run the Model selection procedure ---------------------------------------
#' ============================================================================

# Set things up ---------------------------------------------------------------
# Candidate parameter values 
psiC <- seq(-0.15, 0.8, 0.15) 
tauC <- c(150, seq(300, 10000, 2350)) 

# Effects
effects <- ~ 1 + outdegreeSender(scaling = "std") + 
	indegreeReceiver(scaling = "std") + inertia(scaling = "std") + 
	otp(scaling = "std")
varRate <- c("baseline", "outdegreeSender", "inertia", "otp")

# Saving space
fitR <- loglikR <- list()
fitD <- loglikD <- list()

# Progress indicator
pb <- txtProgressBar(min = 0, max = R*length(psiC)*length(tauC), style = 3)
i <- 1

plan(multisession(workers = 4))
# (r) Loop over data sets ------------------------------------------------------
for (r in 1:R) {
	# Saving space
	fitR[[r]] <- fitD[[r]] <- list()
	loglikR[[r]] <- matrix(NA, nrow = length(psiC), ncol = length(tauC))
	loglikD[[r]] <- matrix(NA, nrow = length(psiC), ncol = length(tauC))
	# Select the generated edgelist
	edgelist <- dat[[r]]
	# Define the support set
	supplist <- atRisk_sim(eventseq = edgelist, actors = 1:N)
	
	# (j) Loop over psi ------------------------------------------------------
	for (j in 1:length(psiC)) {
		# Saving space
		fitR[[r]][[j]] <- list()
		fitD[[r]][[j]] <- list()
		
		# (k) Loop over tau ---------------------------------------------------
		for (k in 1:length(tauC)) {
			
			# Determine the type of memory process
			memory <- ifelse(tauC[k] == Inf, "full", "Brandes")
			
			# Compute statistics
			statsObject <- duration_remstats_sim(effects = effects, 
				edgelist = edgelist, 
				actors = 1:N, origin = 0, memory = memory, 
				memory_value = tauC[k], psi = psiC[j])
			
			### Event rate ### --------------------------------------------	
			# Estimate parameters
			fitR[[r]][[j]][[k]] <- relevent::rem(eventlist = statsObject$evls, 
				statslist = statsObject$statistics[,,varRate], 
				supplist = supplist, 
				estimator = "MLE", timing = "interval")
			# Get the log-likelihood 
			loglikR[[r]][j,k] <- fitR[[r]][[j]][[k]]$loglik
			
			### Event duration ### --------------------------------------------	
			# Select the covariates data.frame
			events <- statsObject$evls[,1]
			statsD <- t(sapply(1:length(events), function(i) {
				statsObject$statistics[i,events[i],]
			}))
			statsD <- statsD[,-1]
			statsD <- data.frame(statsD)
			
			# Fit with survreg
			fitD[[r]][[j]][[k]] <- survreg(Surv(edgelist$duration) ~ 1 + 
					indegreeReceiver + inertia + otp, data = statsD, 
				dist = "exponential")
			# Get the log-likelihood 
			loglikD[[r]][j,k] <- fitD[[r]][[j]][[k]]$loglik[2]
			
			# Progress update
			setTxtProgressBar(pb, i)
			i <- i + 1
		}
	}
}

# Save results
save.image(file = paste0("Results/sim_results", scenario, ".RData"))

#' ============================================================================
#  5. Inspect the results  ----------------------------------------------------
#' ============================================================================

# Event rate
indices <- sapply(loglikR, function(x) {
	which(x == max(x), arr.ind = T) 
})

psiEst <- sapply(indices[1,], function(j) {psiC[j]})
tauEst <- sapply(indices[2,], function(j) {tauC[j]})

mean(psiEst)
mean(tauEst) 

# Event duration
indices <- sapply(loglikD, function(x) {
	which(x == max(x), arr.ind = T) 
})

psiEst <- sapply(indices[1,], function(j) {psiC[j]})
tauEst <- sapply(indices[2,], function(j) {tauC[j]})

mean(psiEst)
mean(tauEst)