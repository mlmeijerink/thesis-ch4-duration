# Required packages 
library(Rcpp)
library(remify)
library(remstats)
library(relevent)
library(progress)
library(survival)
library(SurvRegCensCov)
library(RColorBrewer)
library(parallel)
library(future)
library(future.apply)

#' generateSequence ------------------------------------------------------------
#' 
#' See Algorithm 1. 
#' 
#' @param N number of actors
#' @param M number of events 
#' @param effects_rate remstats formula with effects for the event rate 
#' @param effects_duration remstats formula with effects for the event duration
#' @param beta vector with model parameters for the event rate 
#' @param theta vector with model parameters for the event duration
#' @param tau vector of length 2, the first refers to the halflife for the
#' event rate, the second to the halflife for the event duration 
#' @param psi vector of length 2, the first refers to the duration weight 
#' parameter for the event rate, the second to the duration weight parameter 
#' for the event duration 
#' @param random the amount of randomly sampled events to initialize the 
#' statistics (default is 0)

generateSequence <- function(N, M, effects_rate, effects_duration, beta, theta, 
	tau, psi, verbose = TRUE) {
	
	# Step 1: Prepare input values where necessary 
	beta <- matrix(beta)    # Make sure beta is a matrix
	theta <- matrix(theta)  # Make sure theta is a matrix
	
	tau_rate <- tau[1]      # Half-life for the event rate
	tau_duration <- tau[2]	# Half-life for the event duration
	
	psi_rate <- psi[1]			# Duration weight parameter for the event rate
	psi_duration <- psi[2]	# Duration weight parameter for the event duration
	
	# Step 2: Set t0
	t <- 0
	
	# Step 3: Initialize the risk set
	dummy <- data.frame(time = 1, actor1 = 1, actor2 = 1, 
		duration = 1, weight = 1)
	out <- tomstats(effects = ~ 1, dummy, actors = 1:N)
	riskset <- out$riskset
	# Add a time variable to the risk set that indicates when events end, i.e., 
	# when events are available to be sampled again
	riskset$time <- t
	
	# Step 4: Initialize the statistics 
	# Set endogenous statistics equal to zero
	stats_rate <- array(0, dim = c(1, N*(N-1), length(beta)))
	# Set the baseline statistic equal to one
	stats_rate[,,1] <- 1
	# Set endogenous statistics equal to zero      
	stats_duration <- array(0, dim = c(1, N*(N-1), length(theta)))
	# Set the baseline statistic equal to one
	stats_duration[,,1] <- 1
	
	# Check whether the effects, halflife (tau) and duration weight param (psi) 
	# for the event rate and event duration part of the model are the same 
	effCheck <- (effects_rate == effects_duration & tau_rate == tau_duration & 
			psi_rate == psi_duration)
	
	# Saving space
	edgelist <- data.frame() 
	
	# Progress
	pb <- progress_bar$new(total = M)
	
	## TEMP
	supplist <- list()
	
	# Step 5: Start the loop
	for(m in 1:M) {
		# Determine the dyads that are at risk for experiencing an event
		atrisk <- (riskset$time <= t)
		
		# Step 6: Compute the event rate
		lambda <- exp(stats_rate[1,,]%*%beta)
		
		# Step 7: Sample the waiting time 
		tDelta <- rexp(n = 1, rate = sum(lambda[atrisk]))
		
		# Step 8: Sample the dyad
		probs <- lambda/sum(lambda[atrisk])
		probs[!atrisk] <- 0
		dyad <- sample(x = 1:nrow(riskset), size = 1, prob = probs)
		
		# Step 9: Compute the duration rate
		eta <- exp(stats_duration[1,dyad,]%*%theta)
		
		# Step 10: Sample the duration 
		# Commented Weibull out (see Bender et al., 2003 for the procedure)
		# usample <- runif(n = 1)
		# dm <- (-log(usample)/lambdaD)^(1/gamma)
		dm <- rexp(n = 1, rate = eta)
		
		# Save the sampled values in the edgelist
		edgelist <- rbind(
			edgelist, 
			data.frame(
				time = t + tDelta,
				actor1 = riskset$sender[dyad],
				actor2 = riskset$receiver[dyad],
				duration = dm,
				weight = 1
			))
		
		# Step 11: Update the time
		t <- t + tDelta
		
		# Step 12: Update the statistics 
		# Determine the type of memory process
		memory <- ifelse(tau_rate == Inf, "full", "Brandes")
		# Temporary edgelist to compute statistics
		if (m == 1) {
			temp <- edgelist[m,]
		} else {
			# Select all previous events
			temp <- edgelist[1:(m-1),]
			# Change the time to when the events ended
			temp$time <- temp$time + temp$duration
			# Add the just sampled event
			temp <- rbind(temp, edgelist[m,])
			# Only select events that ended before the just sampled event started 
			temp <- subset(temp, time <= edgelist$time[m])
		}
		temp <- temp[order(temp$time),]
		# Determine the weight of the events
		temp$weight <- temp$duration^psi_rate
		# Compute the statistics for the rate
		out_rate <- tomstats(
			edgelist = temp, 
			effects = effects_rate, 
			actors = 1:N, 
			memory = memory,
			memory_value = tau_rate,
			start = nrow(temp), stop = nrow(temp),
			output = "stats_only")
		stats_rate <- out_rate$statistics
		# Update the statistics for the duration 
		if(effCheck) {
			stats_duration <- stats_rate
		} else {
			# Determine the type of memory process
			memory <- ifelse(tau_duration == Inf, "full", "Brandes")
			# Determine the weight of the events
			temp$weight <- temp$duration^psi_duration
			# Compute the statistics for the duration
			out_duration <- tomstats(
				edgelist = temp, 
				effects = effects_duration, 
				actors = 1:N,
				memory = memory,
				memory_value = tau_duration,
				start = nrow(temp), stop = nrow(temp),
				output = "stats_only")
			stats_duration <- out_duration$statistics
		}
		
		## TEMP
		supplist[[m]] <- atrisk
		
		# Step 13: Update the risk set
		riskset$time[
			(riskset$sender == riskset$sender[dyad] | 
					riskset$receiver == riskset$sender[dyad] | 
					riskset$sender == riskset$receiver[dyad] | 
					riskset$receiver == riskset$receiver[dyad]) & 
				riskset$time < as.numeric(t + dm)] <- t + dm
		
		# Progress indicator
		if (verbose) pb$tick()
	}
	
	# Output
	return(edgelist[,-5])
	
}


#' atRisk
#' 
#' Generates a "supplist" object for relevent::rem, i.e., a matrix with on the 
#' rows the time points and in the columns the pairs, where a FALSE indicates 
#' that a pair is not available for interaction at the respective time because 
#' at least one of the actors is already involved in an ongoing interaction. 
#' 
#' @param eventseq sequence with all the events (time, actor1, actor2, 
#' duration).
#' @param actors vector with all actors that should be included in the risk set.
#' @param directed whether events in the risk set are directed. Can also be 
#' set as an attribute of eventseq.
#' 
#' @return matrix 
#' @export

atRisk <- function(eventseq, actors, directed = TRUE) {
	# Set missing event duration equal to 0
	eventseq$duration[is.na(eventseq$duration)] <- 0
	
	# Create the risk set
	out <- tomstats(effects = ~ 1, edgelist = eventseq, actors = actors, 
		directed = directed)
	riskset <- out$riskset
	names(riskset)[c(1,2)] <- c("sender", "receiver")
	
	# Add a time variable to the risk set that indicates when events end after 
	# they have started, i.e., when events are available for interaction again
	riskset$time <- 0
	
	# Saving space
	out <- matrix(TRUE, nrow(eventseq), nrow(riskset))
	
	# Loop over events
	for(i in 1:nrow(eventseq)) {
		# Are events available for interaction?
		if(i > 1) {out[i,] <- (riskset$time <= eventseq$time[i])}
		
		# Update the time variable
		riskset$time[
			(riskset$sender == eventseq$actor1[i] | 
					riskset$receiver == eventseq$actor1[i] | 
					riskset$sender == eventseq$actor2[i] | 
					riskset$receiver == eventseq$actor2[i]) & 
				riskset$time < as.numeric(eventseq$time[i] + eventseq$duration[i])] <- 
			eventseq$time[i] + eventseq$duration[i]
	}
	
	out
}

atRisk_sim <- function(eventseq, actors, directed = TRUE) {
	# Set missing event duration equal to 0
	eventseq$duration[is.na(eventseq$duration)] <- 0
	
	# Create the risk set
	out <- tomstats(effects = ~ 1, edgelist = eventseq, actors = actors, 
		directed = directed)
	riskset <- out$riskset
	names(riskset)[c(1,2)] <- c("sender", "receiver")
	
	# Add a time variable to the risk set that indicates when events end after 
	# they have started, i.e., when events are available for interaction again
	riskset$time <- 0
	
	# Saving space
	out <- matrix(TRUE, nrow(eventseq), nrow(riskset))
	
	# Loop over events
	for(i in 1:nrow(eventseq)) {
		# Are events available for interaction?
		if(i > 1) {out[i,] <- (riskset$time <= eventseq$time[i])}
		# Remove the actors in the last event (only in the simulation study)
		if(i > 1) {
			ind <- which(riskset$sender == eventseq$actor1[i-1] | 
					riskset$receiver == eventseq$actor1[i-1] | 
					riskset$sender == eventseq$actor2[i-1] | 
					riskset$receiver == eventseq$actor2[i-1])
			out[i,ind] <- FALSE 
		}
		
		# Update the time variable
		riskset$time[
			(riskset$sender == eventseq$actor1[i] | 
					riskset$receiver == eventseq$actor1[i] | 
					riskset$sender == eventseq$actor2[i] | 
					riskset$receiver == eventseq$actor2[i]) & 
				riskset$time < as.numeric(eventseq$time[i] + eventseq$duration[i])] <- 
			eventseq$time[i] + eventseq$duration[i]
	}
	
	out
}

#' compute_lambda
#' 
#' Cpp function to compute lambda

cppFunction("
	
	arma::mat compute_lambda(arma::cube statistics, arma::colvec coef, 
	arma::mat supplist) {
	
		arma::mat lambda(statistics.n_rows, statistics.n_cols, arma::fill::zeros);
		
		for(arma::uword i = 0; i < coef.n_elem; ++i) {
			arma::mat statsslice = statistics.slice(i);
			arma::mat lambdaslice = coef(i)*statsslice;
			lambda += lambdaslice;
		}
		
		lambda = exp(lambda);
		lambda = lambda%supplist;
		
		return lambda;
	
	}", depends = "RcppArmadillo")


#' Model selection function  for the event rate.  
model.select.rate <- function(psi.candidates, tau.candidates, edgelist, 
	effects, actors = NULL, directed = TRUE, origin = NULL, supplist = NULL, 
	cl = cl) {
	# Set things up -------------------------------------------------------------
	if(any(edgelist$duration <= 0, na.rm = TRUE)) {
		stop("Duration values must be larger than zero.")
	}
	if(anyNA(edgelist$duration)) {
		warning("Missing event durations: these are in the computation of the statistics replaced with the minimum event duration")
	}
	
	# Set temporary event duration equal to the minimum event duration
	missing <- is.na(edgelist$duration)
	edgelist$duration[missing] <- min(edgelist$duration, na.rm = T)
	
	# Saving space for the log-likelihood
	loglik.R <- matrix(data = NA, nrow = length(psi.candidates), 
		ncol = length(tau.candidates))
	
	# Multiply time with 1000 (otherwise remify cannot always create 
	# the inter-event time object correctly)
	edgelist$time <- edgelist$time * 1000
	# Run remify::reh to prepare the edgelist
	rehObject <- reh(edgelist = edgelist, actors = actors, 
		directed = directed, origin = origin, model = "tie")
	edgelist$time <- edgelist$time / 1000
	
	# Progress bar
	pb <- txtProgressBar(min = 0, max = length(psi.candidates), 
		style = 3)
	
	# Loop over psi ------------------------------------------------------------
	for(j in 1:length(psi.candidates)) {
		# Loop over tau 
		for(k in 1:length(tau.candidates)) {
			# Determine the type of memory process
			memory <- ifelse(tau.candidates[k] == Inf, "full", "Brandes")
			# Sequence statistics
			statsObject <- dremStats(effects = effects, edgelist = edgelist, 
				actors = actors, directed = directed, origin = origin, 
				memory = memory, memory_value = tau.candidates[k], 
				psi = psi.candidates[j])
			
			# Fix event timing
			statsObject$evls[,2] <- cumsum(rehObject$intereventTime) / 1000
			
			# Prepare the supplist
			if(is.null(supplist)) {
				supplist <- matrix(data = 1, 
					nrow = nrow(edgelist), 
					ncol = ncol(statsObject$statistics))
			}
			# Parameter estimation 
			fit <- relevent::rem(eventlist = statsObject$evls, 
				statslist = statsObject$statistics, supplist = supplist, 
				estimator = "MLE", timing = "interval")
			# Log-likelihood 
			loglik.R[j,k] <- fit$loglik
		}
		# Update progress
		setTxtProgressBar(pb, j)
	}
	
	# Output
	loglik.R
}


#' model.select.duration
#' 
#' Model selection function  for the event duration  
model.select.duration <- function(psi.candidates, tau.candidates,
	edgelist, effects, actors = NULL, directed = TRUE, origin = NULL, cl = cl) {
	# Set things up --------------------------------------------------------------
	if(any(edgelist$duration <= 0, na.rm = T)) {
		stop("Duration values must be larger than zero.")
	}
	if(anyNA(edgelist$duration)) {
		warning("Missing event durations: these are in the computation of the statistics replaced with the minimum event duration and removed for the estimation of the model parameters")
	}
	
	# Set temporary event duration equal to the minimum event duration
	missing <- is.na(edgelist$duration)
	edgelist$duration[missing] <- min(edgelist$duration, na.rm = T)
	
	# Saving space for the log-likelihood 
	loglik.D <- matrix(data = NA, nrow = length(psi.candidates),
		ncol = length(tau.candidates))
	
	# Progress bar
	pb <- txtProgressBar(min = 0, max = length(psi.candidates), 
		style = 3)
	
	# Loop over psi ------------------------------------------------------------
	for(j in 1:length(psi.candidates)) {
		# Loop over tau 
		for(k in 1:length(tau.candidates)) {
			# Determine the type of memory process
			memory <- ifelse(tau.candidates[k] == Inf, "full", "Brandes")
			# Sequence statistics
			statsObject <- dremStats(effects = effects, edgelist = edgelist, 
				actors = actors, directed = directed, origin = origin, 
				memory = memory, memory_value = tau.candidates[k], 
				psi = psi.candidates[j])
			# Observed event indices
			events <- statsObject$evls[,1]
			# Remove events with missing event durations
			if(any(missing)) {
				events <- events[!missing]
				duration <- edgelist$duration[!missing]
			} else {
				duration <- edgelist$duration
			}
			# Select the observed event statistics
			stats.D <- t(sapply(1:length(events), function(i){
				statsObject$statistics[i,events[i],]
			}))
			# Remove the intercept
			stats.D <- stats.D[,-1]
			# Transform to data.frame
			stats.D <- data.frame(stats.D)
			# Fit with survreg 
			fit <- survreg(Surv(duration) ~ ., data = stats.D, 
				dist = "exponential")
			# Log-likelihood
			loglik.D[j,k] <- fit$loglik[2]
		}
		# Update progress
		setTxtProgressBar(pb, j)
	}
	
	# Output
	loglik.D
}

#' inspect.loglik.results

inspect.loglik.results <- function(loglik, psi.candidates, tau.candidates,
	max.tau = NULL) {
	# Find the best loglik
	best <- which(loglik == max(loglik), arr.ind = T)
	best.psi <- psi.candidates[best[1]]
	best.tau <- tau.candidates[best[2]]
	
	# Plotting with tau.candidates that include "inf"
	tau.candidates2 <- tau.candidates
	if(any(tau.candidates2 == Inf)) {
		tau.candidates2[tau.candidates2==Inf] <- max.tau
	}
	
	# (1) Filled contour plot
	#filled.contour(
	#	x = psi.candidates, 
	#	y = tau.candidates2, 
	#	z = loglik,
	#	xlab = expression(psi),
	#	ylab = expression(tau),
	#	cex.lab = 1,
	#	plot.axes = {
	#		axis(1)
	#		axis(2)
	#		contour(psi.candidates, tau.candidates2, loglik, add = TRUE)
	#	}
	#)
	
	# Prepare the next plots 
	# Set parameters
	par(mfrow = c(1,2))
	#par(ask=TRUE)
	# Get colors
	colors <- rep(brewer.pal(n = 12, name = "Paired"), 5)
	
	# (2) Plot one line for each psi.candidate
	plot(x = tau.candidates2, y = loglik[best[1],], type = "b", 
		ylim = c(min(loglik), (max(loglik) + 0.4*(max(loglik)-min(loglik)))),
		xlab = expression(tau), ylab = "Log likelihood",
		lwd = 2)
	for(i in 1:length(psi.candidates)) {
		if(i == best[1]) {next}
		points(tau.candidates2, loglik[i,], col = colors[i])
		lines(tau.candidates2, loglik[i,], col = colors[i])
	}
	col.legend <- colors[1:length(psi.candidates)]
	col.legend[best[1]] <- "black"
		lwd.legend <- rep(1, length(psi.candidates))
		lwd.legend[best[1]] <- 2
		legend("topleft", title = expression(psi),
			legend = round(psi.candidates, 2),
			col = col.legend, lwd = lwd.legend,
			lty = 1, ncol = 4, cex = 0.7)
		
		par(ask=FALSE)
		# (3) Plot one line for each tau.candidate
		plot(x = psi.candidates, y = loglik[,best[2]], type = "b", 
			ylim = c(min(loglik), (max(loglik) + 0.4*(max(loglik)-min(loglik)))),
			xlab = expression(psi), ylab = "Log likelihood", lwd = 2)
		for(i in 1:length(tau.candidates)) {
			if(i == best[2]) {next}
			points(psi.candidates, loglik[,i], col = colors[i])
			lines(psi.candidates, loglik[,i], col = colors[i])
		}
		col.legend <- colors[1:length(tau.candidates)]
		col.legend[best[2]] <- "black"
			lwd.legend <- rep(1, length(tau.candidates))
			lwd.legend[best[2]] <- 2
			legend("topleft", title = expression(tau),
				legend = round(tau.candidates, 2),
				col = col.legend, lwd = lwd.legend,
				lty = 1, ncol = 3, cex = 0.7)
			par(mfrow = c(1,1))
			
			# (4) Print results 
			loglik.table <- 
				array(loglik, dim = c(length(psi.candidates), length(tau.candidates)),
					dimnames = list(round(psi.candidates, 2), round(tau.candidates, 2)))
			loglik.table <- as.table(loglik.table)
			names(attributes(loglik.table)$dimnames) <- c("psi", "tau")
			print(loglik.table)
			
			cat(paste0("Best: psi is ", round(best.psi, 2), ", tau is ", 
				round(best.tau, 2)), "\r")
}


#' both
#' 
#' Function to create a matrix for a remstats "tie"-effect that indicates
#' whether two actors both have some attribute value. 

both <- function(attributes, variable, value) { 
	# Collect the actors
	actors <- sort(attributes$id)
	# Select the relevant variable
	var <- attributes[,variable]
	# Prepare an empty matrix to hold the information
	x <- matrix(data = NA, nrow = length(actors), ncol = length(actors))
	# Loop over actors
	for(i in 1:length(actors)) {
		for(j in 1:length(actors)) {
			# The stat is 1 if both actors have the value and 0 if not
			x[i,j] <- ifelse(var[attributes$id == actors[i]] == value &
					var[attributes$id == actors[j]] == value, 1, 0)
		}
	}
	# Update the row and column names
	rownames(x) <- colnames(x) <- actors
	# Output
	x
}


#' combination
#' 
#' Function to create a matrix for a remstats "tie"-effect that indicates
#' whether two actors score a combination of two attribute values. 

combination <- function(attributes, variable, value1, value2) {
	# Collect the actors
	actors <- sort(attributes$id)
	# Select the relevant variable
	var <- attributes[,variable]
	# Prepare an empty matrix to hold the information
	x <- matrix(data = NA, nrow = length(actors), ncol = length(actors))
	# Loop over actors
	for(i in 1:length(actors)) {
		for(j in 1:length(actors)) {
			# The stat is 1 if the actors score the combination of values on the
			# attribute
			x[i,j] <- ifelse(((var[attributes$id == actors[i]] == value1 &
					var[attributes$id == actors[j]] == value2) |
					(var[attributes$id == actors[i]] == value2 &
							var[attributes$id == actors[j]] == value1)), 1, 0)
		}
	}
	# Update the row and column names
	rownames(x) <- colnames(x) <- actors
	# Output
	x 
}

#' dremStats
#' 
#' A wrapper around remstats::tomstats to compute statistics for a model with 
#' event duration. 
#' @inheritParams remstats::tomstats
#' @param edgelist an object of class "data.frame" or "matrix" characterizing 
#' the relational event history sorted by time with columns "time", "actor1", 
#' "actor2", "duration". 
#' @param psi duration weight parameter. 
dremStats <- function(effects, edgelist, attributes = NULL, actors = NULL, 
	directed = TRUE, origin = NULL, memory = c("full", "window", "Brandes"), 
	memory_value = Inf, psi = 0) {
	
	# Dummy run to get objects
	statsOut <- tomstats(effects = ~ 1, edgelist = edgelist, actors = actors, 
		directed = directed, origin = origin, memory = "full")
	
	# Prepare edgelist
	edgelist$time_end <- edgelist$time + edgelist$duration
	edgelist$weight <- edgelist$duration^psi
	
	# Compute statistics in a loop: at each observed time point, compute based 
	# on the time ends of the previous observed events
	stats <- parLapply(cl, 1:nrow(edgelist), function(i) {
		if (i == 1) {
			altEdges <- edgelist[i,]
		} else {
			altEdges <- edgelist[1:(i-1),]
			altEdges$time <- altEdges$time_end
			altEdges <- rbind(altEdges, edgelist[i,])
			altEdges <- subset(altEdges, time <= edgelist$time[i])
		}
		
		altEdges <- altEdges[order(altEdges$time),]
		
		tomstats(effects = effects, edgelist = altEdges, 		
			attributes = attributes, actors = actors, directed = directed, 
			origin = origin, memory = memory, memory_value = memory_value, 
			start = nrow(altEdges), stop = nrow(altEdges), 
			output = "stats_only")$statistics
	})
	
	statsOut$statistics <- abind::abind(stats, along = 1)
	dimnames(statsOut$statistics)[[3]] <- dimnames(stats[[1]])[[3]]
	
	# Output
	statsOut	
}

duration_remstats_sim <- function(effects, edgelist, attributes = NULL, actors = NULL, 
	directed = TRUE, origin = NULL, memory = c("full", "window", "Brandes"), 
	memory_value = Inf, psi = 0) {
	
	# Dummy run to get objects
	statsOut <- tomstats(effects = ~ 1, edgelist = edgelist, actors = actors, 
		directed = directed, origin = origin, memory = "full")
	
	# Prepare edgelist
	edgelist$time_end <- edgelist$time + edgelist$duration
	edgelist$weight <- edgelist$duration^psi
	
	# Compute statistics in a loop: at each observed time point, compute based 
	# on the time ends of the previous observed events
	plan(multisession)
	stats <- future_lapply(1:nrow(edgelist), function(i) {
		if (i == 1 | i == 2) {
			altEdges <- edgelist[i,]
		} else {
			# Never regard the last event (only in the simulation)
			altEdges <- edgelist[1:(i-2),] 
			# Change the time of previous events to when they ended
			altEdges$time <- altEdges$time_end 
			# Add the previous event
			altEdges <- rbind(altEdges, edgelist[i-1,]) 
			# Only regard events that ended before the last starts
			altEdges <- subset(altEdges, time <= edgelist$time[i-1]) 
		}
		
		altEdges <- altEdges[order(altEdges$time),]
		
		tomstats(effects = effects, edgelist = altEdges, 		
			attributes = attributes, actors = actors, directed = directed, 
			origin = origin, memory = memory, memory_value = memory_value, 
			start = nrow(altEdges), stop = nrow(altEdges), 
			output = "stats_only")$statistics
	}, future.seed = TRUE)
	
	statsOut$statistics <- abind::abind(stats, along = 1)
	dimnames(statsOut$statistics)[[3]] <- dimnames(stats[[1]])[[3]]
	
	# Output
	statsOut	
}

#' get_psOngoing 
#' 
#' Compute ongoing event p-shift statistic. This statistic is equal to one at 
#' time t for all dyads with at least one actor currently involved in another 
#' event. 
#'
#' @param statsObject remstats output object
#' 
get_psOngoing <- function(statsObject) {
	
	# Elements of the statsObject
	evls <- statsObject$evls 
	riskset <- statsObject$riskset 
	edgelist <- statsObject$edgelist
	
	# Add a time variable to the risk set that indicates when events with these # actors end 
	riskset$time <- 0
	
	# Saving space
	stat <- matrix(NA, nrow = nrow(evls), ncol = nrow(riskset))
	
	# Loop over events to fill the statistic
	for(m in 1:nrow(evls)) {
		# Do dyads have actors currently in an interaction?
		stat[m,] <- ifelse(riskset$time <= evls[m,2], 0, 1)
		
		# Event attributes
		actor1 <- riskset$actor1[evls[m,1]]
		actor2 <- riskset$actor2[evls[m,1]]
		time <- edgelist$time[m]
		duration <- edgelist$duration[m]
		
		# Update the time variable (NOTE: to use with directed events, the 
		# column names of the risk set are different. Fix later.)
		riskset$time[(riskset$actor1 == actor1 | riskset$actor2 == actor1 | 
				riskset$actor1 == actor2 | riskset$actor2 == actor2) & 
				riskset$time < as.numeric(time + duration)] <- time + duration
	}    
	
	# Output
	return(stat)
}

