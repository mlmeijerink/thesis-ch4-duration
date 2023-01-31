source("functions.R")

library(remify)
library(remstats)
library(prodlim)
library(parallel)

# Should we run the analysis on a small sample of the data?
sample <- FALSE

# Load the data ---------------------------------------------------------------
load("Data/hospital.RData")

# Add date and time information ----------------------------------------------
# Start date and time of the data collection
start <- as.POSIXct("2010-12-06 13:00 pm")
# Add the data and time of the events to the edgelist
edgelist$date <- start + edgelist$time

# Take a sample ---------------------------------------------------------------
if(sample) {
	# Random seed
	set.seed(742416)
	# Probabilities to be observed per status type
	probs <- data.frame(prop.table(table(attributes$status)))
	# Relate these probabilities to the actors
	probs <- probs$Freq[match(attributes$status, probs$Var1)]
	# Let the probabilities sum to one
	probs <- probs/sum(probs)
	# Sample 25 actors with these probabilities
	ids <- sample(x = attributes$id, size = 25, prob = probs)
	# Select these actors
	attributes <- subset(attributes, id %in% ids)
	# Select events with the sampled actors
	edgelist <- subset(edgelist, actor1 %in% ids & actor2 %in% ids)
}

# Interpolate events ---------------------------------------------------------
# Collect the observed dyads 
obs.dyads <- edgelist[,c(2,3)]
obs.dyads <- obs.dyads[!duplicated(obs.dyads),]
# Connect each event with the dyad id
edgelist$edge <- prodlim::row.match(edgelist[,c(2,3)], obs.dyads)
# Add when the events end
edgelist$time.end <- edgelist$time + edgelist$duration
# Add event numbers
edgelist$event.number <- 1:nrow(edgelist)
# Collect the events per observed dyad in a list
dyad.events <- lapply(1:nrow(obs.dyads), function(i) {
	subset(edgelist, edge == i)
})
# Per dyad, interpolate events within 75 seconds 
ip.dyad.events <- lapply(dyad.events, function(x) {
	# First observed event with the dyad
	m <- 1
	# Loop over all observed events with the dyad
	while(m < nrow(x)) {
		# Check if the next event starts within 75 seconds from the end of the 
		# current event
		if(x$time[m+1] <= x$time.end[m] + 75) {
			# Change the end time of the current event with the end time of the next
			# event
			x$time.end[m] <- x$time.end[m+1]
			# Remove the next event from the edgelist 
			x <- x[-(m+1),]
		}
		# Continue with the next event 
		m <- m + 1
	}
	# Update the event duration 
	x$duration <- x$time.end-x$time
	# Output
	x
})
# Bind the events per dyad in a new edgelist
ip.edgelist <- do.call(rbind, ip.dyad.events)
# Order the events based on their original event number
ip.edgelist <- ip.edgelist[order(ip.edgelist$event.number),]
# Replace the old edgelist
edgelist <- ip.edgelist

# Collect actors presence in the hospital (i.e., their shifts) ----------------
# Collect shift information per actor
shifts <- lapply(attributes$id, function(id) {
	# Collect the events with the respective actor
	actor.events <- subset(edgelist, actor1 == id | actor2 == id)
	# If the actor is a patient, let their stay correspond to the first 
	# event until the last event with this actor
	if(attributes$status[attributes$id == id] == "PAT") {
		x <- data.frame(start = 1, end = nrow(actor.events))
	} else {
		# Define a new shift if the time until the next event is 7 hours or longer
		cutpoints <- which(diff(actor.events$time) >= 7*60*60)
		x <- lapply(1:(length(cutpoints) + 1), function(i) {
			# Start of the shift
			start <- ifelse(i == 1, 1, cutpoints[i - 1]+1)
			# End of the shift
			end <- ifelse(i == length(cutpoints)+1, nrow(actor.events), cutpoints[i])
			# Collect in a data.frame
			data.frame(start = start, end = end)
		})
		# Bind shifts
		x <- data.frame(do.call(rbind, x))
	}
	# Add date and time
	x$start_date <- actor.events$date[x$start]
	x$end_date <- actor.events$date[x$end]
	# Add actor id and status
	x$id <- id
	x$status <- attributes$status[attributes$id == id]
	# Output
	x
})
# Bind rows
shifts <- do.call(rbind, shifts)
# Add length of the shift
shifts$hours <- as.numeric((shifts$end_date - shifts$start_date)/(60*60))

# Create the risk set ---------------------------------------------------------
# Collect the actors
actors <- attributes$id
# Dummy edgelist 
dummy <- data.frame(time = 1, actor1 = actors[1], actor2 = actors[2])
# Run remstats 
statsObject <- tomstats(edgelist = dummy, effects = ~ 1, actors = actors,
	directed = FALSE, origin = 0)
# Extract the risk set
riskset <- statsObject$riskset
# Change the column names
colnames(riskset)[c(1,2)] <- c("actor1", "actor2")
# Add time-variable
riskset$time <- 0
# Update edge
edgelist$edge <- prodlim::row.match(edgelist[,c(2,3)], riskset[,c(1,2)])

# Create the supplist ---------------------------------------------------------
# Empty supplist
supplist <- matrix(data = TRUE, nrow = nrow(edgelist), ncol = nrow(riskset))
# Loop over events
for(i in 1:nrow(edgelist)) {
	# Are events available for interaction based on their risk set time?
	if(i > 1) {supplist[i,] <- (riskset$time <= edgelist$time[i])}
	# Add whether events are available for interaction based on their shift:
	# Indicate the available actors
	aa <- subset(shifts, start_date <= edgelist$date[i] &
			end_date >= edgelist$date[i])$id
	# Indicate the available dyads
	ad <- subset(riskset, actor1 %in% aa & actor2 %in% aa)$id
	# Indicate whether dyads are available
	supplist[i,which(!(riskset$id %in% ad))] <- FALSE
	# Update the time variable
	riskset$time[riskset$actor1 == edgelist$actor1[i] & 
			riskset$actor2 == edgelist$actor2[i]] <- 
		edgelist$time[i] + edgelist$duration[i]
}

# Prepare tie effects ---------------------------------------------------------
# Both actors are administrative staff
bothADM <- both(attributes, "status", "ADM") 
# Both actors are nurses
bothNUR <- both(attributes, "status", "NUR") 
# Both actors are medical staff
bothMED <- both(attributes, "status", "MED") 
# Both actors are patients
bothPAT <- both(attributes, "status", "PAT") 
# One actor is nurse, the other administrative staff
combNurAdm <- combination(attributes, "status", "NUR", "ADM")
# One actor is nurse, the other medical staff
combNurMed <- combination(attributes, "status", "NUR", "MED")
# One actor is patient, the other administrative staff
combPatAdm <- combination(attributes, "status", "PAT", "ADM")
# One actor is patient, the other medical staff
combPatMed <- combination(attributes, "status", "PAT", "MED")
# One actor is medical staff, the other administrative staff
combMEDAdm <- combination(attributes, "status", "MED", "ADM")

# Create the pshift join ------------------------------------------------------
dummy <- remstats(edgelist = edgelist, attributes = attributes$id, 
	directed = FALSE, tie_effects = ~ 1)
psOngoing <- get_psOngoing(dummy)

# Specify the remstats effects ------------------------------------------------
form <- ~ 1 + 
	(totaldegreeDyad(scaling = "std") + 
			inertia(scaling = "std") + 
			sp(scaling = "std") + 
			userStat(psOngoing, "join"))*
	(tie(bothADM, "bothADM") +
			tie(bothNUR, "bothNUR") +
			tie(bothMED, "bothMED") +
			tie(bothPAT, "bothPAT") +
			tie(combNurAdm, "combNurAdm") +
			tie(combNurMed, "combNurMed") +
			tie(combPatAdm, "combPatAdm") +
			tie(combPatMed, "combPatMed") +
			tie(combMEDAdm, "combMEDAdm"))

# Model selection ------------------------------------------------------------
# psi candidates
psi.candidates <- seq(-1, 2, length.out = 10)
# Tau candidates
tau.candidates <- seq(5, 48*60, length.out = 5)*60
# Parallel computing
cl <- makeCluster(4)
clusterEvalQ(cl, library(remstats))
clusterExport(cl, c("form", "psOngoing", "bothADM", "bothNUR", "bothMED", 
	"bothPAT", "combNurAdm", "combNurMed", "combPatAdm", "combPatMed", 
	"combMEDAdm"))
# Select the model for the event rate
loglikR <- model.select.rate(psi.candidates, tau.candidates, edgelist, 
	effects = form, actors = attributes$id, directed = FALSE, 
	supplist = supplist, cl = cl)
# Parallel computing
cl <- makeCluster(4)
clusterEvalQ(cl, library(remstats))
clusterExport(cl, c("form", "psOngoing", "bothADM", "bothNUR", "bothMED", 
	"bothPAT", "combNurAdm", "combNurMed", "combPatAdm", "combPatMed", 
	"combMEDAdm"))
# Select the model for the event duration
loglikD <- model.select.duration(psi.candidates, tau.candidates, edgelist, 
	effects = form, actors = attributes$id, directed = FALSE, cl = cl)

# Model selection run 2 ------------------------------------------------------	
# Refine psi and tau 
psi.candidates2 <- seq(0, 0.6, length.out = 11)
tau.candidates2 <- c(2.5, 5, 15, 30, 60)*60
# Parallel computing
cl <- makeCluster(4)
clusterEvalQ(cl, library(remstats))
clusterExport(cl, c("form", "psOngoing", "bothADM", "bothNUR", "bothMED", 
	"bothPAT", "combNurAdm", "combNurMed", "combPatAdm", "combPatMed", 
	"combMEDAdm"))
# Event rate
loglikR2 <- model.select.rate(psi.candidates2, tau.candidates2, edgelist, 
	effects = form, actors = attributes$id, 
	directed = FALSE, supplist = supplist)
# Refine psi
psi.candidates3 <- seq(0.5, 1.5, length.out = 11)
# Parallel computing
cl <- makeCluster(4)
clusterEvalQ(cl, library(remstats))
clusterExport(cl, c("form", "psOngoing", "bothADM", "bothNUR", "bothMED", 
	"bothPAT", "combNurAdm", "combNurMed", "combPatAdm", "combPatMed", 
	"combMEDAdm"))
# Event duration
loglikD2 <- model.select.duration(psi.candidates3, tau.candidates2, edgelist, 
	effects = form, actors = attributes$id, 
	directed = FALSE)

# Save results
save.image(file = "Results/application1_modelselection.RData")

# Model estimation (event rate) -----------------------------------------------
# Set the parameters
psi.rate <- 0.42
tau.rate <- 150
# Set up parallel computing
cl <- makeCluster(4)
clusterEvalQ(cl, library(remstats))
clusterExport(cl, c("form", "psOngoing", "bothADM", "bothNUR", "bothMED", 
	"bothPAT", "combNurAdm", "combNurMed", "combPatAdm", "combPatMed", 
	"combMEDAdm", "tau.rate", "psi.rate", "edgelist", "attributes"))
# Compute the statistics
statsObject <- dremStats(effects = form, edgelist = edgelist, 
	actors = attributes$id, directed = FALSE, memory = "Brandes", 
	memory_value = tau.rate, psi = psi.rate)
# Fix event timing
rehObject <- reh(edgelist = edgelist, actors = attributes$id, 
	directed = FALSE, origin = 0, model = "tie")
statsObject$evls[,2] <- cumsum(rehObject$intereventTime)
# Parameter estimation 
fitR <- relevent::rem(eventlist = statsObject$evls, 
	statslist = statsObject$statistics, supplist = supplist,
	estimator = "MLE", timing = "interval")
# Compute lambda
lambda <- compute_lambda(statsObject$statistics, fitR$coef, supplist)
# Rank of the observed event (lower is better)
rank <- sapply(1:nrow(statsObject$evls), function(i) {
	nrow(statsObject$riskset) -
		rank(lambda[i,], ties.method = "min")[statsObject$evls[i,1]]
})

# Compute the statistics
statsObject.Null <- dremStats(effects = form, edgelist = edgelist, 
	actors = attributes$id, directed = FALSE)
# Fix event timing
rehObject.Null <- reh(edgelist = edgelist, actors = attributes$id, 
	directed = FALSE, origin = 0, model = "tie")
statsObject.Null$evls[,2] <- cumsum(rehObject.Null$intereventTime)
# Parameter estimation 
fitR.Null <- relevent::rem(eventlist = statsObject.Null$evls, 
	statslist = statsObject.Null$statistics, supplist = supplist,
	estimator = "MLE", timing = "interval")
# Compute lambda
lambda.Null <- compute_lambda(statsObject.Null$statistics, fitR.Null$coef, supplist)
# Rank of the observed event (lower is better)
rank.Null <- sapply(1:nrow(statsObject.Null$evls), function(i) {
	nrow(statsObject.Null$riskset) -
		rank(lambda.Null[i,], ties.method = "min")[statsObject.Null$evls[i,1]]
})

# Save results
save.image(file = "Results/application1_modelselection.RData")

# Model estimation (event duration) -----------------------------------------------
# Set the parameters
psi.duration <- 1
tau.duration <- 300
# Sequence statistics
cl <- makeCluster(4)
clusterEvalQ(cl, library(remstats))
clusterExport(cl, c("form", "psOngoing", "bothADM", "bothNUR", "bothMED", 
	"bothPAT", "combNurAdm", "combNurMed", "combPatAdm", "combPatMed", 
	"combMEDAdm", "tau.duration", "psi.duration", "edgelist", "attributes"))
statsObject <- dremStats(effects = form, edgelist = edgelist, 
	actors = attributes$id, directed = FALSE, memory = "Brandes", 
	memory_value = tau.duration, psi = psi.duration)
# Observed event indices
obsEvents <- statsObject$evls[,1]
# Select the observed event statistics
stats.D <- t(sapply(1:length(obsEvents), function(i){
	statsObject$statistics[i,obsEvents[i],]
}))
# Remove the intercept
stats.D <- stats.D[,-1]
# Transform to data.frame
stats.D <- data.frame(stats.D)
# Fit with survreg 
fitD <- survreg(Surv(edgelist$duration) ~ ., data = stats.D, 
	dist = "exponential")

# Save results
save.image(file = "Results/application1_modelselection.RData")