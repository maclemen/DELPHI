# Simulation program for novel design
# Clement Ma
# Febrary 24, 2018
# Modified March 4, 2018
# Modified July 3, 2018
# Modified August 31, 2018 - Require full cohort of X patients to complete before proceeding to next cohort.
# Modified January 28, 2019
# Modified May 15, 2019
# Modified September 21, 2020 - for DELPHI app
# Modified November 3, 2020 - for DELPHI app


library(dfcrm)

#####################
# TARGET.CRM FUNCTION: requires cohort size=3 regardless of intra-patient de-escalation
# Jan 28, 2019: include indicator "target.crm" 
	# target.crm=0: NO enrollment of patients at one dose below
	# target.crm=1: Enrollment of patients at one dose below
	# target.crm=2: Enrollment of patients at current best dose based on available information, cannot be higher than current dose

# May 15, 2019: stratified design "min.cohortB" - require minimum number of cohort B patients

target.crm <- function(prior, target.tox, number.trials, true.tox, arrival.rate, prop.B, target.crm, 
                       min.cohortB=0, cycle.length, cohort.size, max.N, start.level) {
# information of interest
total.patients <- rep(0, number.trials)
num.cohortA.patients <- rep(0, number.trials)
num.cohortB.patients <- rep(0, number.trials)
num.group1.patients <- rep(0, number.trials)
num.group2.patients <- rep(0, number.trials)
MTD.selection <- rep(0,number.trials)
study.duration <- rep(0,number.trials)
result.num.dose.changes <- rep(0,number.trials)

observe.tox <- mat.or.vec(nr=length(true.tox), nc=number.trials)
patient.allocation <- mat.or.vec(nr=length(true.tox), nc=number.trials)

start.time <- proc.time()
for (i in 1:number.trials) {
	#print (c("Trial ", i))
	# Enroll first patient
	num.slots <- cohort.size # counter to be in the main cohort
	num.waiting <- 0 # counter for number of cohort A patients waiting
	wait.cohort <- numeric(0) # array for the cohort of patients in wait list
	study.end <- 0 # flag to end study
	current.dose <- start.level # current dose level
	timeline.time <- 0 # overall study timeline
	num.dose.changes <- 0
	PID <- 1
	#print(PID)
	arrive.time <- timeline.time + rpois(1, arrival.rate)
	timeline.time <- arrive.time
	cohortB <- rbinom(1, 1, prob=prop.B)
	d <- current.dose
	# Determine toxicities
	dose.tox <- rbinom(1,1,prob=true.tox[d])
	if(dose.tox==1) {
		end.time <- timeline.time + runif(n=1, min=0, max=cycle.length)
	} else {
		end.time <- timeline.time + cycle.length
	}
	
	include <- 0
	group <- 1
	wait.list <- 0
	prev.cohort.end.time <- 0
	num.slots <- num.slots - 1 # Used up one slot
	pt <- data.frame (PID, arrive.time, cohortB, d, dose.tox, end.time, group, wait.list, include)
	
	# Update counts
	observe.tox[d,i] <- observe.tox[d,i]+dose.tox
	patient.allocation[d,i] <- patient.allocation[d,i]+1
	
	while (study.end == 0 | length(pt$PID[pt$cohortB==1]) < min.cohortB) { # NEW require at least min.cohortB patients
		# Patient arrives
		#print(PID)
		#print(pt)
		#print(timeline.time)

		if (study.end == 0) {
			PID <- PID+1
			arrive.time <- timeline.time + rpois(1, arrival.rate)
			timeline.time <- arrive.time
			cohortB <- rbinom(1, 1, prob=prop.B)
		} else if (study.end == 1) { # reached study end but not sufficient cohort B patients
			# print("study end == 1 but not sufficient cohort B patients")
			# sample patients until a cohort B patient arrives
			repeat {
				#print(PID)
				PID <- PID+1
				arrive.time <- timeline.time + rpois(1, arrival.rate)
				timeline.time <- arrive.time
				cohortB <- rbinom(1, 1, prob=prop.B)
				if (cohortB == 1) {
					break
				}
			}
		}			

		# Update which patients to include
		pt$include <- ifelse(pt$end.time <= timeline.time, 1, 0) # Modified 8/31/2018 to "<="
	
		# Check if new cohort can be enrolled [only include main cohort timeline]
		if (num.slots == 0 & timeline.time > tail(pt$end.time[pt$group==1],n=1)) {
			num.slots <- cohort.size
		}
	
		if (num.slots == cohort.size) { # start of main cohort
			wait.list <- 0
			# Check for waiting list patients
			if (num.waiting > 0) {
				#print ("Enrolling waitlist patient")
				arrive.time <- prev.cohort.end.time
				timeline.time <- arrive.time
				cohortB <- wait.cohort[1] # take first wait list patient
				wait.cohort <- wait.cohort[-1] # remove from waitlist
				# Update which patients to include
				pt$include <- ifelse(pt$end.time <= timeline.time, 1, 0) # Modified 8/31/2018 to "<="
				wait.list <- 1
				num.waiting <- num.waiting - 1
			}
	
			if (length(pt$PID[pt$include==1])==0) { # no completed observations yet
				d <- current.dose
			} else {
				run.crm <- crm(prior=prior, target=target.tox, tox=pt$dose.tox, level=pt$d, include=which(pt$include==1))
				#print(c("CRM dose", run.crm$mtd))
				if (run.crm$mtd > current.dose & current.dose < length(prior)) {
					current.dose <- current.dose + 1 # Can only escalate one dose higher
					d <- current.dose
					num.dose.changes <- num.dose.changes+1
				} else {
					if (run.crm$mtd < current.dose) {
						num.dose.changes <- num.dose.changes+1
					}
					current.dose <- run.crm$mtd
					d <- current.dose
					##num.slots <- cohort.size # NEW: if de-escalate, then reset num.slots to cohort.size
				}
				#print(c("Recommended dose", current.dose))
			}
			# Determine toxicities
			dose.tox <- rbinom(1,1,prob=true.tox[d])
			if(dose.tox==1) {
				end.time <- timeline.time + runif(n=1, min=0, max=cycle.length)
			} else {
				end.time <- timeline.time + cycle.length
			}

			group <- 1
			num.slots <- num.slots - 1 # Used up one slot
	
			# Add to pt tracker matrix
			pt <- rbind(pt, c(PID, arrive.time, cohortB, d, dose.tox, end.time, group, wait.list, include))
			
			# Update counts
			observe.tox[d,i] <- observe.tox[d,i]+dose.tox
			patient.allocation[d,i] <- patient.allocation[d,i]+1
	
		} else if (num.slots > 0 & num.slots < cohort.size) { # enroll in main cohort
			wait.list <- 0
			# Check for waiting list patients
			if (num.waiting > 0) {
				#print ("Enrolling waitlist patient")
				arrive.time <- prev.cohort.end.time
				timeline.time <- arrive.time
				cohortB <- wait.cohort[1] # take first wait list patient
				wait.cohort <- wait.cohort[-1] # remove from waitlist
				# Update which patients to include
				pt$include <- ifelse(pt$end.time <= timeline.time, 1, 0) # Modified 8/31/2018 to "<="
				wait.list <- 1
				num.waiting <- num.waiting - 1
			}

			if (length(pt$PID[pt$include==1])==0) { # no completed observations yet
				d <- current.dose
			} else {
				# Calculate recommended dose
				run.crm <- crm(prior=prior, target=target.tox, tox=pt$dose.tox, level=pt$d, include=which(pt$include==1))
				if (run.crm$mtd < current.dose) { # Within a cohort, can only stay at same dose or de-escalate
					current.dose <- run.crm$mtd
					#print(c("CRM dose", current.dose))
					d <- current.dose
					num.dose.changes <- num.dose.changes+1
					##num.slots <- cohort.size # NEW: if de-escalate, then reset num.slots to cohort.size
				} else {
					d <- current.dose
				}
				#print(c("Recommended dose", current.dose))
			}
			# Determine toxicities
			dose.tox <- rbinom(1,1,prob=true.tox[d])
			if(dose.tox==1) {
				end.time <- timeline.time + runif(n=1, min=0, max=cycle.length)
			} else {
				end.time <- timeline.time + cycle.length
			}
	
			group <- 1
			num.slots <- num.slots - 1 # Used up one slot
	
			# Add to pt tracker matrix
			pt <- rbind(pt, c(PID, arrive.time, cohortB, d, dose.tox, end.time, group, wait.list, include))
	
			# Update counts
			observe.tox[d,i] <- observe.tox[d,i]+dose.tox
			patient.allocation[d,i] <- patient.allocation[d,i]+1
	
		} else if (num.slots == 0) { # cohort B patients enrolled on one dose below current
			# Capture previous cohort end time
			prev.cohort.end.time <- tail(pt$end.time[pt$group==1],n=1)
	
			 	# YES capture Cohort B patients at one dose below
				if (cohortB == 1 & current.dose > 1 & target.crm==1) {
					d <- current.dose-1

					# Determine toxicities
					dose.tox <- rbinom(1,1,prob=true.tox[d])
					if(dose.tox==1) {
						end.time <- timeline.time + runif(n=1, min=0, max=cycle.length)
					} else {
						end.time <- timeline.time + cycle.length
					}
					group <- 2
					wait.list <- 0
					# Add to pt tracker matrix
					pt <- rbind(pt, c(PID, arrive.time, cohortB, d, dose.tox, end.time, group, wait.list, include))
		
					# Update counts
					observe.tox[d,i] <- observe.tox[d,i]+dose.tox
					patient.allocation[d,i] <- patient.allocation[d,i]+1
			 	# YES capture Cohort B patients at current recommended dose
				} else if (cohortB == 1 & target.crm==2) { # NEW: allow capture of Cohort B patients even at lowest dose level
					d <- current.dose

					# Determine toxicities
					dose.tox <- rbinom(1,1,prob=true.tox[d])
					if(dose.tox==1) {
						end.time <- timeline.time + runif(n=1, min=0, max=cycle.length)
					} else {
						end.time <- timeline.time + cycle.length
					}
					group <- 2
					wait.list <- 0
					# Add to pt tracker matrix
					pt <- rbind(pt, c(PID, arrive.time, cohortB, d, dose.tox, end.time, group, wait.list, include))
		
					# Update counts
					observe.tox[d,i] <- observe.tox[d,i]+dose.tox
					patient.allocation[d,i] <- patient.allocation[d,i]+1

				} else { # patients added to waiting list with 50% chance of actual enrollment
					if(rbinom(1,1,0.5)==1) {
						num.waiting <- num.waiting + 1
						wait.cohort <- cbind(wait.cohort, cohortB) # Capture which cohort it belongs to
						#print (c("Adding patient to waitlist: ", PID))
						#print (c("waitlist cohort: ", wait.cohort))
					}
				}

		} else {
			print (c("Error, number of slots is", numslots))
			break
		}

		# Check for study end
		if (length(pt$PID)==max.N) { 
			study.end<-1 
			# Study end, finish observing remaining patients
			timeline.time <- max(pt$end.time)
			pt$include <- 1
			run.crm.final <- crm(prior=prior, target=target.tox, tox=pt$dose.tox, level=pt$d, include=which(pt$include==1))
			#print(pt)
			#print (c("MTD is", run.crm.final$mtd))
			MTD.selection[i] <- run.crm.final$mtd
			total.patients[i] <- length(pt$PID)
			num.cohortA.patients[i] <- length(pt$PID[pt$cohortB==0])
			num.cohortB.patients[i] <- length(pt$PID[pt$cohortB==1])
			num.group1.patients[i] <- length(pt$PID[pt$group==1])
			num.group2.patients[i] <- length(pt$PID[pt$group==2])
			study.duration[i] <- timeline.time
			result.num.dose.changes[i] <- num.dose.changes
		}
	}
}

MTD.selection.table <- table(MTD.selection)
true.MTD <- which.min(round(abs(target.tox-true.tox),10))
PCS <- MTD.selection.table[true.MTD] / sum(MTD.selection.table)
obs.tox.overall <- sum(observe.tox)/sum(patient.allocation)

patient.allocation.table <- rowSums(patient.allocation)/sum(patient.allocation)
obs.tox.table <- rowSums(observe.tox)/sum(patient.allocation)

mean.cohortB = mean(num.group2.patients)
sd.cohortB = sd(num.group2.patients)
mean.duration = mean(study.duration)
sd.duration = sd(study.duration)

print(proc.time() - start.time)
result <- list(prior=prior, target.tox=target.tox, number.trials=number.trials, true.tox=true.tox, arrival.rate=arrival.rate, 
prop.B=prop.B, target.crm=target.crm, min.cohortB=min.cohortB, cycle.length=cycle.length, cohort.size=cohort.size, max.N=max.N, start.level=start.level, 
total.patients=total.patients, num.cohortA.patients=num.cohortA.patients, num.cohortB.patients=num.cohortB.patients,
num.group1.patients=num.group1.patients, num.group2.patients=num.group2.patients, results.num.dose.changes=result.num.dose.changes, MTD.selection=MTD.selection,
study.duration=study.duration, observe.tox=observe.tox, patient.allocation=patient.allocation,

MTD.selection.table = MTD.selection.table, true.MTD=true.MTD, PCS=PCS, obs.tox.overall=obs.tox.overall, patient.allocation.table=patient.allocation.table, obs.tox.table=obs.tox.table,
mean.cohortB=mean.cohortB, sd.cohortB=sd.cohortB, mean.duration=mean.duration, sd.duration=sd.duration)
return(result)
}
############################

#######
# Scenario 1
set.seed(7652)
crm1 <- target.crm(prior=c(0.05,0.1,0.2,0.3), target.tox=0.2, 
number.trials=100, true.tox=c(0.05,0.12,0.20,0.30), 
arrival.rate=15, prop.B=0.1, target.crm=1, min.cohortB=0, 
cycle.length=28, cohort.size=3, max.N=18, start.level=2)

crm1$MTD.selection.table
crm1$max.N
crm1$PCS
crm1$patient.allocation.table
crm1$obs.tox.table
crm1$mean.cohortB
crm1$sd.cohortB
crm1$mean.duration
crm1$sd.duration

# Scenario 2
set.seed(7652)
crm2 <- target.crm(prior=c(0.05,0.1,0.2,0.3), target.tox=0.2, 
number.trials=100, true.tox=c(0.20,0.30,0.40,0.50), 
arrival.rate=15, prop.B=0.1, target.crm=1, min.cohortB=0, 
cycle.length=28, cohort.size=3, max.N=18, start.level=1)

crm2$MTD.selection.table
crm2$max.N
crm2$PCS
crm2$patient.allocation.table
crm2$obs.tox.table
crm2$mean.cohortB
crm2$sd.cohortB
crm2$mean.duration
crm2$sd.duration


## Scenario 3
set.seed(45453)
crm3 <- target.crm(prior=c(0.05,0.1,0.2,0.3,0.4), target.tox=0.2, 
number.trials=100, true.tox=c(0.05, 0.08, 0.15, 0.22, 0.30), 
arrival.rate=30, prop.B=0.2, target.crm=1, min.cohortB=0, 
cycle.length=28, cohort.size=3, max.N=18, start.level=2)

crm3$MTD.selection.table
crm3$max.N
crm3$PCS
crm3$patient.allocation.table
crm3$obs.tox.overall
crm3$obs.tox.table
crm3$mean.cohortB
crm3$sd.cohortB
crm3$mean.duration
crm3$sd.duration

