#### Code Repository for statistical analyses conducted in Van Mierlo et al. 2022 "Occupancy of invasive Northern 
#### Crayfish (Faxonius virilis) in northern systems is driven primarily by tributary water temperature". Code was 
#### first written under R version 3.6.3 and was modified to run in the latest version of R (as of May 26th, 2022)
#### which is 4.2.0. 

#### 0. Load Required Packages ####
library(ecodist)
library(vegan)
library(ape)
library(magrittr)
library(vctrs)
library(pillar)
library(ellipsis)
require(tidyverse)
library(dplyr)
library(zoo)
library(MuMIn)
library(AICcmodavg)
require(unmarked)
library(ggplot2)
library(glmm)
library(readr)
library(glmmTMB)
library(readr)
library(DHARMa)

#### 1. Creation of Physical Complexity Gradient using PCoA ####
crayenv = read_csv("cray_env_modified1.csv")
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
sites.to.leave.out = c('', 'BUC1', 'CON1', 'POI1', 'VER2', 'VER3', 'SLA1')
crayenv = crayenv[!crayenv$site_code %in% sites.to.leave.out, ] #Excludes sites from data set
variables.to.keep = c('site_code','perc_macro', 'perc_rocky', 'woody', 'in_bould', 'o_h_b')
# note the above only includes physical complexity variables as that is all that is needed. 
cray.site.vars1 = crayenv[, variables.to.keep]
#Change all these to numeric for mean calculation
cray.site.vars1 = cbind(site_code = crayenv$site_code, cray.site.vars1 %>% mutate_all(as.numeric))
#Group by site and take the mean value over each trapline
cray.site.vars1 = aggregate(. ~ site_code, data = cray.site.vars1, FUN = mean, na.rm = TRUE, na.action = function(x) x)
cray.site.vars1[is.nan(cray.site.vars1)] = NA
data = cray.site.vars1
data
rownames(data)=data$site_code; # creates rownames and define the objects for ordination (sites)
sitelabels=data[,1] # creates labels that we need later
dat2=data[,2:6] # keeps numerical variables

# PCoA with Gower's Calculation and Plot
gower.dist = vegdist(dat2, method = "gower")
res = pcoa(gower.dist)
res$values # axis 1 explains 54.1% variance in the data.
res$vectors # values to be made into complexity gradient.Use Axis 1.  
par(mfrow=c(1,1))
biplot(res)
dat2 = apply(dat2, 2, scale, center=TRUE, scale=TRUE)
biplot(res, dat2)

# scree plot for PCoA
height = as.data.frame(res$values)
height = height$Eigenvalues
barplot(height, names = paste ('Axis', 1:37), las = 3, ylab = "eigenvalues")
# Produced 16 axes with positive eigen values. 
pcoa_vect = as.data.frame(res$vectors)

#### 2. Creation of unmarkedFrameOccu object for use in unmarked occuRN modeling ####
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
#Sites excluded from analysis due to too few environmental measurements. 
sites.to.leave.out = c('', 'BUC1', 'CON1', 'POI1', 'VER2', 'VER3', 'SLA1')
crayenv = crayenv[!crayenv$site_code %in% sites.to.leave.out, ] #Excludes sites from data set
cray = read_csv("cray_modified.csv")
cray = cray[!cray$site_code %in% sites.to.leave.out, ]#Excludes sites from data set
#Group observations by unique combinations of site + trapline, then add a new column indicating whether that trapline had a recorded crayfish and # crayfish per trapline
cray.sum.traps = cray %>% group_by(site_code, trap_line) %>% summarise(occupied = as.numeric(sum(number_cray) > 0),
                                                                       abundance= sum(number_cray)) 
cray.abundances = cray.sum.traps[, c('site_code', 'trap_line', 'abundance')] # adds column titles
#Manipulate the dataframe so each trapline is now an extra column
cray.abundances = cray.abundances %>% spread(trap_line, abundance) %>% as.data.frame
cray.sum.traps = cray.sum.traps[, c('site_code', 'trap_line', 'occupied')]
#Manipulate the dataframe so each trapline is now an extra column
cray.sum.traps = cray.sum.traps %>% spread(trap_line, occupied) %>% as.data.frame # creates a site and trpline matrix with just binary pres/absence data.
Y_occu = cray.sum.traps[, -1]# removes site code column so it can be stitched to dataframe later. 


variables.to.keep = c('site_code','temp', 'flow_ms', 'depth_cm', 'j_turb')
cray.site.vars = crayenv[, variables.to.keep]
#Change all these to numerics to calculate the mean
cray.site.vars<-transform(cray.site.vars,temp=as.numeric(as.character(temp)))
cray.site.vars<-transform(cray.site.vars,flow_ms=as.numeric(as.character(flow_ms)))
cray.site.vars<-transform(cray.site.vars,depth_cm=as.numeric(as.character(depth_cm)))
cray.site.vars<-transform(cray.site.vars,j_turb=as.numeric(as.character(j_turb)))
#Group by site and take the mean value over each trapline
cray.site.vars = aggregate(. ~ site_code, data = cray.site.vars, FUN = mean, na.rm = TRUE, na.action = function(x) x)
cray.site.vars = cbind(cray.site.vars, pcoa_vect) # Adds in the complexity gradient variable. 
cray.site.vars = cray.site.vars[,c(1:6)]
cray.site.vars = cray.site.vars %>% 
  rename(
    complex = Axis.1)
cray.site.names = cray.site.vars[,1]
cray.site.vars.wnames = cray.site.vars
cray.site.vars = cray.site.vars[, -1]# removes the site code column so it can be stitched to dataframe later.
#Convert NaN's to NA's
cray.site.vars[is.nan(cray.site.vars)] = NA
siteCovs_occu = cray.site.vars


cray.trapline.vars = crayenv[, variables.to.keep]
#Change all to numerics to calculate the mean. Additionally, "fill down" NA values for a site with the new value
cray.trapline.vars = cbind(site_code = crayenv$site_code, trapline = crayenv$trapline, 
                           cray.trapline.vars %>% mutate_all(as.numeric))
#Group by site and trapline, and take the mean value over each trapline
cray.trapline.vars = aggregate(. ~ site_code + trapline, data = cray.trapline.vars, FUN = mean, na.rm = TRUE, na.action = function(x) x)
#Convert NaN's to NA's
cray.trapline.vars[is.nan(cray.trapline.vars)] = NA
obsCovs_occu = lapply(X = variables.to.keep, FUN = function(var) {
  this.df = cbind(cray.trapline.vars[, c('site_code', 'trapline')], xvar = cray.trapline.vars[, var])
  this.df = this.df %>% group_by(site_code, trapline) %>% summarise(xvar = mean(xvar, na.rm = TRUE)) %>% spread(trapline, xvar) %>% as.data.frame
  #Convert NaN's to NA's
  this.df[is.nan(this.df)] = NA
  as.matrix(this.df[,-1])
})
names(obsCovs_occu) = variables.to.keep
#Add abundance data 
obsCovs_occu = c(obsCovs_occu, abundances = list(as.matrix(cray.abundances[,-1])))

occu.data = unmarkedFrameOccu(Y_occu, siteCovs_occu, obsCovs_occu)

#### 3. Global model creation & Goodness of Fit test ####
model = occuRN( ~ flow_ms+depth_cm ~ temp+flow_ms+j_turb+complex, data = occu.data, control = list(maxit = 5000))
summary(model)

mb.gof.test(model, nsim = 100, plot.hist = FALSE,
            report = NULL, parallel = TRUE, ncores,
            cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, maxK = NULL)

# Fit (p-value =0.04, c-hat = 2.22) is reasonable but needs to be accounted for using QAICc (to adjust for overdispersion/lack of fit)
# for model selection and inflating SE by root(c-hat) = root(2.22). (Mackenzie & Bailey , 2004)



#### 4. Iteration of Global model and Identification of best models using QAICc ####

#all iterations with QAICc and root c-hat inflated SE
all.models = dredge(model, trace = 2, rank ="QAICc", chat = 2.22) 
best.models.table = all.models[all.models$delta < 3, ]

best.models = model.avg(all.models, subset = (delta < 3))
sds = best.models$coefArray[, 2, ]
colnames(sds) = paste0(colnames(sds), "_se")
best.models.table= cbind(best.models.table, sds)

# Compute all top models from above table output
best1st = occuRN( ~ 1 ~ temp, data = occu.data, control = list(maxit = 5000))
best2nd = occuRN( ~ flow_ms ~ temp, data = occu.data, control = list(maxit = 5000))
best3rd = occuRN( ~ 1 ~ temp+j_turb, data = occu.data, control = list(maxit = 5000))
best4th = occuRN( ~ 1 ~ temp+complex, data = occu.data, control = list(maxit = 5000))
best5th = occuRN( ~ 1 ~ temp+flow_ms, data = occu.data, control = list(maxit = 5000))
best6th = occuRN( ~ depth_cm ~ temp, data = occu.data, control = list(maxit = 5000))



#### 5. Model averaging of P(lam) as a function of each occupancy covariate ####
##set up candidate models  (within deltaQAICc = 3)        
Cand.mod <- list( )
##Best 1st-6th models  
Cand.mod[[1]] <-occuRN( ~ 1 ~ temp, data = occu.data, control = list(maxit = 5000))
Cand.mod[[2]] <-occuRN( ~ flow_ms ~ temp, data = occu.data, control = list(maxit = 5000))
Cand.mod[[3]] <-occuRN( ~ 1 ~ temp+j_turb, data = occu.data, control = list(maxit = 5000))
Cand.mod[[4]] <-occuRN( ~ 1 ~ temp+complex, data = occu.data, control = list(maxit = 5000))
Cand.mod[[5]] <- occuRN( ~ 1 ~ temp+flow_ms, data = occu.data, control = list(maxit = 5000))
Cand.mod[[6]] <- occuRN( ~ depth_cm ~ temp, data = occu.data, control = list(maxit = 5000))

##assign names to each model
Modnames <- c("best1st", "best2nd", "best3rd", "best4th", "best5th", "best6th") 

##compute model-averaged estimates for parameters appearing on top models
temppred= modavg(parm = "temp
                 ", cand.set = Cand.mod, modnames = Modnames, c.hat = 2.22, parm.type ="lambda")
temppred$Parameter 
temppred$Mod.avg.beta 
temppred$Mod.avg.table

## Computing and graphing model averaged lam predictions as a function of temp
# tutorial followed can be found here: https://jamesepaterson.github.io/jamespatersonblog/2020-11-09_occupancy_part2.html

# First, set-up a new dataframe to predict along a sequence of the covariate.
# Predicting requires all covariates, so let's hold the other covariates constant at their mean value
mod.avg.data <- data.frame(temp = seq(min(cray.site.vars$temp), 
                                      max(cray.site.vars$temp), by = 0.5),
                           flow_ms = mean(cray.site.vars$flow_ms), # hold other variables constant
                           depth_cm = mean(cray.site.vars$depth_cm), # hold other variables constant
                           complex = mean(cray.site.vars$complex),
                           j_turb = mean(cray.site.vars$j_turb)) # hold other variables constant

# Model-averaged prediction of occupancy and confidence interval
mod_avg_pred <- modavgPred(Cand.mod,
                           c.hat = 2.22,   # to change variance inflation factor, default = 1) 
                           parm.type = "lambda", second.ord=TRUE, # lambda = occupancyRN , second.ord = TRUE for QAICc
                           newdata = mod.avg.data)[c("mod.avg.pred",
                                                     "lower.CL",
                                                     "upper.CL")]

# Put prediction, confidence interval, and covariate values together in a data frame
mod_avg_pred_table <- data.frame(Predicted = mod_avg_pred$mod.avg.pred,
                                 lower = mod_avg_pred$lower.CL,
                                 upper = mod_avg_pred$upper.CL,
                                 mod.avg.data)

# Plot the relationship
mod_avg_pred_plot_temp <- ggplot(mod_avg_pred_table, aes(x = temp, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, linetype = "solid", fill = "grey70") +
  geom_path(size = 1, linetype = "solid" ) +
  labs(x = "Temperature C", y = "Occupancy") +
  geom_hline(aes(yintercept=0.5),colour="grey40", linetype = "dotted") +
  geom_vline(aes(xintercept=18.6),colour="grey50", linetype = "longdash") +
  theme_classic() +
  coord_cartesian(ylim = c(0,1), xlim=c(9.99,21.5)) +
  theme(text = element_text( size=12, colour = "black"),
        axis.text = element_text(colour = "black"))
mod_avg_pred_plot_temp 

## Lam as a function of physical complexity
comppred= modavg(parm = "complex", cand.set = Cand.mod, modnames = Modnames, c.hat = 2.22, parm.type ="lambda")
comppred$Parameter 
comppred$Mod.avg.beta
comppred$Mod.avg.table

mod.avg.data <- data.frame(complex = seq(min(cray.site.vars$complex), 
                                         max(cray.site.vars$complex), by = 0.05), # change the "by" number according the variables range.
                           flow_ms = mean(cray.site.vars$flow_ms), # hold other variables constant
                           depth_cm = mean(cray.site.vars$depth_cm), # hold other variables constant
                           temp = mean(cray.site.vars$temp),
                           j_turb = mean(cray.site.vars$j_turb)) # hold other variables constant

# Model-averaged prediction of occupancy and confidence interval
mod_avg_pred <- modavgPred(Cand.mod,
                           c.hat = 2.22,   # to change variance inflation factor, default = 1) 
                           parm.type = "lambda", second.ord=TRUE, # lambda = occupancyRN , second.ord = TRUE for QAICc
                           newdata = mod.avg.data)[c("mod.avg.pred",
                                                     "lower.CL",
                                                     "upper.CL")]

# Put prediction, confidence interval, and covariate values together in a data frame
mod_avg_pred_table <- data.frame(Predicted = mod_avg_pred$mod.avg.pred,
                                 lower = mod_avg_pred$lower.CL,
                                 upper = mod_avg_pred$upper.CL,
                                 mod.avg.data)

# Plot the relationship
mod_avg_pred_plot <- ggplot(mod_avg_pred_table, aes(x = complex, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, linetype = "solid", fill = "grey70") +
  geom_path(size = 1, linetype = "solid" ) +
  labs(x = "Physical Complexity", y = "Occupancy") +
  geom_hline(aes(yintercept=0.5),colour="grey40", linetype = "dotted") +
  theme_classic() +
  coord_cartesian(ylim = c(0,1), xlim=c(-0.45,0.35)) +
  theme(text = element_text( size=12, colour = "black"),
        axis.text = element_text(colour = "black"))
mod_avg_pred_plot 
# note that the aes(x= , xlime=c(, and lans(x= needed to be changed to suit NMDS1

# Lam as a function of flow velocity 
comppred= modavg(parm = "flow_ms", cand.set = Cand.mod, modnames = Modnames, c.hat = 2.22, parm.type ="lambda")
comppred$Parameter 
comppred$Mod.avg.beta
comppred$Mod.avg.table

mod.avg.data <- data.frame(flow_ms = seq(min(cray.site.vars$flow_ms), 
                                         max(cray.site.vars$flow_ms), by = 0.01), # note: hve the change the "by" number according the varibles range.
                           complex = mean(cray.site.vars$complex), # hold other variables constant
                           depth_cm = mean(cray.site.vars$depth_cm), # hold other variables constant
                           temp = mean(cray.site.vars$temp),
                           j_turb = mean(cray.site.vars$j_turb)) # hold other variables constant

# Model-averaged prediction of occupancy and confidence interval
mod_avg_pred <- modavgPred(Cand.mod,
                           c.hat = 2.22,   # to change variance inflation factor, default = 1) 
                           parm.type = "lambda", second.ord=TRUE, # lambda = occupancyRN , second.ord = TRUE for QAICc
                           newdata = mod.avg.data)[c("mod.avg.pred",
                                                     "lower.CL",
                                                     "upper.CL")]

# Put prediction, confidence interval, and covariate values together in a data frame
mod_avg_pred_table <- data.frame(Predicted = mod_avg_pred$mod.avg.pred,
                                 lower = mod_avg_pred$lower.CL,
                                 upper = mod_avg_pred$upper.CL,
                                 mod.avg.data)

# Plot the relationship
mod_avg_pred_plot <- ggplot(mod_avg_pred_table, aes(x = flow_ms, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, linetype = "solid", fill = "grey70") +
  geom_path(size = 1, linetype = "solid" ) +
  labs(x = "Flow Velocity (m/s)", y = "Occupancy") +
  geom_hline(aes(yintercept=0.5),colour="grey40", linetype = "dotted") +
  theme_classic() +
  coord_cartesian(ylim = c(0,1), xlim=c(0.02,0.36)) +
  theme(text = element_text( size=12, colour = "black"),
        axis.text = element_text(colour = "black"))
mod_avg_pred_plot

### Lam as a function of turbidity
comppred= modavg(parm = "j_turb", cand.set = Cand.mod, modnames = Modnames, c.hat = 2.22, parm.type ="lambda")
comppred$Parameter 
comppred$Mod.avg.beta
comppred$Mod.avg.table
mod.avg.data <- data.frame(j_turb= seq(min(cray.site.vars$j_turb), 
                                       max(cray.site.vars$j_turb), by = 5), # note: have to change the "by" number according the variables' range.
                           complex = mean(cray.site.vars$complex), # hold other variables constant
                           depth_cm = mean(cray.site.vars$depth_cm), # hold other variables constant
                           temp = mean(cray.site.vars$temp),
                           flow_ms = mean(cray.site.vars$flow_ms)) # hold other variables constant

# Model-averaged prediction of occupancy and confidence interval
mod_avg_pred <- modavgPred(Cand.mod,
                           c.hat = 2.22,   # to change variance inflation factor, default = 1) 
                           parm.type = "lambda", second.ord=TRUE, # lambda = occupancyRN , second.ord = TRUE for QAICc
                           newdata = mod.avg.data)[c("mod.avg.pred",
                                                     "lower.CL",
                                                     "upper.CL")]

# Put prediction, confidence interval, and covariate values together in a data frame
mod_avg_pred_table <- data.frame(Predicted = mod_avg_pred$mod.avg.pred,
                                 lower = mod_avg_pred$lower.CL,
                                 upper = mod_avg_pred$upper.CL,
                                 mod.avg.data)

# Plot the relationship
mod_avg_pred_plot <- ggplot(mod_avg_pred_table, aes(x = j_turb, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, linetype = "solid", fill = "grey70") +
  geom_path(size = 1, linetype = "solid" ) +
  labs(x = "Turbidity (NTU)", y = "Occupancy") +
  geom_hline(aes(yintercept=0.5),colour="grey40", linetype = "dotted") +
  theme_classic() +
  coord_cartesian(ylim = c(0,1), xlim=c(0,170)) +
  theme(text = element_text( size=12, colour = "black"),
        axis.text = element_text(colour = "black"))
mod_avg_pred_plot 


#### 6. Model Validation using Leave-One-Out Cross Validation Method ####
null = occuRN( ~ 1 ~ 1, data = occu.data, control = list(maxit = 5000)) # null model to compare best models to.
fl <- fitList(best1st,best2nd, best3rd, best4th, best5th, best6th,null)
(cvlist <- crossVal(fl, method='leaveOneOut'))
# Result: Method: leave-one-out

#Root mean square error:
#  Estimate     SD
#best1st   0.3250 0.2351
#best2nd   0.3279 0.2402
#best3rd   0.3269 0.2406
#best4th   0.3378 0.2513
#best5th   0.3311 0.2375
#best6th   0.3282 0.2383
#null      0.4160 0.1885  # Best models are approx. 67% accurate in predicting occupancy and are 
# approx. 9% better at predicting occupancy than the null model. 

#### 7. Temperature Standardized Occupancy model to Investigate if Water tempertature is masking other relationships ####
standtemp18 <- subset(cray.site.vars, temp >= 18) # 19 sites and temp range of 4 degrees (18.1 -22.2)
#create new occu data frame with temp >+18 degree temp subset.
pcoa_axis1 = as.data.frame(pcoa_axis1)
pcoa_axis1

is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
#Sites excluded from analysis due to too few environmental measurments. 
sites.to.leave.out = c('', 'BUC1', 'CON1', 'POI1', 'VER2', 'VER3', 'SLA1')
crayenv = crayenv[!crayenv$site_code %in% sites.to.leave.out, ] #Excludes sites from data set
temp.sites.keep =c('BEA1', 'BMD1','BMD2','BMD3', 'EGG2','LIT1','MIS1','MOO1', 'RSS3','SITE45','SITE7','SLA2', 'STU2','VER4','VER5','WAS1','WMD2','WOU1','WSK5')
cray = read_csv("cray_modified.csv")
cray = cray[!cray$site_code %in% sites.to.leave.out, ]#Excludes sites from data set
#Group observations by unique combinations of site + trapline, then add a new column indicating whether that trapline had a recorded crayfish and # crayfish per trapline
cray.sum.traps = cray %>% group_by(site_code, trap_line) %>% summarise(occupied = as.numeric(sum(number_cray) > 0),
                                                                       abundance= sum(number_cray)) 
cray.abundances = cray.sum.traps[, c('site_code', 'trap_line', 'abundance')] # adds column titles
#Manipulate the dataframe so each trapline is now an extra column
cray.abundances = cray.abundances %>% spread(trap_line, abundance) %>% as.data.frame
cray.abundances =cray.abundances[cray.abundances$site_code %in% temp.sites.keep, ]
cray.sum.traps = cray.sum.traps[, c('site_code', 'trap_line', 'occupied')]
#Manipulate the dataframe so each trapline is now an extra column
cray.sum.traps = cray.sum.traps %>% spread(trap_line, occupied) %>% as.data.frame # creates a site and trpline matrix with just binary pres/absence data.
cray.sum.traps =cray.sum.traps[cray.sum.traps$site_code %in% temp.sites.keep, ]
Y_occu = cray.sum.traps[, -1]# removes site code column so it can be stitched to dataframe later. 

variables.to.keep = c('site_code','temp', 'flow_ms', 'depth_cm', 'j_turb')
cray.site.vars = crayenv[, variables.to.keep]
#Change all these to numerics so I can calculate the mean
cray.site.vars<-transform(cray.site.vars,temp=as.numeric(as.character(temp)))
cray.site.vars<-transform(cray.site.vars,flow_ms=as.numeric(as.character(flow_ms)))
cray.site.vars<-transform(cray.site.vars,depth_cm=as.numeric(as.character(depth_cm)))
cray.site.vars<-transform(cray.site.vars,j_turb=as.numeric(as.character(j_turb)))
#cray.site.vars = cbind(site_code = crayenv$site_code, cray.site.vars %>% mutate_all(as.numeric))
#Group by site and take the mean value over each trapline
cray.site.vars = aggregate(. ~ site_code, data = cray.site.vars, FUN = mean, na.rm = TRUE, na.action = function(x) x)
cray.site.vars = cbind(cray.site.vars, pcoa_vect) # Adds in the complexity gradient variable. 
cray.site.vars = cray.site.vars[,c(1:6)]
cray.site.vars = cray.site.vars %>% 
  rename(
    complex = Axis.1)
cray.site.names = cray.site.vars[,1]
cray.site.vars = cray.site.vars[, -1]# # removes the site code column so it can be stitched to dataframe later.
#Convert NaN's to NA's
cray.site.vars[is.nan(cray.site.vars)] = NA
cray.site.vars <- subset(cray.site.vars, temp >= 18)
siteCovs_occu = cray.site.vars

cray.trapline.vars = crayenv[, variables.to.keep]
#Change all these to numerics so I can calculate the mean. Additionally, "fill down" NA values for a site with the new value
cray.trapline.vars = cbind(site_code = crayenv$site_code, trapline = crayenv$trapline, 
                           cray.trapline.vars %>% mutate_all(as.numeric))
cray.trapline.vars =cray.trapline.vars[cray.trapline.vars$site_code %in% temp.sites.keep, ]
#Group by site and trapline, and take the mean value over each trapline
cray.trapline.vars = aggregate(. ~ site_code + trapline, data = cray.trapline.vars, FUN = mean, na.rm = TRUE, na.action = function(x) x)
#Convert NaN's to NA's
cray.trapline.vars[is.nan(cray.trapline.vars)] = NA
obsCovs_occu = lapply(X = variables.to.keep, FUN = function(var) {
  this.df = cbind(cray.trapline.vars[, c('site_code', 'trapline')], xvar = cray.trapline.vars[, var])
  this.df = this.df %>% group_by(site_code, trapline) %>% summarise(xvar = mean(xvar, na.rm = TRUE)) %>% spread(trapline, xvar) %>% as.data.frame
  #Convert NaN's to NA's
  this.df[is.nan(this.df)] = NA
  as.matrix(this.df[,-1])
})
names(obsCovs_occu) = variables.to.keep

#Add abundance data
obsCovs_occu = c(obsCovs_occu, abundances = list(as.matrix(cray.abundances[,-1])))

occu.data = unmarkedFrameOccu(Y_occu, siteCovs_occu, obsCovs_occu)



### 18 C+ Temperature standardized model
model = occuRN( ~ flow_ms+depth_cm ~ temp+flow_ms+j_turb+complex, data = occu.data, control = list(maxit = 5000))
summary(model)

mb.gof.test(model, nsim = 100, plot.hist = FALSE,
            report = NULL, parallel = TRUE, ncores,
            cex.axis = 1, cex.lab = 1, cex.main = 1,
            lwd = 1, maxK = NULL) # c-hat = 1.91

all.models = dredge(model, trace = 2, rank ="QAICc", chat = 1.91) 
best.models.table = all.models[all.models$delta < 2, ]

# models within 2 QAICc top model.
best.models = model.avg(all.models, subset = (delta < 2))
sds = best.models$coefArray[, 2, ]
colnames(sds) = paste0(colnames(sds), "_se")
best.models.table= cbind(best.models.table, sds)



#### 8. Dataset Creation for GLMM modeling####

# Create fixed variable data subset
crayenv = read_csv("cray_env_modified1.csv")
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
sites.to.leave.out = c('', 'BUC1', 'CON1', 'POI1', 'VER2', 'VER3', 'SLA1')
crayenv = crayenv[!crayenv$site_code %in% sites.to.leave.out, ] #Excludes sites from data set
variables.to.keep = c('site_code', 'temp', 'flow_ms', 'depth_cm', 'j_turb')
cray.site.vars = crayenv[, variables.to.keep]
#Change all these to numeric so I can calculate the mean
cray.site.vars<-transform(cray.site.vars,temp=as.numeric(as.character(temp)))
cray.site.vars<-transform(cray.site.vars,flow_ms=as.numeric(as.character(flow_ms)))
cray.site.vars<-transform(cray.site.vars,depth_cm=as.numeric(as.character(depth_cm)))
cray.site.vars<-transform(cray.site.vars,j_turb=as.numeric(as.character(j_turb)))

#Group by site and take the mean value over each trapline
cray.site.vars = aggregate(. ~ site_code, data = cray.site.vars, FUN = mean, na.rm = TRUE, na.action = function(x) x)
cray.site.vars = cbind(cray.site.vars, pcoa_vect) # Adds in the complexity gradient variable. 
cray.site.vars = cray.site.vars[,c(1:6)]
cray.site.vars = cray.site.vars %>% 
  rename(
    complex = Axis.1)
cray.site.vars[is.nan(cray.site.vars)] = NA
summary(cray.site.vars)

# Create dependent variable (abundance) data subset 
cray = read_csv("cray_modified.csv")
sites.to.leave.out = c('', 'BUC1', 'CON1', 'POI1', 'VER2', 'VER3', 'SLA1')
cray = cray[!cray$site_code %in% sites.to.leave.out, ]#Excludes sites from data set
cray<-transform(cray,cluster=as.factor(as.character(cluster)))
cray<-transform(cray,trap_line=as.factor(as.character(trap_line)))
cray<-transform(cray,trap=as.factor(as.character(trap)))
abun =aggregate(cray$number_cray, by=list(category = cray$site_code), FUN = sum, na.rm=TRUE)
summary(abun)
abun = abun %>% 
  rename(
    site_code = category,
    abundance= x
  )

# Combine 3 data subsets to form final data set
full <- abun
full<-merge(full, cray.site.vars, by="site_code")

#### 9. GLMM Modeling with Spatial (Site) Random Effect ####
global.pois <- glmmTMB(abundance ~ temp + j_turb + flow_ms + complex + (1|site_code), ziformula =~0, family=poisson(), data=full)
c_hat(global.pois) # c-hat' = 0.15 (method: pearson estimator).. so good fit and no overdispersion.

## Following Bolker interpretation advice and Noelle's test

# first check random effect structure
none <- glmmTMB(abundance ~ temp + j_turb + flow_ms + complex , ziformula =~0, family=poisson(), data=full) 
site <- glmmTMB(abundance ~ temp + j_turb + flow_ms + complex + (1|site_code) , ziformula =~0, family=poisson(), data=full)
AIC(none, site)
#      df      AIC
#none  5 353.0392
#site  6 144.4128
## Best model includes site 

# Final global model 
global <- glmmTMB(abundance ~ temp + j_turb + flow_ms + complex + (1|site_code) , ziformula =~0, family=poisson(), data=full)

# DHARMa work through: https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#workflow-in-dharma
library(DHARMa)

### Check residuals to determine if model is well specified 
# Simulate residuals: (n=250 simulations)
simres.global <- simulateResiduals(global)
# Check for overdispersion in the simulated residuals: 
testDispersion(simres.global) 
# Plot residuals: 
plot(simres.global) 
# In general: looks good.. like the model is well specified
# for the left plot which described over/under dispersion we have no significant curves in the data and n.s test result = good dispersion
# for the right plot which looks at dispersion around the quartiles is good and n.s for all. 

### Zero-inflation tests. 
par(mfrow = c(1,2))
plot(full$temp, full$abundance, xlab = "temp", ylab = "abun") 
plot(full$flow_ms, full$abundance, xlab = "flow", ylab = "abun") 
plot(full$j_turb, full$abundance, xlab = "turb", ylab = "abun") 
plot(full$complex, full$abundance, xlab = "NMDS1", ylab = "abun") 
# all of these plots show the line of y=0 observations along the enviro variable... which is normally an indicator
# of zero-infaltion... however zero inflation normally causes overdispersion which is not show to be an issue as 
# investigated above. Also the test below indicates only very slight zeroinflation. 
# Check for 0 inflation in the simulated residuals:
testZeroInflation(simres.global)# ratioObsSim = 1.0868 means there is very slight (even negligible) zero-inflation

# model output interpretation 
summary(global)
# Same relationships exist in the model.

##All Iterations of global model #
options(na.action = "na.fail")
all.models = dredge(global, trace = 2, rank ="AICc") # AICc because no overdispersion but small sample size.
best.models.table = all.models[all.models$delta < 2, ]

# models within 2 AICc fo top model.
best.models = model.avg(all.models, subset = (delta < 2))
sds = best.models$coefArray[, 2, ]
colnames(sds) = paste0(colnames(sds), "_se")
best.models.table= cbind(best.models.table, sds)

first <- glmmTMB(abundance ~ temp + (1|site_code) , ziformula =~0, family=poisson(), data=full)
second <- glmmTMB(abundance ~ temp + complex + (1|site_code) , ziformula =~0, family=poisson(), data=full)
third <- glmmTMB(abundance ~ temp + flow_ms + (1|site_code) , ziformula =~0, family=poisson(), data=full)

summary(first)
summary(second)
summary(third)

## Random Effect Interpretations 
global <- glmmTMB(abundance ~ temp + j_turb + flow_ms + complex + (1|site_code) , ziformula =~0, family=poisson(), data=full)
global.no <- glmmTMB(abundance ~ temp + j_turb + flow_ms + complex , ziformula =~0, family=poisson(), data=full)

cc <- fixef(global)[1] ## intercept
ss <- sqrt(unlist(VarCorr(global))) ## random effects SD
plogis(qnorm(c(0.025,0.975),mean=cc,sd=ss))

pchisq(2*(logLik(global)-logLik(global.no)),
       df=1,lower.tail=FALSE)/2


first <- glmmTMB(abundance ~ temp + (1|site_code) , ziformula =~0, family=poisson(), data=full)
first.no <- glmmTMB(abundance ~ temp , ziformula =~0, family=poisson(), data=full)
second <- glmmTMB(abundance ~ temp + complex + (1|site_code) , ziformula =~0, family=poisson(), data=full)
second.no <- glmmTMB(abundance ~ temp + complex , ziformula =~0, family=poisson(), data=full)
third <- glmmTMB(abundance ~ temp + flow_ms + (1|site_code) , ziformula =~0, family=poisson(), data=full)
third.no <- glmmTMB(abundance ~ temp + flow_ms , ziformula =~0, family=poisson(), data=full)

pchisq(2*(logLik(first)-logLik(first.no)),
       df=1,lower.tail=FALSE)/2
# 'log Lik.' 2.189011e-64 (df=3) - significant 
pchisq(2*(logLik(second)-logLik(second.no)),
       df=1,lower.tail=FALSE)/2
# 'log Lik.' 2.6635e-61 (df = 4)  - significant 
pchisq(2*(logLik(third)-logLik(third.no)),
       df=1,lower.tail=FALSE)/2
# 'log Lik.' 8.650524e-59 (df=4) - significant
#### 10. Investigation of VIF for coliniarity in models ####
# run section 1-4 before running next lines of code
vif(model, 'state') # temp = 1.449163, flow = 1.178380, turb = 1.448678, complex = 1.250033
vif(best1st, 'state') # won't calculate because only 1 variable (temp) so there cannot be any covariance.
vif(best2nd, 'state') # won't calculate because only 1 variable (temp) so there cannot be any covariance.
vif(best3rd, 'state') # temp 1.144791, turb 1.144791 
vif(best4th, 'state') # temp = 1.004115, complex = 1.004115 
vif(best5th, 'state') # temp = 1.000561, flow = 1.000561
vif(best6th, 'state') # won't calculate because only 1 variable so there cannot be any covariance.

# run section 6 before the next lines of code
best1stTempCont = occuRN( ~ 1 ~ 1, data = occu.data, control = list(maxit = 5000))
best2ndTempCont = occuRN( ~ 1 ~ temp + complex, data = occu.data, control = list(maxit = 5000))
best3rdTempCont = occuRN( ~ depth_cm ~ temp, data = occu.data, control = list(maxit = 5000))
vif(model, 'state') # Global model: temp = 1.449163, flow = 1.178380, turb = 1.448678, complex = 1.250033 
vif(best1stTempCont, 'state') # won't calculate because no variables variables in model, so there cannot be any covariance.
vif(best2ndTempCont, 'state') # temp = 1.013979, complex = 1.013979 
vif(best3rdTempCont, 'state') # won't calculate because only 1 variable (temp) so there cannot be any covariance.

#### 11. Other Small Analyses - Temp 2020 vs. Historic & Comparison to Lake Winnipeg Data ####
## 2020. Vs. Historic Water Temp
# run sections 1 -4 before running the following code.
#create dataset
nsrtemp = read_csv("nsr.temp.csv")
summary(nsrtemp)
temp.vars.keep = c('RIVER_BASIN_CODE', 'RIVER_SUB_BASIN_CODE', 'STATION_NO', 'STATION_NAME', 'STATION_DESCRIPTION', 'M_LATITUDE', 'M_LONGITUDE', 'month', 'SAMPLE_DATETIME', 'temp1', 'temp2','temp3')
nsrtemp = nsrtemp[, temp.vars.keep]
nsrtemp$location <- paste(nsrtemp$M_LATITUDE, nsrtemp$M_LONGITUDE, sep="_")
nsrtemp$location = as.factor(nsrtemp$location)
duplicated(nsrtemp$location)
indv.sites =unique(nsrtemp$location)
indv.sites # 226 indvidual sites in the data set
summer<-subset(nsrtemp, month == c(6:8)) # summer measurments data subset, observations n=726
july<-subset(summer, month == 7) # july measurments data subset, observations n=274

### Now I can calculate average summer and july water temperature at each sites within each dataset and then matchup 
### my actual sites to those averages afterward. 
summer =summer %>% mutate(temp = coalesce(temp1,temp2,temp3)) %>%
  select(RIVER_BASIN_CODE, RIVER_SUB_BASIN_CODE, STATION_NO, STATION_NAME, STATION_DESCRIPTION, M_LATITUDE, M_LONGITUDE, month, SAMPLE_DATETIME, location, temp)
july =july %>% mutate(temp = coalesce(temp1,temp2,temp3)) %>%
  select(RIVER_BASIN_CODE, RIVER_SUB_BASIN_CODE, STATION_NO, STATION_NAME, STATION_DESCRIPTION, M_LATITUDE, M_LONGITUDE, month,SAMPLE_DATETIME, location,  temp)

# calculate the date range in which these water temp samples were collected
library(stringr)
summer.dates=str_split_fixed(summer$SAMPLE_DATETIME, "/", 3)
summer.dates =as.data.frame(summer.dates)
summary(summer.dates$V3) #summer samples taken 1954-2019 
july.dates=str_split_fixed(july$SAMPLE_DATETIME, "/", 3)
july.dates =as.data.frame(july.dates)
summary(july.dates$V3) #july samples taken 1954-2019

summary(summer$temp)
summary(july$temp)
meansummer = summer %>% 
  group_by(location) %>% 
  summarise(average = mean(temp)) # Mean summer water temperature at each unique site n = 111
summary(meansummer$average) # min = 5.65  , mean = 16.15  , max = 22.27 for whole NSR June-Aug 1954-2019 
meanjuly = july %>% 
  group_by(location) %>% 
  summarise(average = mean(temp)) # Mean July water temperature at each unique site n = 83
summary(meanjuly$average) # min = 5.65  , mean = 17.26  , max = 24.80  for whole NSR july 1954-2019 

summary(cray.site.vars$temp) # min = 9.41  , mean = 17.23  , max = 22.21 snapshot 2020 temps
par(mfrow = c(2,2))
hist(meansummer$average)
hist(cray.site.vars$temp)
hist(meanjuly$average)

### The meansummer mean and max temps are very close to the mean and max of my snapshot dataset. 
### Also, when the histograms are compared, they have very similar slightly left skewed but mostly normal 
### distributions. 1. this is a valid comparison and 2. that the snapshot 
### temperature dataset is a good representation of water temperature.

### t-test of means b/t snapshot and NSR temps.
# do the stats t-test
t.test(summer$temp,cray.site.vars$temp) # insignificant. 

# create dataset so barchar can be plotted. 
keepme = c("site_code","temp")
cray.site.temp = cray.site.vars.wnames[,keepme]
cray.site.temp$group <- ("2020 NSR")
keep = c("location","temp")
summer1 = summer[,keep]
summer1 = summer1 %>%  rename(site_code = location)
summer1$group <- ("Historic NSR")
comb=rbind(summer1,cray.site.temp)

#bar chart
library(ggpubr)
ggboxplot(comb, x = "group", y = "temp", 
          color = "group", palette = "grey",
          ylab = "Temperature", xlab = " ")

### Comparison to Lake Winnipeg
#Load datas set
doug_report = read_csv("model_val.csv") # data set from Winnipeg river report Apendix F from Doug. 
ones = subset(doug_report, vir_pres_abs == "1")
summary(doug_report$temp)


## Overlay Validation data set on top of model averaged temp plot
mod_avg_valid_plot2 <-mod_avg_pred_plot_temp +geom_point() +geom_point(data=doug_report, shape=8, aes(x = temp, y = vir_pres_abs),colour='black')
mod_avg_valid_plot2


