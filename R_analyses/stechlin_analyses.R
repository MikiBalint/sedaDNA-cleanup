library(vegan)
library(knitr)
library(boral)
library(mvabund)
library(corrplot)
library(leaps)
library(car)
library(mgcv)
library(effects)
# library(geosphere)
# library(ape) # installed with ctv, infos here: http://www.phytools.org/eqg/Exercise_3.2/
# library(vegan3d)
# library(bvenn)

##### Read in all data
# Read all abundance data
EmblAssign = read.csv(file="../../Data/stechlin_assigned_190915.tab",
                      header=T, sep='|', row.names = 1)

# Read POP and elements data
POP_elem = read.csv(file = "../../../Data/Stechlin_pre-analysis/stechlin_pop_elements.csv", 
                    header = T, row.names = 1)

# Read experimental setup data
ExpSet = read.csv(file = "../../../Data/Stechlin_pre-analysis/sample_infos.csv",
                  header = T, row.names = 1)

##### Filter the DNA abundance data
# which columns have the status info for head, internal, singletons?
StatusEmbl = EmblAssign[,grepl("obiclean.status", names(EmblAssign))]

# summarize the status info for each seq. variant
# redundant with the kept columns
# EmblAssign = data.frame(EmblAssign,
#                         h_count = rowSums(StatusEmbl == "h", na.rm=T),
#                         s_count = rowSums(StatusEmbl == "s", na.rm=T), 
#                         i_count = rowSums(StatusEmbl == "i", na.rm=T))

# Keep only sequence variants that were seen at least once as 'head' or 
# the 'singleton' count is higher than the 'intermediate' count
EmblHead = EmblAssign[(EmblAssign$obiclean_headcount) > 0 | 
                        EmblAssign$obiclean_singletoncount > EmblAssign$obiclean_internalcount,]

##### Clean up negative controls
# Remove the maximum read number of a sequence variant 
# found in a negative control from every sample that 
# contains that variant

# Extraction controls
ExtCont = grep("sample.EXT", names(EmblHead))

# PCR controls
PCRCont = grep("sample.PCR", names(EmblHead))

# Multiplexing controls
MPXCont = grep("sample.MPX", names(EmblHead))

# Write out negative control assignments for more analysis
NegCont = cbind(name = EmblHead$scientific_name, 
                EmblHead[,c(ExtCont, PCRCont, MPXCont)])

# Aggregate according to taxon
NegRegate = aggregate(. ~ name, NegCont, sum, na.action = na.exclude)
rownames(NegRegate) = NegRegate$name
NegRegate = NegRegate[,2:length(names(NegRegate))]

# Remove taxa that were not seen in negatives
NegRegate = NegRegate[apply(NegRegate,1,sum) > 0,]

# Write for checking negatives
write.csv(file = "negative_control_identities.csv", NegRegate)

# Maximum number of reads in any control sample
MaxControl = apply(EmblHead[,c(ExtCont,PCRCont,MPXCont)], 1, max)

# Extract the highest read number of a sequence variant in a control 
# from every sample
EmblControlled = EmblHead

EmblControlled[,grep("sample", names(EmblHead))] <- 
  sweep(EmblHead[,grep("sample", names(EmblHead))], 1, MaxControl, "-")

# Set negative values to 0. Warnings are because the non-numeric cells
EmblControlled[EmblControlled < 0] <- 0

# Remove sequence variants with no reads left
EmblControlled = EmblControlled[apply(EmblControlled[,grep("sample", names(EmblHead))],1,sum) > 0,]

# I will work with OTUs, so this may be not useful
# # Categories for sequence variants
# # Species abundance table from samples
# SamGregate = data.frame(name = EmblControlled$scientific_name,
#                         EmblControlled[,grep("sample.ST", names(EmblControlled))])
# SamCounts = aggregate(. ~ name, SamGregate, sum, na.action = na.exclude)
# rownames(SamCounts) = SamCounts$name
# SamCounts = SamCounts[,2:85]
# SamCounts = SamCounts[apply(SamCounts,1,sum) > 0,]

# Write sequence variants of taxa in table
# write.csv(file = "stechlin_taxon_abund_matrix.csv", SamCounts)
# write.csv(file = "stechlin_seq-variant_abund_matrix.csv", SamGregate)
write.csv(file = "variant_sequences.csv", data.frame(name = rownames(EmblControlled),
                                                     sci = EmblControlled$scientific_name,
                                                     seq = EmblControlled$sequence))

##### Methodology predictors of DNA data
OTUCounts = EmblControlled[,grep("sample.ST", names(EmblControlled))]

# Keep OTUs with at least some read
OTUCounts = OTUCounts[apply(OTUCounts,1,sum) > 0,]

# The data for this includes only the samples, but not the controls,
# as these have special diversiy and compositional contraints.

# Transpose the abundance matrix
OTUCountsT = t(OTUCounts)

# The same for experimental predictors
ExpPredictor = ExpSet[grep("ST01", rownames(ExpSet)),]

# Check sample correspondence
data.frame(rownames(OTUCountsT), rownames(ExpPredictor))

# PCNM vectors of the plate coordinates
PlatePCNM = pcnm(dist(ExpPredictor$column, ExpPredictor$row, method="euc"))

# PCNM visualizations
par(mfrow = c(2,4), mar = c(4,4,3,1))
for (i in colnames(scores(PlatePCNM))){
  ordisurf(ExpPredictor[,c(7,8)], scores(PlatePCNM, choi=i), 
           bubble = 4, main = i)
}

# Add the PCNM vectors to the experimental predictors
ExpPredictor = data.frame(ExpPredictor, scores(PlatePCNM))

# Forward selection infos from here: http://www.statmethods.net/stats/regression.html

## 1. DNA concentrations
ConcModel <- regsubsets(conc ~ kit + person + extract_order + weight,
                        data=ExpPredictor, nbest=10)

# Plot the model r2-s
par(mfrow = c(1,2))
plot(ConcModel,scale="r2")
subsets(ConcModel, statistic="rsq")

# R2 of best model
ConcModel.summary = summary(ConcModel)
max(ConcModel.summary$rsq)
ConcModel.summary$outmat[ConcModel.summary$rsq == max(ConcModel.summary$rsq),]

# Test the final model
conc.lm1 = lm(conc ~ kit + person + extract_order + weight,
              data=ExpPredictor[5:96,])

conc.lm2 = lm(conc ~ kit + extract_order + person + weight,
              data=ExpPredictor[5:96,])

# ANOVA tables
summary(conc.lm1)
kable(anova(conc.lm1))
kable(anova(conc.lm2))

# Plots
par(mar=c(4,4,1,1))

# Personnel differences
boxplot(conc ~ person, data=ExpPredictor[5:96,], 
        ylab = "DNA concentrations (ng/ul)", col="grey", boxwex=.5, notch=T)

# Kit differences
boxplot(conc ~ kit, data=ExpPredictor[5:96,], 
        ylab = "DNA concentrations (ng/ul)", col="grey", boxwex=.5, notch=T)

# Weight differences
plot(effect("weight", conc.lm1, multiline=T))

# Extraction order differences
plot(effect("extract_order", conc.lm1, multiline=T))

# plot(ExpPredictor[,"weight"], ExpPredictor[,"conc"], pch = 19,
#      xlab = "Sediment weight (mg)", ylab = "DNA concentration (ng/ul)")
# 
# abline(lm(conc ~ weight, data = ExpPredictor[5:96,]), lwd = 2)
# 
# plot(ExpPredictor[5:96,"extract_order"], ExpPredictor[5:96,"conc"], pch = 19,
#      xlab = "Extraction order", ylab = "DNA concentration (ng/ul)")
# abline(lm(conc ~ extract_order, data = ExpPredictor[5:96,]), lwd = 2)

## 2. Taxon richness and diversity

# Calculate Hill's series
CountsHill = renyi(OTUCountsT, hill = T)

# 2.1. Models of Hill's H0
H0.Model <- regsubsets(CountsHill$`0` ~ reads + kit + person + extract_order + 
                         pcr_order + weight + conc +
                         PCNM1 + PCNM2 + PCNM3 + PCNM4 + 
                         PCNM5 + PCNM6 + PCNM7,
                       data=ExpPredictor, nbest=10)

# Plot the model r2-s
summary.H0Model = summary(H0.Model)
plot(H0.Model,scale="r2")
subsets(H0.Model, statistic="rsq")

# Model with highest R2
summary.H0Model$outmat[summary.H0Model$rsq == max(summary.H0Model$rsq),]

# This model is way too complex, and many other models have similar rsq.
# I take the model with only 4 predictors that has rsq > 0.65
H0.compare = data.frame(subset = rownames(summary.H0Model$which), summary.H0Model$rsq, summary.H0Model$which)
# reads + kit + conc + PCNM3 + PCNM4

# # Similarly good models
# summary.H0Model$outmat[summary.H0Model$rsq > max(summary.H0Model$rsq) - 0.02,]

# Selected model:
H0.lm1 = lm(CountsHill$`0` ~ reads + kit + conc + PCNM3 + PCNM4, 
            data=ExpPredictor)

# H0 model statistics
summary(H0.lm1)
anova(H0.lm1)
kable(anova(H0.lm1))

# Plots of effects
plot(allEffects(H0.lm1))

# PCNM and richness
par(mar = c(2,2,1,1), mfrow = c(1,1))
plot(ExpPredictor$column, ExpPredictor$row, cex=CountsHill$`0`/150, pch=19)
ordisurf(ExpPredictor[,c(7,8)], scores(PlatePCNM, choi="PCNM4"), add = T)

coef(H0.lm1)
confint(H0.lm1)

# 2.2. Models of Hill's H1
H1.Model <- regsubsets(CountsHill$`1` ~ reads + kit + person + extract_order + 
                         pcr_order + weight + conc +
                         PCNM1 + PCNM2 + PCNM3 + PCNM4 + 
                         PCNM5 + PCNM6 + PCNM7,
                       data=ExpPredictor, nbest=10)

# Plot the model r2-s
par(mfrow = c(1,1), mar = c(4,4,1,1))
summary.H1Model = summary(H1.Model)
# plot(H1.Model,scale="r2")
subsets(H1.Model, statistic="rsq")

# Model with highest R2
summary.H1Model$outmat[summary.H1Model$rsq == max(summary.H1Model$rsq),]

# I take the model with only 4 predictors. rsq = 0.35877466
H1.compare = data.frame(subset = rownames(summary.H1Model$which), summary.H1Model$rsq, summary.H1Model$which)
# reads + person + pcr_order + weight + PCNM2

# Final model
H1.lm1 = lm(CountsHill$`1` ~ reads + person + pcr_order + weight + PCNM2, data=ExpPredictor)

# H2 model statistics
summary(H1.lm1)
anova(H1.lm1)
kable(anova(H1.lm1))

# Plots of effects
plot(allEffects(H1.lm1))

# plot(effect("pcr_order", H1.lm1, multiline=TRUE), ylab = "Hill's H1")

# 2.3. Models of Hill's H2
H2.Model <- regsubsets(CountsHill$`2` ~ reads + kit + person + extract_order + 
                         pcr_order + weight + conc +
                         PCNM1 + PCNM2 + PCNM3 + PCNM4 + 
                         PCNM5 + PCNM6 + PCNM7,
                       data=ExpPredictor, nbest=10)

# Plot the model r2-s
# par(mfrow = c(1,2), mar = c(4,4,1,1))
summary.H2Model = summary(H2.Model)
# plot(H2.Model,scale="r2")
subsets(H2.Model, statistic="rsq")

# Model with highest R2
# summary.H2Model$outmat[summary.H2Model$rsq == max(summary.H2Model$rsq),]

# Model with few variables but high Rsq Rsq = 
H2.compare = data.frame(subset = rownames(summary.H2Model$which), summary.H2Model$rsq, summary.H2Model$which)
# kit + person + extract_order + weight + conc

H2.lm1 = lm(CountsHill$`2` ~ kit + person + extract_order + weight + conc, data=ExpPredictor)

# H2 model statistics
summary(H2.lm1)
anova(H2.lm1)
kable(anova(H2.lm1))

# Plots of effects
# plot(effect("extract_order", H2.lm1, multiline=TRUE), ylab = "Hill's H2")
plot(allEffects(H2.lm1))

## 3. Community composition

# ### 3. Define core OTUs
# ## Summarize reads
# TotCount = apply(OTUCountsT,2,sum)
# 
# ## The average read number of OTUs
# MeanCount=apply(OTUCountsT,2,function(vec) mean(vec[vec>0]))
# 
# ## In how many samples is an OTU present?
# TotPresent = apply(OTUCountsT,2,function(vec) sum(vec>0))
# 
# ## The highest read number of an OTU in a sample
# MaxCount=apply(OTUCountsT,2,max)
# 
# ## Plotting incidence against abundance
# plot(TotPresent, MaxCount, xlab="Incidence",
#      ylab="Maximum Abundance", pch=20)
# 
# plot(TotPresent, log(MaxCount), xlab="Incidence",
#      ylab="log(Maximum Abundance)", pch=20)
# 
# ## Create a smoothed trendline
# gam1 = gam(log(MaxCount)~s(TotPresent))
# 
# plot(gam1, residuals=T, shade=T, rug=F, cex=2.6,
#      xlab="Incidence", ylab="logMean Abundance") # , xaxp=c(0,150,15)
# 
# ## consider OTUs as core if present in at least 20 samples
# IsFreq = TotPresent >= 20
# OTU.some = OTUCountsT[,IsFreq]

core.mvabund = mvabund(OTUCountsT)

# Manual stepwise model selection. Predictor order established with
# assumptions about potential laboratory error importance, or the
# experiences from the previous tests. Predictors considered: only those that were at least 
# once significant in the previous models
# reads + kit + person + extract_order + weight + conc + PCNM4

# 3.1. Stepwise commmunity model selection

core.m1 = manyglm(core.mvabund ~ reads + kit + person + extract_order + weight + conc + PCNM4, 
                  data = ExpPredictor, family = "negative.binomial")

core.m2 = manyglm(core.mvabund ~ reads + kit + person + extract_order + weight + conc, 
                  data = ExpPredictor, family = "negative.binomial")

anova(core.m1, core.m2, nBoot = 10)

core.m3 = manyglm(core.mvabund ~ reads + kit + person + extract_order + weight, 
                  data = ExpPredictor, family = "negative.binomial")

anova(core.m2, core.m3, nBoot = 10)

core.m4 = manyglm(core.mvabund ~ reads + kit + person + extract_order, 
                  data = ExpPredictor, family = "negative.binomial")

anova(core.m3, core.m4, nBoot = 10)

core.m5 = manyglm(core.mvabund ~ reads + kit + person, 
                  data = ExpPredictor, family = "negative.binomial")

anova(core.m4, core.m5, nBoot = 10)

core.m6 = manyglm(core.mvabund ~ reads + kit + extract_order, 
                  data = ExpPredictor, family = "negative.binomial")

anova(core.m4, core.m6, nBoot = 10)

core.m7 = manyglm(core.mvabund ~ reads + person + extract_order, 
                  data = ExpPredictor, family = "negative.binomial")

anova(core.m4, core.m7, nBoot = 10)

core.m8 = manyglm(core.mvabund ~ person + extract_order, 
                  data = ExpPredictor, family = "negative.binomial")

anova(core.m4, core.m8, nBoot = 10)

# A promising community model has reads + kit + person + extract_order
m4.anova = anova(core.m4, nBoot = 100, p.uni = "adjusted")
m4.anova$table

# 
# 
# methods.manyglm2 = manyglm(Samples.mvabund ~ reads + kit, 
#                            data = ExpPredictor, family = "negative.binomial")
# 
# methods.manyglm1$AICsum
# methods.manyglm2$AICsum
# 
# methods.manyglm3 = manyglm(Samples.mvabund ~ reads + person, 
#                            data = ExpPredictor, family = "negative.binomial")
# methods.manyglm3$AICsum
# 
# methods.manyglm4 = manyglm(Samples.mvabund ~ reads + person + extract_order, 
#                            data = ExpPredictor, family = "negative.binomial")
# methods.manyglm4$AICsum
# 
# methods.manyglm5 = manyglm(Samples.mvabund ~ reads + person + extract_order +
#                              pcr_order, 
#                            data = ExpPredictor, family = "negative.binomial")
# methods.manyglm5$AICsum
# 
# methods.manyglm6 = manyglm(Samples.mvabund ~ reads + person + extract_order +
#                              weight, 
#                            data = ExpPredictor, family = "negative.binomial")
# methods.manyglm6$AICsum
# 
# methods.manyglm7 = manyglm(Samples.mvabund ~ reads + person + extract_order +
#                              weight + conc, 
#                            data = ExpPredictor, family = "negative.binomial")
# methods.manyglm7$AICsum
# 
# methods.manyglm8 = manyglm(Samples.mvabund ~ reads + person + extract_order +
#                              weight + conc + PCNM1, 
#                            data = ExpPredictor, family = "negative.binomial")
# methods.manyglm8$AICsum
# 
# methods.manyglm9 = manyglm(Samples.mvabund ~ reads + person + extract_order +
#                              weight + conc + PCNM2, 
#                            data = ExpPredictor, family = "negative.binomial")
# methods.manyglm9$AICsum
# 
# methods.manyglm10 = manyglm(Samples.mvabund ~ reads + person + extract_order +
#                              weight + conc + PCNM3, 
#                            data = ExpPredictor, family = "negative.binomial")
# methods.manyglm10$AICsum
# 
# methods.manyglm11 = manyglm(Samples.mvabund ~ reads + person + extract_order +
#                              weight + conc + PCNM4, 
#                            data = ExpPredictor, family = "negative.binomial")
# methods.manyglm11$AICsum
# 
# methods.manyglm12 = manyglm(Samples.mvabund ~ reads + person + extract_order +
#                               weight + conc + PCNM5, 
#                             data = ExpPredictor, family = "negative.binomial")
# methods.manyglm12$AICsum
# 
# methods.manyglm13 = manyglm(Samples.mvabund ~ reads + person + extract_order +
#                               weight + conc + PCNM6, 
#                             data = ExpPredictor, family = "negative.binomial")
# methods.manyglm13$AICsum
# 
# methods.manyglm14 = manyglm(Samples.mvabund ~ reads + person + extract_order +
#                               weight + conc + PCNM7, 
#                             data = ExpPredictor, family = "negative.binomial")
# methods.manyglm14$AICsum
# 
# # The best model is 
# methods.manyglm7$call
# methods.manyglm7$AICsum
# plot(methods.manyglm7, which=c(1:4))
# 
# # 3.2. Model test
# community.anova = anova(methods.manyglm7, nBoot = 100, p.uni = "adjusted")
# kable(community.anova$table)
# community.summary = summary(methods.manyglm7, nBoot = 100, test = "LR")
# kable(community.summary$coefficients)

# model-based ordination
# Summary of overdispersion parameters
summary(methods.manyglm7$theta)
hist(methods.manyglm7$theta)

# Set overdispersion prior
set.prior = list(type = c("normal","normal","normal","uniform"),
                 hypparams = c(100, 20, 100, 20))

# LV ordination done on species with relatively low overdispersion
comm.ord = boral(OTUCountsT[,core.m4$theta < 50],
                 family = "negative.binomial", 
                 prior.control = set.prior, num.lv = 2, n.burnin = 10, 
                 n.iteration = 100, n.thin = 1)


# Plot the ordinations
par(mfrow = c(2,2), mar = c(2,2,2,1))
predictors = c("reads", "extract_order", "weight", "conc")
for (i in predictors) {
  ordicomm = ordiplot(comm.ord$lv.median, choices = c(1,2), type = "none", cex =0.5,
                      display = "sites", xlim = c(-0.3,0.3))
  points(ordicomm,"sites", pch=20, col=as.numeric(ExpPredictor$person))
  ordisurf(ordicomm, ExpPredictor[,i], add=T, col = "black", main=i)
}

# eddig

# ordiellipse(ordi.comm, ExpPredictor$person,cex=.5, 
#             draw="polygon", col="orange",
#             alpha=200,kind="se",conf=0.95, 
#             show.groups=(c("Miki")))
# ordiellipse(ordi.comm, ExpPredictor$person,cex=.5, 
#             draw="polygon", col="purple",
#             alpha=200,kind="se",conf=0.95, 
#             show.groups=(c("Orsi")))

##### Pesticide, elements and eDNA visualizations
# depths corresponding the samples
depths = read.csv(file="depths.csv", header=T)

# Transposed species abundance matrix
SamCountsT = t(SamCounts)

# Z-scores of counts
SamZ = scale(SamCountsT, center = T, scale = T)

# Center the POP and element data to get z-scores
PE_scaled = cbind(POP_elem[,1:2], 
                  scale(POP_elem[,3:length(names(POP_elem))], 
                        center = T, scale = T))

# Plot together with elements and POPs
par(mfrow = c(3,1), mar = c(2,4,1,1))
# POP
plot(PE_scaled$depth, seq(min(PE_scaled[,3:ncol(PE_scaled)], na.rm=T), 
                          max(PE_scaled[,3:ncol(PE_scaled)], na.rm=T), 
                          (max(PE_scaled[,3:ncol(PE_scaled)], na.rm=T) - 
                             min(PE_scaled[,3:ncol(PE_scaled)], na.rm=T))/(nrow(PE_scaled)-1)), 
     type = "n",
     ylab = "Pesticides", xlab = "Depth (cm)", main = "")
for (i in names(PE_scaled[3:14])){
  points(PE_scaled$depth, PE_scaled[,i], type = "l", col = sample(1:58), lwd=2)
}

# Elements
plot(PE_scaled$depth, seq(min(PE_scaled[,3:ncol(PE_scaled)], na.rm=T), 
                          max(PE_scaled[,3:ncol(PE_scaled)], na.rm=T), 
                          (max(PE_scaled[,3:ncol(PE_scaled)], na.rm=T) - 
                             min(PE_scaled[,3:ncol(PE_scaled)], na.rm=T))/(nrow(PE_scaled)-1)), 
     type = "n",
     ylab = "Elements", xlab = "Depth (cm)", main = "")
for (i in names(PE_scaled[15:31])){
  points(PE_scaled$depth, PE_scaled[,i], type = "l", col = sample(1:58), lwd=2)
}

# Taxa
plot(PE_scaled$depth, seq(min(SamZ, na.rm=T), # keeping PE_scaled$depth: to use the same depth range
                          max(SamZ, na.rm=T), 
                          (max(SamZ, na.rm=T) - 
                             min(SamZ, na.rm=T))/(nrow(PE_scaled)-1)), 
     type = "n",
     ylab = "Taxa", xlab = "Depth (cm)", main = "")
for (i in names(as.data.frame(SamZ))) {
  points(depths$depth, as.data.frame(SamZ)[,i], type = "p", col = sample(1:66), lwd=2)
}

corrplot(cor(cbind(depth = depths$depth, SamCountsT), method="spearman"),
         tl.cex = 0.5, tl.col = "black")



##### Analysis of community time series
# Only credible reads in samples are used for ecological analyses. 







# frogs, humans, other aquatic stuff, etc.
# Amphi, HighGroup, "Homo sapiens", FarmAnim, Bird, Fish, Insect, Mammal
Amphi = c("Dendropsophus leucophyllatus","Dendropsophus melanargyreus",
          "Dendropsophus minutus","Dendropsophus nanus","Dermatonotus muelleri",
          "Elachistocleis sp.","Eupemphix nattereri","Hylinae","Hyloidea",
          "Hypsiboas punctatus","Hypsiboas raniceps","Leptodactylus",
          "Leptodactylus chaquensis","Leptodactylus elenae","Leptodactylus fuscus",
          "Leptodactylus podicipinus","Leptodactylus syphax","Leptodactylus vastus",
          "Osteocephalus","Osteocephalus taurinus","Phyllomedusa azurea",
          "Phyllomedusa boliviana","Physalaemus albonotatus","Physalaemus centralis",
          "Physalaemus cuvieri","Pseudis limellum","Pseudis paradoxa",
          "Pseudopaludicola mystacalis","Rhinella schneideri","Scinax fuscomarginatus",
          "Scinax fuscovarius","Scinax nasicus","Scinax ruber","Scinax sp. FB-2014a",
          "Bufonidae","Elachistocleis","Gastrophryninae","Hylidae",
          "Leptodactylus latinasus","Melanophryniscus","Pseudis","Pseudis laevis","Scinax")

HighGroup = c("root", "Bilateria", "Amniota")

FarmAnim = c("Canis", "Sus scrofa", "Gallus", "Mus", "Bos", "Capreolus capreolus")
# write.csv(file = "other_taxa.csv", 
#           MetaHead$sci_name[!MetaHead$sci_name %in% 
#                               c(Amphi, HighGroup, "Homo sapiens", FarmAnim)])

Bird = c("Meleagris gallopavo", "Jacana", "Gallus")

Fish = c("Astyanax","Brachyplatystoma rousseauxii","Callichthys callichthys","Characidae",
         "Characiphysae","Opisthonema oglinum","Salmonidae","Sardinops","Sardinops sagax",
         "Siluriformes","Siluroidei","Synbranchus marmoratus","Thraupidae","Thunnini")

Insect = c("Micrathyria ocellata", "Libellulidae", "Micrathyria", "Tabanus", "Tramea")

Mammal = c("Homo sapiens", "Canis", "Sus scrofa", "Mus", "Bos", "Capreolus capreolus",
           "Myotis", "Laurasiatheria")

# Proportions of different groups
par(mar=c(1,1,1,1))
pie(c(sum(SummedControlled[MetaHead$sci_name %in% Amphi,]), 
      sum(SummedControlled[MetaHead$sci_name %in% Mammal,]),
      sum(SummedControlled[MetaHead$sci_name %in% Fish,]),
      sum(SummedControlled[MetaHead$sci_name %in% Insect,]),
      sum(SummedControlled[MetaHead$sci_name %in% Bird,]),
      sum(SummedControlled[MetaHead$sci_name %in% HighGroup,])),
    labels = paste(c("Frogs\n2 198 212", "Mammals\n38 568", 
                     "Fish\n1 873 599", "Insects\n314 910", 
                     "Birds\n4 785", "Higher groups\n1 384 940")),
    col = c(gray(0.9), gray(0.75), gray(0.60), gray(0.45), gray(0.30), gray(0.15)))

sum(SummedControlled[MetaHead$sci_name %in% Amphi,])+ 
  sum(SummedControlled[MetaHead$sci_name %in% Mammal,])+
  sum(SummedControlled[MetaHead$sci_name %in% Fish,])+
  sum(SummedControlled[MetaHead$sci_name %in% Insect,])+
  sum(SummedControlled[MetaHead$sci_name %in% Bird,])+
  sum(SummedControlled[MetaHead$sci_name %in% HighGroup,])

# Farm animals:
sum(SummedControlled[MetaHead$sci_name %in% FarmAnim,])

# Humans:
sum(SummedControlled[MetaHead$sci_name == "Homo sapiens",])

# Libellulidae: 
sum(SummedControlled[MetaHead$sci_name %in% c("Micrathyria ocellata", "Libellulidae", 
                                              "Micrathyria", "Tramea"),])







# Separate the assignment info and the abundances, so the sequence variants
# can be aggregated according to species

# Sequence assignment info
# AssignmentInfo = data.frame(best_id = EmblControlled$best_identity.db_i18S_V9_embl125,
#                             best_match = EmblControlled$best_match, 
#                       count = EmblControlled$count, 
#                       family = EmblControlled$family_name,
#                       genus = EmblControlled$genus_name, 
#                       rank = EmblControlled$rank,
#                       sci_name = EmblControlled$scientific_name, 
#                       taxid = EmblControlled$taxid, 
#                       sequence = EmblControlled$sequence,
#                       seq_names = rownames(EmblControlled),
#                       row.names=10)
# 
# # Abundance data: get all columns that have "sample" in column names
# AbundControlled = EmblControlled[, grepl("sample.ST", names(EmblControlled))]
# 
# # Controlled positive controls
# PosControlled = EmblControlled[, grepl("sample.POS", names(EmblControlled))]


# # combine the replicates of samples
# # get sample names, code from here: http://stackoverflow.com/questions/9704213/r-remove-part-of-string
# SampleNames = levels(as.factor(sapply(strsplit(names(AbundControlled), 
#                                                split='16S', fixed=TRUE), 
#                                       function(x) (x[1]))))
# 
# # Sum the replicates for each sample
# SummedReps = data.frame(row.names = rownames(AbundControlled))
# for (i in 1:length(SampleNames)){
#   ActualSet = grep(SampleNames[i], names(AbundControlled)) # grep the columns of interest
#   SummedReps = cbind(SummedReps, apply(AbundControlled[ActualSet], 1, sum))
# }
# colnames(SummedReps) = SampleNames
# # write.csv(file="test.csv", AbundControlled)
# 
# # In how many replicates observed per sample?
# PresentReps = data.frame(row.names = rownames(AbundControlled))
# for (i in 1:length(SampleNames)){
#   ActualSet = grep(SampleNames[i], names(AbundControlled))
#   Selected = AbundControlled[ActualSet]
#   Selected[Selected > 0] <- 1 # set the read numbers to 1
#   PresentReps = cbind(PresentReps, apply(Selected, 1, sum))
# }
# colnames(PresentReps) = SampleNames
# 
# # Set read numbers to 0 in a sample if the sequence variant was not observed in at least
# # two PCR replicates
# SummedControlled = SummedReps
# SummedControlled[PresentReps < 2] <- 0

# Aggregate frog sequence variants according to the species 
FrogAggregate = data.frame(name = MetaHead$sci_name[MetaHead$sci_name %in% Amphi], 
                           SummedControlled[MetaHead$sci_name %in% Amphi,])
FrogCounts = aggregate(. ~ name, FrogAggregate, sum, na.action = na.exclude)
rownames(FrogCounts) = FrogCounts$name
FrogCounts = FrogCounts[,2:96]

# Correct the Scinax madeirae: Scinax sp. FB-2014a is actually S. madeirae
# Species list
rownames(FrogCounts)[34] <- "Scinax madeirae"

# Remove taxa with no observations 
# (can happen due to the cleanup steps on the sequence variants)
FrogCounts = FrogCounts[apply(FrogCounts, 1, sum) != 0, ]

# Generate site metadata from column names
# Lake codes
Lakes = substring(names(FrogCounts), 8, 9)

# Site metadata: lakes, preservation/extraction methods, read numbers
SiteMeta = data.frame(Lakes = Lakes, 
                      Methods = c(rep(c("Control"),31), 
                                  rep(c("A","B","C"),12),
                                  rep("A",4), rep("B",4), rep("C",3),
                                  rep("A",3), rep("B",3), rep("C",3),
                                  rep("A",3), rep("B",3), rep("C",2)),
                      Reads = apply(FrogCounts, 2, sum))
LakeCodes = c("T1", "T2", "T3", "T4", "T5")

# Sum up species counts by sites
FrogLakes = data.frame(row.names = rownames(FrogCounts))
for (i in 1:length(LakeCodes)) {
  ActualSet = grep(LakeCodes[i], names(FrogCounts))
  FrogLakes = cbind(FrogLakes, apply(FrogCounts[ActualSet], 1, sum))
}
colnames(FrogLakes) = LakeCodes

# Nice species table
kable(cbind(FrogLakes, "Total (with controls)" = apply(FrogCounts, 1, sum)))

#######
# Venn diagrams (http://stackoverflow.com/questions/11722497/how-to-make-venn-diagrams-in-r)
RegioList = c("Ameerega picta",
              "Ceratophrys sp.",
              "Chiasmocleis albopunctata",
              "Dendropsophus cachimbo",
              "Dendropsophus leali",
              "Dendropsophus leucophyllatus",
              "Dendropsophus melanargyreus",
              "Dendropsophus minutus",
              "Dendropsophus nanus",
              "Dendropsophus salli",
              "Dermatonotus muelleri",
              "Elachistocleis sp.",
              "Eupemphix nattereri",
              "Hypsiboas geographicus",
              "Hypsiboas punctatus",
              "Hypsiboas raniceps",
              "Leptodactylus cf. didymus",
              "Leptodactylus cf. diptyx",
              "Leptodactylus elenae",
              "Leptodactylus fuscus",
              "Leptodactylus leptodactyloides",
              "Leptodactylus mystacinus",
              "Leptodactylus syphax",
              "Leptodactylus vastus",
              "Oreobates heterodactylus",
              "Osteocephalus taurinus",
              "Phyllomedusa azurea",
              "Phyllomedusa boliviana",
              "Physalaemus centralis",
              "Physalaemus albonotatus",
              "Physalaemus cuvieri",
              "Pseudis paradoxa",
              "Pseudopaludicola mystacalis",
              "Rhinella cf. paraguayensis",
              "Rhinella schneideri",
              "Rhinella mirandaribeiroi",
              "Scinax fuscomarginatus",
              "Scinax fuscovarius",
              "Scinax madeirae",
              "Scinax nasicus",
              "Scinax ruber",
              "Sphaenorhynchus lacteus",
              "Trachycephalus cf. typhonius")
#####           

SurveyList = c("Dendropsophus leucophyllatus",
               "Dendropsophus melanargyreus",
               "Dendropsophus minutus",
               "Dendropsophus nanus",
               "Dendropsophus salli",
               "Elachistocleis sp.",
               "Eupemphix nattereri",
               "Hypsiboas geographicus",
               "Hypsiboas punctatus",
               "Hypsiboas raniceps(d)",
               "Leptodactylus fuscus",
               "Leptodactylus syphax",
               "Leptodactylus vastus",
               "Phyllomedusa azurea",
               "Phyllomedusa boliviana",
               "Physalaemus albonotatus",
               "Physalaemus centralis",
               "Pseudis paradoxa(d)",
               "Pseudopaludicola mystacalis",
               "Sphaenorhynchus lacteus",
               "Scinax fuscomarginatus",
               "Scinax fuscovarius",
               "Scinax madeirae",
               "Scinax ruber",
               "Scinax nasicus")

eDNAList = rownames(FrogLakes)[apply(FrogLakes,1,sum) > 0]
#####           

# eDNA and multi-year regional species list
bvenn(list(eDNA = eDNAList, "Multi-year\nobservations" = RegioList), 
      fontsize=10)

# eDNA and audio-visual survey list
bvenn(list(eDNA = eDNAList, "Audio-visual\nsurvey" = SurveyList), 
      fontsize=10)

# Positive controls
# PCE: DNA from 12 species in equal concentrations
C.PCE = grep("sample.P.PCE", names(FrogCounts))

# PCUNE: DNA from 12 species in stepwise doubled concentrations
C.PCUNE = grep("sample.P.PCUNE", names(FrogCounts))

# Species in the positiv control
PosList = sort(c("Phyllomedusa azurea",
                 "Physalaemus centralis",
                 "Hypsiboas geographicus",
                 "Pseudopaludicola mystacalis",
                 "Hypsiboas punctatus",
                 "Leptodactylus fuscus",
                 "Hypsiboas raniceps",
                 "Pseudis limellum",
                 "Leptodactylus podicipinus",
                 "Scinax madeirae",
                 "Dendropsophus leucophyllatus",
                 "Dendropsophus nanus"))

# One of the 12 species (Hypsiboas geographicus) is not present in the PCE results.
IsInPCE = cbind(PosList, PosList %in% rownames(FrogCounts))

# Visualize conncentrations in the PCE. Lines are PCE samples.
palette(colors())
par(mar=c(8,4,1,1))
plot(c(1:12), seq(0.1, (max(FrogCounts[,C.PCE])), 
                  (max(FrogCounts[,C.PCE]))/12), type="n",
     xaxt="n", xlab="", ylab="Read numbers") 
axis(1, at = c(1:12), labels = PosList, las = 2, cex.axis=0.6)
for (i in C.PCE) {
  lines(c(1:12), (FrogCounts[PosList,i]), col=i+10, lwd=2)
}

# The seven PCE samples are positively correlated at R > 0.7.

# Correction factors for PCE (obsolate)
PCEConc = 5 #ng/ul
PCEUNEConc = read.csv(file="PCUNE_conc.csv", header=T, row.names = 1)
# 
# #  the mean abundance of the species in PCE with 1x dilution 
# # in PCUNE should correspond the  5 ng/ul conc: Phyllomedusa azurea
# mean(as.numeric(FrogCounts["Phyllomedusa azurea", C.PCE]))
# 
# # Divide this by all abundances of control species to get the correction factor
# CorrectFactor = mean(as.numeric(FrogCounts["Phyllomedusa azurea", C.PCE])) /
#   apply(FrogCounts[PosList,],1,mean)
# 
# # Corrected PCE abundances
# FrogCountsCorrected = FrogCounts[PosList,] * CorrectFactor
# apply(FrogCountsCorrected, 1, sum)

# Evaluation of non-equimolar concentrations and read numbers
# Correlations with the original DNA concentrations
# PCUNE abundances VS original concentrations

RangeX = range(PCEUNEConc[,"conc"])
RangeY = range(apply(FrogCounts[,C.PCUNE],1,mean))

par(mfrow = c(1,1), mar=c(4,4,3,1))
plot(RangeX, RangeY, type="n", xlab="DNA concentrations", ylab = "Read numbers", 
     log="x", main = "Positiv non-equimolar controls")
for (i in PosList[PosList != "Hypsiboas geographicus"]){
  points(PCEUNEConc[i,"conc"], mean(as.numeric(FrogCounts[i,C.PCUNE])), 
         pch=19)
}

# Correlation of concentrations and read numbers in the PCUNE
ConcMeanRead = data.frame()
for (i in PosList[PosList != "Hypsiboas geographicus"]){
  ConcRead = c(PCEUNEConc[i,"conc"], mean(as.numeric(FrogCounts[i,C.PCUNE])),
               mean(as.numeric(FrogCounts[i,C.PCE])))
  ConcMeanRead = rbind(ConcMeanRead, ConcRead)
}
colnames(ConcMeanRead) = c("conc", "MeanReadPCUNE", "MeanReadPCE")
rownames(ConcMeanRead) = PosList[PosList != "Hypsiboas geographicus"]

cor(ConcMeanRead, method = "pearson")


# # Visualize conncentrations in the PCE. Lines are PCE samples.
# palette(colors())
# par(mar=c(8,4,1,1))
# plot(c(1:11), seq(0.1, (max(FrogCountsCorrected[,C.PCE])), 
#                   (max(FrogCountsCorrected[,C.PCE]))/11), type="n",
#      xaxt="n", xlab="", ylab="Read numbers") 
# axis(1, at = c(1:11), labels = PosList[PosList != "Hypsiboas geographicus"], las = 2, cex.axis=0.6)
# for (i in C.PCE) {
#   lines(c(1:11), (FrogCountsCorrected[PosList != "Hypsiboas geographicus",i]), col=i+10, lwd=2)
# }



#####
# Ecological signal: differences among the lakes VS preservation/extraction methods
# Transpose the frog counts: species in columns, samples in lines
FrogCountsT = t(FrogCounts)

# Filter for lake samples
Samples = SiteMeta$Methods != "Control"
FrogCountsLake = FrogCountsT[Samples,]


# Some samples don't have reads at all, some species were not seen in lakes
summary(apply(FrogCountsLake,1,sum))
summary(apply(FrogCountsLake,2,sum))

# Meaningful species matrix for lakes
FrogCountsNoZero = FrogCountsLake[apply(FrogCountsLake,1,sum) > 0,
                                  apply(FrogCountsLake,2,sum) > 0]

# Corresponding metadata
MetaLakesNoZero = SiteMeta[rownames(SiteMeta) %in% rownames(FrogCountsNoZero),]

# Remove species present in only one sample
FrogCountsNoOne = FrogCountsNoZero[,apply(FrogCountsNoZero,2,function(vec) sum(vec>0)) > 1]

# still species in every sample?
summary(apply(FrogCountsNoOne,1,sum))

# Model-based comparison of the lakes
FrogsMvabund = mvabund(FrogCountsNoOne)

# Multispecies GLM and ANOVA. Reads control for differences in sequencing depth.
m1 = manyglm(FrogsMvabund ~ Reads + Lakes, 
             data=MetaLakesNoZero)

# Increase nBoot=1000, and include p.uni for species-level statistics
m1.anova = anova(m1, nBoot = 1000, p.uni = "adjusted")

# nice anova table
kable(m1.anova$table)

# Visualize the community results
# Model-based ordination with reads as covariates

# Overdispersion parameter
summary(m1$theta)
hist(m1$theta)

set.prior = set.prior <- list(type = c("normal","normal","normal","uniform"),
                              hypparams = c(100, 20, 100, 1))

# Unconstrained ordination
# ord.m.noconst <- boral(FrogCountsNoOne, family = "negative.binomial", num.lv = 2,
#                        prior.control = set.prior,
#                        n.burnin = 10, n.iteration = 100, n.thin = 1)
# 
# par(mfrow = c(2,2))
# plot(ord.m.noconst)

# Constrained ordination with sequencing depth
# Just keep the old 2 LV plot
ord.m.const <- boral(FrogCountsNoOne, X=MetaLakesNoZero$Reads, 
                     prior.control = set.prior,
                     family = "negative.binomial", num.lv = 2, 
                     n.burnin = 10, n.iteration = 100, n.thin = 1)

# Diagnostics
par(mfrow = c(2,2))
plot(ord.m.const)

# Colors for ordination
palette(colors())
LakeCol = gsub("T1", "green", MetaLakesNoZero$Lakes)
LakeCol = gsub("T2", "orange", LakeCol)
LakeCol = gsub("T3", "red", LakeCol)
LakeCol = gsub("T4", "purple", LakeCol)
LakeCol = gsub("T5", "blue", LakeCol)
MyColors = data.frame(cols = c("green","orange","red","purple","blue"),
                      lakes = levels(factor(MetaLakesNoZero$Lakes)))

# Constrained ordination plot
par(mfrow=c(1,1), mar=c(4,4,3,1), oma=c(1,1,0,0))
ordifrog= ordiplot(ord.m.const$lv.median, type = "none", cex =0.5 
                   ,display = "sites", xlim = c(-4,3),
                   main="Constrained LV ordination")
points(ordifrog,"sites", pch=20 ,col=LakeCol)
legend(2, 5, paste(c("Small Croco", "Vastus", "Wetland basin", "Lacha Susanna", 
                     "Wetland Centro")), fill=as.vector(MyColors$cols), 
       border="white", bty="n", cex=0.7)
ordiellipse(ordifrog, factor(MetaLakesNoZero$Lakes), cex=.5, 
            draw="polygon", col="green",
            alpha=100,kind="se",conf=0.95, 
            show.groups=(c("T1")))
ordiellipse(ordifrog, factor(MetaLakesNoZero$Lakes), cex=.5, 
            draw="polygon", col="orange",
            alpha=100,kind="se",conf=0.95, 
            show.groups=(c("T2")))
ordiellipse(ordifrog, factor(MetaLakesNoZero$Lakes), cex=.5, 
            draw="polygon", col="red",
            alpha=100,kind="se",conf=0.95, 
            show.groups=(c("T3")))
ordiellipse(ordifrog, factor(MetaLakesNoZero$Lakes), cex=.5, 
            draw="polygon", col="purple",
            alpha=100,kind="se",conf=0.95, 
            show.groups=(c("T4")))
ordiellipse(ordifrog, factor(MetaLakesNoZero$Lakes), cex=.5, 
            draw="polygon", col="blue",
            alpha=100,kind="se",conf=0.95, 
            show.groups=(c("T5")))

# Unconstrained ordination plot
# plot(ord.m.noconst$lv.median, col=LakeCol, 
#      pch=19, main="Unconstrained LV ordination", las=1, xlab = "LV1", ylab = "LV2")
# legend(-6.5, 2, c("Small Croco", "Vastus", "Wetland basin", "Lacha Susanna", 
#                   "Wetland Centro"), fill=as.vector(MyColors$cols), 
#        border="white", bty="n")
# for (i in levels(factor(MetaLakesNoZero$Lakes))) {
#   ellipse(center = c(mean(ord.m.noconst$lv.median[MetaLakesNoZero$Lakes == i,"LV1"]), 
#                      mean(ord.m.noconst$lv.median[MetaLakesNoZero$Lakes == i,"LV2"])),
#           shape = cov(ord.m.noconst$lv.median[MetaLakesNoZero$Lakes == i,]),
#           radius = 0.95*mean(std.error(ord.m.noconst$lv.median[MetaLakesNoZero$Lakes == i,])), 
#           add=T, col=as.vector(MyColors$cols[MyColors$lakes == i]))
# }


# Distribution of read abundances
hist(apply(FrogCountsT,2,sum), nclass=20, col="grey", 
     main="Distribution of read abundances",
     xlab="Read counts per species", ylab="Frequency")

# Simple NMDS plot
# LakeNMDS = metaMDS(FrogCountsT[SamplesNoZero & Samples,])
# 
# par(mar=c(4,4,1,1))
# plot(LakeNMDS$points, type="n", xlab="NMDS1", ylab="NMDS2")
# ordispider(LakeNMDS, SiteMeta[SamplesNoZero & Samples,"Lakes"], col="grey")
# points(LakeNMDS, pch=20)
# mylegend = legend(-3.3, 3, c("Small Croco", "Vastus", "Wetland basin", "Lacha Susanna", "Wetland Centro"), 
#                   fill=c("orange","green","purple","red","blue"), border="white", bty="n")
# ordiellipse(LakeNMDS, SiteMeta[SamplesNoZero & Samples,"Lakes"],cex=.5, draw="polygon", 
#             col=c("orange"), alpha=170,kind="se",conf=0.95, show.groups=(c("T1")))
# ordiellipse(LakeNMDS, SiteMeta[SamplesNoZero & Samples,"Lakes"],cex=.5, draw="polygon", col=c("green"),
#             alpha=170,kind="se",conf=0.95, show.groups=(c("T2")))
# ordiellipse(LakeNMDS, SiteMeta[SamplesNoZero & Samples,"Lakes"],cex=.5, draw="polygon", col=c("purple"),
#             alpha=170,kind="se",conf=0.95, show.groups=(c("T3")))
# ordiellipse(LakeNMDS, SiteMeta[SamplesNoZero & Samples,"Lakes"],cex=.5, draw="polygon", col=c("red"),
#             alpha=170,kind="se",conf=0.95, show.groups=(c("T4")))
# ordiellipse(LakeNMDS, SiteMeta[SamplesNoZero & Samples,"Lakes"],cex=.5, draw="polygon", col=c("blue"),
#             alpha=170,kind="se",conf=0.95, show.groups=(c("T5")))
# 
# 
# 
# 



# Preservation-extraction method evaluation with species occupancy models
# Script from Thiery comes here.

# Format the data for SOM
# !!! Should be updated.
# data frame of non-control samples
# positive contols: 
# C.PCE = grep("sample.P.PCE", names(SpeCounts))
# C.PCUNE = grep("sample.P.PCUNE", names(SpeCounts))
# 
# # negative controls: 
# C.PNC = grep("sample.P.NC", names(SpeCounts))
# C.NTC = grep("sample.NTC", names(SpeCounts))
# C.NC = grep("sample.NC", names(SpeCounts))
# C.MPX = grep("sample.MPX", names(SpeCounts))
# 
# Controls = c(C.PCE,C.PCUNE,C.PNC,C.NTC,C.NC,C.MPX)
# 
# # reads from lake observations
# AbundLakes = SpeCounts[-Controls]
# 
# # Transform to presence-absences for species occupancy models
# AbundSOM = AbundLakes
# AbundSOM[AbundSOM > 0] <- 1
# 
# # write.csv(file = "SOM_data.csv", AbundSOM)

# Results from Thierry
# SOM1: all species modelled together - for a general comparison of the 
# three preservation / extraction approaches.
# psi = Pr(occupancy)
# p = Pr(detection)
# fp = Pr(flase pos.)
# b = Pr(ambiguous detection occurs) == only 1 PCR positive / 4
SOM1 = read.csv(file = "SOM/SOM_results1.csv", row.names = 1, header=T)

# Plot method effects
# Detection probabilities per method
pDetect = grep('^p\\(meth.\\)', rownames(SOM1))
par(mfrow=c(2,1), mar=c(4,10,1,1))
plot(SOM1[pDetect,"estimate"], c(3:1), xlim=c(min(SOM1$X2.5..CI[pDetect]),
                                              max(SOM1$X97.5..CI[pDetect])), 
     yaxt = "n", xlab = "Detection probability", ylab="", 
     ylim=c(0.5,3.5), pch=19)
segments(SOM1[pDetect,2], c(3:1), SOM1[pDetect,3], c(3:1))
axis(2, at=c(3:1), label=c("GFF (2 um, in liquid: A)", 
                           "GFF (2 um, dried: B)",
                           "Nylon (0.2 um, dried: C)"), las=1)

# Plot method false positives
# False positives per method
pFPos = grep('^fp\\(meth.\\)', rownames(SOM1))
plot(SOM1[pFPos,"estimate"], c(3:1), xlim=c(min(SOM1$X2.5..CI[pFPos]),
                                            max(SOM1$X97.5..CI[pFPos])), 
     yaxt = "n", xlab = "False positives", ylab="", 
     ylim=c(0.5,3.5), pch=19)
segments(SOM1[pFPos,2], c(3:1), SOM1[pFPos,3], c(3:1))
axis(2, at=c(3:1), label=c("GFF (2 um, in liquid: A)", 
                           "GFF (2 um, dried: B)",
                           "Nylon (0.2 um, dried: C)"), las=1)

# SOM2: species modelled separately - for species-specific methodological differences
# Detection probabilities for each species per method
SOM2.p = read.csv(file = "SOM/SOM2_detection_prob.csv", row.names = 1, header=T)

# False positives, each species per method
SOM2.fp = SOM2.p = read.csv(file = "SOM/SOM2_false_pos.csv", row.names = 1, header=T)

# Coefficient plots
par(mfrow=c(1,3), mar=c(4.1,1,0.5,0.5))
plot(rep(0,nrow(SOM2.p)), c(nrow(SOM2.p):1), type="n", 
     bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
plot(SOM2.p[,"pA_estimate"], c(nrow(SOM2.p):1)-0.2, 
     xlim=c(0,0.5), yaxt="n",
     xlab="Detection probability", ylab="", pch=19, cex=0.5, 
     ylim=c(nrow(SOM2.p),1))
points(SOM2.p[,"pB_estimate"], c(nrow(SOM2.p):1), pch=17, cex=0.5)
points(SOM2.p[,"pC_estimate"], c(nrow(SOM2.p):1)+0.2, pch=15, cex=0.5)
segments(SOM2.p$pA_2.5..CI, c(nrow(SOM2.p):1)-0.2, 
         SOM2.p$pA_97.5..CI, c(nrow(SOM2.p):1)-0.2)
segments(SOM2.p$pB_2.5..CI, c(nrow(SOM2.p):1), 
         SOM2.p$pB_97.5..CI, c(nrow(SOM2.p):1))
segments(SOM2.p$pC_2.5..CI, c(nrow(SOM2.p):1)+0.2, 
         SOM2.p$pC_97.5..CI, c(nrow(SOM2.p):1)+0.2)
for(i in 1:nrow(SOM2.p)){
  lines(c(-0,1), c(i,i)-0.5, lty="dotted")
}
axis(2, at=c(nrow(SOM2.p):1), label=rownames(SOM2.p), las=1, cex.axis=0.9,
     tick=F)
plot(SOM2.fp[,"pA_estimate"], c(nrow(SOM2.p):1)-0.2, 
     xlim=c(0,0.5), yaxt="n",
     xlab="False positives", ylab="", pch=19, cex=0.5, 
     ylim=c(nrow(SOM2.p),1))
points(SOM2.fp[,"pB_estimate"], c(nrow(SOM2.fp):1), pch=17, cex=0.5)
points(SOM2.fp[,"pC_estimate"], c(nrow(SOM2.fp):1)+0.2, pch=15, cex=0.5)
segments(SOM2.fp$pA_2.5..CI, c(nrow(SOM2.fp):1)-0.2, 
         SOM2.fp$pA_97.5..CI, c(nrow(SOM2.fp):1)-0.2)
segments(SOM2.fp$pB_2.5..CI, c(nrow(SOM2.fp):1), 
         SOM2.fp$pB_97.5..CI, c(nrow(SOM2.fp):1))
segments(SOM2.fp$pC_2.5..CI, c(nrow(SOM2.fp):1)+0.2, 
         SOM2.fp$pC_97.5..CI, c(nrow(SOM2.fp):1)+0.2)
for(i in 1:nrow(SOM2.p)){
  lines(c(-0,1), c(i,i)-0.5, lty="dotted")
}





# Amphibians
Amphi = c("Dendropsophus leucophyllatus","Dendropsophus melanargyreus",
          "Dendropsophus minutus","Dendropsophus nanus","Dermatonotus muelleri",
          "Elachistocleis sp.","Eupemphix nattereri","Hylinae","Hyloidea",
          "Hypsiboas punctatus","Hypsiboas raniceps","Leptodactylus",
          "Leptodactylus chaquensis","Leptodactylus elenae","Leptodactylus fuscus",
          "Leptodactylus podicipinus","Leptodactylus syphax","Leptodactylus vastus",
          "Osteocephalus","Osteocephalus taurinus","Phyllomedusa azurea",
          "Phyllomedusa boliviana","Physalaemus albonotatus","Physalaemus centralis",
          "Physalaemus cuvieri","Pseudis limellum","Pseudis paradoxa",
          "Pseudopaludicola mystacalis","Rhinella schneideri","Scinax fuscomarginatus",
          "Scinax fuscovarius","Scinax nasicus","Scinax ruber","Scinax sp. FB-2014a",
          "Bufonidae","Elachistocleis","Gastrophryninae","Hylidae",
          "Leptodactylus latinasus","Melanophryniscus","Pseudis","Pseudis laevis","Scinax")
nrow(AbundSOM[Amphi,])

InputSOM = AbundSOM[Amphi,grep("sample.T", names(AbundSOM))]

# create SOM alternative input for export
cbind(rownames(InputSOM)[1], t(InputSOM[1,]))

ArrangedSOM = matrix(ncol = 2)
for (i in 1:nrow(InputSOM)) {
  ArrangedSOM = rbind(ArrangedSOM, cbind(rownames(InputSOM)[i], t(InputSOM[i,])))
}

# This writes out the file that can be edited to get the input format required by Thierry
# write.csv(file="input_som.csv", ArrangedSOM)

# Clarify vague identifications
# Osteocephalus
grep("Osteocephalus", MetaHead$sci_name)
osteo = MetaHead[grep("Osteocephalus", MetaHead$sci_name),]

# Elachistocleis
elach = MetaHead[grep("Elachistocleis", MetaHead$sci_name),]

# Pseudis
pseud = MetaHead[grep("Pseudis", MetaHead$sci_name),]

write.csv(file = "species-to-clear.csv", rbind(osteo,elach,pseud))

# tree checked in seaview, I would believe all variants.

# 















## Hominidae sequences
sum(apply(AbundHead, 1, sum)[MetaHead$family == "Hominidae"], na.rm=T)


summary(MetaAll$sci_name)


# Read distributions
hist(apply(AllData[,8:679], 2, sum), nclass = 20)
sum(AllData[,8:679])

# Calculate how often a sequence is head, internal or singleton
AllData = cbind(AllData, 
                no_sing = apply(AllData[,680:1351], 1, function(x) length(which(x == "s"))),
                no_inte = apply(AllData[,680:1351], 1, function(x) length(which(x == "i"))),
                no_head = apply(AllData[,680:1351], 1, function(x) length(which(x == "h"))))

HasHead = AllData$no_head != 0

# Keep sequence variants that were at least once observed as head 
DataWitHead = AllData[HasHead,]
sum(DataWitHead[,8:679])
hist(apply(DataWitHead[,8:679], 2, sum), nclass = 20)
summary(apply(DataWitHead[,8:679], 2, sum))

# retain only 16S samples
All_16S = grep("16S.$", names(DataWitHead))
Heads16S = cbind(DataWitHead[,1:7], DataWitHead[,All_16S], DataWitHead[,1352:1362])
sum(Heads16S[,8:343])

write.csv(file="16S_R_filtered.csv", Heads16S)





TotPresent = apply(LakeReads,2,function(vec) sum(vec>0))
hist(TotPresent)


# POPs and elements
Meta = read.csv(file="../../../Data/Stechlin_pre-analysis/Stechlin_POP_elements.csv",
                header=T, row.names = 1)

