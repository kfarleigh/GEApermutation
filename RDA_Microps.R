

library(psych)
library(vegan)
library(trio)
library(LEA)
library(dplyr)
library(vcfR)
library(caret)
library(robust)
library(SoDA)
library(adespatial)
library(gtools)
library(qvalue)

###########################
##### !!! Read ME !!! #####
###########################

# This script was adapted from Forester et al., 2018 
# https://popgen.nescent.org/2018-03-27_RDA_GEA.html
# For rda analysis the genotype matrix has to been complete (No NA's) and values can be 0, 1, or 2.
# This is the same format as in LEA
# So what we'll do is take the imputed lfmm file and use that 
# See the LEA tutorial by Olivier Francois or Frichot et al., 2015 for more information  regarding formatting/imputation

# For rda, we will use an enviromental csv with populations and coordinates for downstream visualization purposes

Gtypes <- 'Microps_Ref_C.lfmm'
Env <- 'Microps_EnvData.csv'
Outputprefix <- 'Microps'
###############################################
##### Bring in files and check formatting #####
###############################################

Genotypes<-read.lfmm(Gtypes)
EnvVars<-read.csv(Env)

sum(is.na(Genotypes))
9 %in% Genotypes
# or
any(Genotypes == 9)
# checks for any NA's or 9's in your dataframe

str(EnvVars)
EnvVars$Samples<-as.character(EnvVars$Samples)
# check the structure of your environmental file and change Samples to characters

Samples<-EnvVars$Samples
row.names(Genotypes)<- Samples
row.names(Genotypes)
identical(rownames(Genotypes), EnvVars$Samples)
# makes sure that your samples match each other in your genotypes and environmental variables


################################################################
##### Check for Correlation between Environmental Variables#####
################################################################

# Find variable set that retains the most variables
# First we identify which variables we should remove 
Cor <- findCorrelation(cor(EnvVars[,6:25]), cutoff = 0.7)
SelectedVars <- EnvVars[,6:25]
SelectedVars <- SelectedVars[-Cor]
# with RDA we have to check for correlation between environmental predicts as to not bias our results

Correlation<-cor.plot(SelectedVars, xlas =  2, cex.axis = 0.75, MAR = 8)
# look at correlation between environmental variables, when using 10 or less you can use the paris plot
# otherwise use the correlation plot as it is easier to visualize with more variables


colnames(SelectedVars) <- c("NDVI","MDR","IsoT","TS", "MTWetQ", "MTWarmQ","PWarmQ","PCQ")

write.csv(SelectedVars, file = paste(Outputprefix, "_SelectedVars", ".csv", sep = ""), row.names = FALSE)
# write it out as a file

###########################
##### Model Selection #####
###########################
mod0 <- rda(Genotypes ~ 1, data = SelectedVars, scale = TRUE)
mod1 <- rda(Genotypes ~ ., data = SelectedVars, scale = TRUE)

mod <- ordistep(mod0, scope = formula(mod1))

# Which variables should you include?
mod 

########################
##### Run your RDA #####
########################

rda<- rda(Genotypes ~ ., data = SelectedVars, scale = TRUE)
rda
# runs the rda and looks at it

Rsquared<-RsquareAdj(rda)
Rsquared
# look at the rsquared and adjusted rsquared values

summary(eigenvals(rda, model = "constrained"))
# provides a summary of the eigenvalues for each axis

screeplot(rda, main = "Eigenvalues from RDA", col = "navy")
# plots a scree plot of the eigenvalues

Fullsig <- anova.cca(rda, parallel=getOption("mc.cores"))
Fullsig
# test the significance of the rda
# displays results

Axissig <- anova.cca(rda, by="axis", parallel=getOption("mc.cores"))
Axissig
# find which constrained axes are significant
# this can take awhile 

vif.cca(rda)
# checks the variance inflation factors
# the lower the better


#############################
# lets make some cool plots

EnvVars$Habitat <- as.factor(EnvVars$Habitat)
Pop<- as.factor(EnvVars$Pops)
bg <- c("#FF0000", "#222568")
# pull out populations and get some colors

###############
# axes 1 & 2
plot(rda, type="n", scaling=3)
points(rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3) #displays SNPs
points(rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[EnvVars$Habitat]) #displays individuals
text(rda, scaling=3, display="bp", col="356689", cex=1.25)  #displays predictors
legend("bottomright", legend=levels(EnvVars$Habitat), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg) #adds legend

###############
# axes 1 & 3

plot(rda, type="n", scaling=3, choices = c(1,3))
points(rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices = c(1,3)) #displays SNPs
points(rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[EnvVars$Habitat], choices = c(1,3)) #displays individuals
text(rda, scaling=3, display="bp", col="#0868ac", cex=1, choices = c(1,3))  #displays predictors
legend("bottomright", legend=levels(EnvVars$Habitat), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg) #adds legend

###########################
##### Capblanq method #####
###########################

env = SelectedVars

# How many significant axes
sigaxes <- 5


##################################################
##### Call candidates using Capblancq Method #####
##################################################

# K is the number of RDA axes you want to include (include significant ones)
rdadapt<-function(rda,K)
{
  loadings<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(loadings, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

res_rdadapt<-rdadapt(rda, 5)
res_rdadapt$BH <- p.adjust(res_rdadapt$p.values, method = 'BH')
res_rdadapt$Bon <- p.adjust(res_rdadapt$p.values, method = 'bonferroni')
outliers <- res_rdadapt[which(res_rdadapt$BH < 0.1),]

# Correlation env ~ PC or env ~ RDA
r2_RDA<-NULL
for (i in 1:5)
{
  r2_RDA<-c(r2_RDA,summary(lm(env[,i] ~ rda$CCA$u[,i]))$r.squared)
}
# How correlated is each predictor with each RDA axis?
# Nrow is the number of predictors, ncol is the number of retained RDA axes
# i represents predictors, j = RDA axes
var_RDA<-matrix(nrow = 5, ncol = 5)
for (i in 1:5)
{for (j in 1:5)
{
  var_RDA[i,j]<-summary(lm(env[,i] ~ rda$CCA$u[,j]))$r.squared
}
}
row.names(var_RDA) <- colnames(SelectedVars)
colnames(var_RDA) <- colnames(rda$CCA$u)

##################################
##### Investigate Candidates #####
##################################

colnames(Genotypes) = SNPs$V1

foo <- matrix(nrow=ncol(var_RDA), ncol=sigaxes)
rownames(foo) <- colnames(var_RDA)1
colnames(foo) <- c("NDVI","MTWetQ","MTWarmQ","PWarmQ","PColdQ")
for (i in seq(1:ncol(var_RDA))) {
  nam <- SNPs[i,]
  snp.gen <- Genotypes[,nam]
  foo[i,] <- apply(SelectedVars,2,function(x) cor(x,snp.gen))
}

# Get the environmental data for each of the outliers
Cands <-  data.frame(foo[as.numeric(rownames(outliers)),])
ncolCand <- ncol(Cands)
CandVar <- (ncol(Cands) + 1)
CandVarCor <- (ncol(Cands) + 2)
for (i in rownames(Cands)) {
  bar <- Cands[i,]
  Cands[i,CandVar] <- names(which.max(abs(bar[1:ncolCand]))) # gives the variable
  Cands[i,CandVarCor] <- max(abs(bar[1:ncolCand]))              # gives the correlation
}

colnames(Cands)[6:7] <- c('Variable', 'Correlation')
remove(nam)
remove(snp.gen)
remove(foo)
remove(bar)

############################
##### Permutation Test #####
############################
# This is a test to determine if an adaptive SNP is overrepresented by a specific habitat type
# Use the identify adaptive snps script to determine which allele is adaptive at each snp first
# Then use that information as input here 

# Read in files from identifying the adaptive SNP
# If the major allele was adaptive then this means that the reference allele is of interest
# It thats the case then we will recode those genotypes, warning this will mess up your genotypes 
# So make it a seperate object as done here and do not use it for anything else

Majoradapt <- read.delim('Microps_Majoradapt.txt')
Minoradapt <- read.delim('Microps_Minoradapt.txt')

# Get the names of the candidates
candidates <- rownames(Cands)

# Get the genotypes at each candidate
# Rename the rows to show the corresponding sample
cand_gtypes <- data.frame(Genotypes[,candidates])
cand_gtypes$Samples <- rownames(cand_gtypes)

# Are the samples in a blackbrush (1) or atriplex (0) locality
cand_gtypes$Habitat <- EnvVars$Habitat
cand_gtypes$Habitat[cand_gtypes$Habitat == 'Atriplex'] <- 0
cand_gtypes$Habitat[cand_gtypes$Habitat == 'BlackBrush'] <- 1

# Change genotype coding for the adaptive major alleles
Major_adaptive <- which(colnames(cand_gtypes) %in% Majoradapt$x)

# Get the genotypes to change and remove them from the original data frame
Maj_adapt <- cand_gtypes[,Major_adaptive]
Maj_adapt[Maj_adapt == 2] <- 9
Maj_adapt[Maj_adapt == 0] <- 2
Maj_adapt[Maj_adapt == 9] <- 0

cand_gtypes <- cand_gtypes[,-Major_adaptive]

# Combine the altered data frame with the original ones       ### Warning, your data is now altered (cand_gtypes) ###
cand_gtypes <- cbind(Maj_adapt, cand_gtypes)
cand_names <- cand_gtypes[,-c((ncol(cand_gtypes)- 2),(ncol(cand_gtypes)-1))]

# For each candidate loci which individuals possess the adaptive state?
ncands <- ncol(cand_names)
Loc_ind <- list()
for (i in 1:ncands) {
  Inds <- which(cand_gtypes[,i] == 2)
  Loc_ind[[i]] <- cand_gtypes$Habitat[Inds]
  
}

names(Loc_ind) <- colnames(cand_names)

# Calculate our empirical statistic (proportion of individuals from either habitat type)
# Atriplex 
Atriplex_prop <- list()
for (i in 1:ncands) {
  Loc <- as.numeric(Loc_ind[[i]])
  obs_prop <- length(which(Loc == 0))/length(Loc)
  Atriplex_prop[[i]] <- obs_prop
}
names(Atriplex_prop) <- colnames(cand_names)

# Blackbrush 
Blackbrush_prop <- list()
for (i in 1:ncands) {
  Loc <- as.numeric(Loc_ind[[i]])
  obs_prop <- length(which(Loc == 1))/length(Loc)
  Blackbrush_prop[[i]] <- obs_prop
}
names(Blackbrush_prop) <- colnames(cand_names)


# Calculate our permuted statistic
# First we need to shuffle the tags 
permutations <- 1000 # number of permutations

tags_permutations <- list()

for (i in 1:permutations) {
  tags_permutations[[i]] <- permute(as.numeric(cand_gtypes$Habitat))
  }


# Retain Unique permutations
tags_permutations <- unique(tags_permutations)


# Set up the progress bar 
cat(paste("Running Permutations"),
    fill=1); flush.console()
# Set progress bar 
pb <- txtProgressBar(min = 0, max = (ncol(cand_gtypes)-2), initial = 0, style = 3)

# Run the permutations, we use a progress bar since this can take a while  
Perm_Atri_prop <- list()
Atri_pval <- list()
Perm_prop <- NULL
for (j in 1:(ncol(cand_gtypes)-2)) {
  for (i in 1:permutations) {
  PermSNPi <- cand_gtypes[, c(j, ncol(cand_gtypes))]
  PermSNPi$Habitat <- tags_permutations[[i]]
  PermSNP_adapt <- PermSNPi[which(PermSNPi[,1] == 2),]
  Perm_prop <- c(Perm_prop, length(which(PermSNP_adapt[,2] == 0))/length(PermSNP_adapt[,2]))
  }
P <- mean(Perm_prop > Atriplex_prop[[j]])
Perm_Atri_prop[[j]] <- Perm_prop
Atri_pval[[j]] <- P
Perm_prop <- NULL
P <- NULL
setTxtProgressBar(pb,j)
}

Atriplex_pvals <- data.frame(do.call(rbind, Atri_pval))
Atriplex_pvals[Atriplex_pvals == 0] <- (1/permutations)
rownames(Atriplex_pvals) <- colnames(cand_names)
colnames(Atriplex_pvals) <- 'P-value'

Atriplex_pvals$BH <- p.adjust(Atriplex_pvals$`P-value`, method = 'BH')
Atriplex_pvals$Bon <- p.adjust(Atriplex_pvals$`P-value`, method = 'bonferroni')
Atriplex_overreps <- Atriplex_pvals[which(Atriplex_pvals$BH < 0.1),]

### Blackbrush 
Perm_Blackbrush_prop <- list()
Blackbrush_pval <- list()
for (j in 1:(ncol(cand_gtypes)-2)) {
  for (i in 1:permutations) {
    PermSNPi <- cand_gtypes[, c(j, ncol(cand_gtypes))]
    PermSNPi$Habitat <- tags_permutations[[i]]
    PermSNP_adapt <- PermSNPi[which(PermSNPi[,1] == 2),]
    Perm_prop <- c(Perm_prop, length(which(PermSNP_adapt[,2] == 1))/length(PermSNP_adapt[,2]))
  }
  P <- mean(Perm_prop > Blackbrush_prop[[j]])
  Perm_Blackbrush_prop[[j]] <- Perm_prop
  Blackbrush_pval[[j]] <- P
  Perm_prop <- NULL
  P <- NULL
  setTxtProgressBar(pb,j)
}

Blackbrush_pvals <- data.frame(do.call(rbind, Blackbrush_pval))
Blackbrush_pvals[Blackbrush_pvals == 0] <- (1/permutations)
rownames(Blackbrush_pvals) <- colnames(cand_names)
colnames(Blackbrush_pvals) <- 'P-value'

Blackbrush_pvals$BH <- p.adjust(Blackbrush_pvals$`P-value`, method = 'BH')
Blackbrush_pvals$Bon <- p.adjust(Blackbrush_pvals$`P-value`, method = 'bonferroni')
Blackbrush_overreps <- Blackbrush_pvals[which(Blackbrush_pvals$BH < 0.1),]

write.csv(Atriplex_overreps, file = 'Overrepresented_Atriplex.csv')
write.csv(Blackbrush_overreps, file = 'Overrepresented_Blackbrush.csv')



# !!!!!!!!!!!!!!!!!!!!!!! Your Done !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# for more information go to https://popgen.nescent.org/2018-03-27_RDA_GEA.html
