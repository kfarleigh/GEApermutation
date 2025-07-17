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