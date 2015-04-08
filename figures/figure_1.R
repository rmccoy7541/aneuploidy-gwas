#################################################
# File: figure_1.R
#################################################
# Author: Rajiv McCoy
# This script was used to generate Fig. 1B,C.
# The purpose of these figures is to describe
# characteristics of aneuploidies affecting  
# paternal chromosome copies, which we propose
# are predominantly mitotic in origin.
#################################################

# load packages
library(ggplot2)
library(gridExtra)

source("~/Desktop/aneuploidy_functions.R")

# the aneuploidy calls were published along with the paper

data <- read.table("~/Desktop/aneuploidy_calls.csv", sep = ",", header=T) # import the data

data_filtered <- filterData(data)

data_blastomere <- selectSampleType(data_filtered, blastomere)

data_blastomere <- callPloidy(data_blastomere)

# Read in de-identified genotype codes for single, common 
# associated SNP. This is necessary to limit figures to only 
# the set of cases used in the GWAS (i.e., no repeat cases).

gt <- read.table("/sandbox/rs2305957_unrelated_geno05_mind05_f.gt");
names(gt) <- c("case", "sample_id", "genotype")

data_gt <- merge(data_blastomere, gt, "case")

# Count different forms of aneuploidy

paternal_error_frame <- data.frame(matrix(ncol = 23, nrow = nrow(data_gt)))
for (i in 7:29) {
	new<-((data_gt[,i] == "H100" | data_gt[,i] == "H120" | data_gt[,i] == "H020") & data_gt[,i + 46]==0 & data_gt[,i + 69]==0)
	paternal_error_frame[,i - 6]<-new
}
maternal_error_frame<-data.frame(matrix(ncol = 23, nrow = nrow(data_gt)))
for (i in 7:29) {
	new<-((data_gt[,i] == "H010" | data_gt[,i] == "H210" | data_gt[,i] == "H200") & data_gt[,i + 69]==0)
	maternal_error_frame[,i - 6]<-new
}
maternal_monosomy_frame<-data.frame(matrix(ncol = 23, nrow = nrow(data_gt)))
for (i in 7:29) {
	new<-((data_gt[,i] == "H010" | data_gt[,i] == "H001") & data_gt[,i + 69]==0)
	maternal_monosomy_frame[,i - 6]<-new
}
paternal_monosomy_frame<-data.frame(matrix(ncol = 23, nrow = nrow(data_gt)))
for (i in 7:29) {
	new<-(data_gt[,i] == "H100" & data_gt[,i + 69]==0)
	paternal_monosomy_frame[,i - 6]<-new
}
paternal_trisomy_frame<-data.frame(matrix(ncol = 23, nrow = nrow(data_gt)))
for (i in 7:29) {
	new<-(data_gt[,i] == "H120" & data_gt[,i + 69]==0)
	paternal_trisomy_frame[,i - 6]<-new
}
nullisomy_frame<-data.frame(matrix(ncol = 23, nrow = nrow(data_gt)))
for (i in 7:29) {
	new<-(data_gt[,i] == "H000" & data_gt[,i + 69]==0)
	nullisomy_frame[,i - 6]<-new
}

sum(apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
sum(apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)
sum(apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 )

sum(apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)

sum(apply(maternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)

sum(apply(paternal_monosomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)]==TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)]==TRUE)) > 0 & apply(nullisomy_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0)


#################################################

### generate figures ###

paternal_hist <- hist(data_gt[apply(paternal_error_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0,]$chroms_affected, breaks=0:23)
no_paternal_hist <- hist(data_gt[data_gt$chroms_affected > 0 & apply(paternal_error_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) == 0,]$chroms_affected, breaks=0:23)

results <- data.frame(rbind(cbind(paternal_hist$density, "Paternal chrom. affected"), cbind(no_paternal_hist$density, "No paternal chrom. affected")))
results$X1 <- as.numeric(as.character(results$X1))
a <- ggplot(data = results, aes(y = X1, fill = factor(X2), x = rep(1:23, 2))) + geom_bar(stat = "identity", position = "dodge") + xlab("No. of aneuploid chroms.") + labs(fill = "") + ylab('Density') + theme(legend.justification = c(0,0), legend.position = c(.05, 0.65))



data_gt$paternal <- FALSE
data_gt[apply(paternal_error_frame[,1:23], 1, function(x) sum(x[!is.na(x)] == TRUE)) > 0,]$paternal <- TRUE

aneuploidyByAge(data_gt[data_gt$paternal == TRUE,], "Paternal chrom. affected")
aneuploidyByAge(data_gt[data_gt$paternal == TRUE,], "No paternal chrom. affected")


#aggregate over chromosomes stratified by maternal age
results_all_chroms<- data.frame(matrix(ncol = 4))
minAge <- round(min(data_gt$maternal_age[!is.na(data_gt$maternal_age)]))
maxAge <- round(max(data_gt$maternal_age[!is.na(data_gt$maternal_age)]))
for (k in minAge:maxAge) {	
	age_subset <- data_gt[round(data_gt$maternal_age) == k & !is.na(round(data_gt$maternal_age)),]
	paternal_prop <- sum(age_subset$paternal == TRUE) / nrow(age_subset)
	paternal_se <- sqrt((paternal_prop * (1 - paternal_prop)) / nrow(age_subset))
	paternal <- c(paternal_prop, paternal_se, k, "Paternal chrom. affected")
	no_paternal_prop <- sum(age_subset$paternal == FALSE & age_subset$ploidy == FALSE) / nrow(age_subset)
	no_paternal_se <- sqrt((no_paternal_prop * (1 - no_paternal_prop)) / nrow(age_subset))
	no_paternal <- c(no_paternal_prop, no_paternal_se, k, "No paternal chrom. affected")
	results_all_chroms<-rbind(results_all_chroms, paternal, no_paternal)
}
results_all_chroms<-results_all_chroms[-1,]
results_all_chroms$X1<-as.numeric(results_all_chroms$X1)
results_all_chroms$X2<-as.numeric(results_all_chroms$X2)
results_all_chroms$X3<-as.numeric(results_all_chroms$X3)
names(results_all_chroms) <- c("prop", "se", "age", "category")

### generate figures ### 

limits <- aes(ymax = (prop + se), ymin = (prop - se))

b <-ggplot(data = results_all_chroms, aes(x = age, y = prop, col = factor(category))) + geom_line() + geom_point() + theme(legend.title = element_blank()) + xlab("Maternal age") + ylab("Prop. aneuploid blastomeres") + theme(legend.justification = c(0, 0), legend.position=c(.0, 0.7)) + geom_errorbar(limits, width = 0.5)

# put them in a single multi-panel figure

a<-ggplot_gtable(ggplot_build(a))
b<-ggplot_gtable(ggplot_build(b))
b$widths <- a$widths
grid.newpage()
grid.arrange(a, b, nrow = 1)
