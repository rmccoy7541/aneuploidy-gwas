#################################################
# File: figure_3.R
#################################################
# Author: Rajiv McCoy
# The purpose of this figure is to describe 
# various aspects of the effect of genotype at 
# the associated locus on various phenotypes 
# related to aneuploidy. These figures require
# genotype data to produce.
#################################################

# load packages
library(ggplot2)
library(gridExtra)
library(gtable)

source("~/Desktop/aneuploidy_functions.R")

# the aneuploidy calls were published along with the paper

data <- read.table("~/Desktop/aneuploidy_calls.csv", sep = ",", header=T) # import the data

# don't do any QC filtering, since we only care about number of embryos submitted, not their ploidy calls
# data_filtered <- filterData(data)  

data_te <- selectSampleType(data, TE)

data_te <- callPloidy(data_te)

#################################################

# read in genotype data at the associated locus; combine "discovery" and "validation" sets from GWAS
gt <- read.table("/sandbox/rs2305957_f.gt"); gt2 <- read.table("/sandbox/rs2305957_validate_f.gt"); gt <- rbind(gt, gt2)
names(gt) <- c("case", "sample_id", "genotype")
data_te <- merge(data_te, gt, "case")

means <- data.frame(table(data_te$genotype) / table(data_te[!duplicated(data_te$case),]$genotype))

fam_counts <- data.frame(table(data_te$case))
names(fam_counts) <- c("case", "counts")
results <- merge(fam_counts, gt, "case")
results$gt <- NA
results[results$genotype == "GG",]$gt <- 0
results[results$genotype == "AG",]$gt <- 1
results[results$genotype == "AA",]$gt <- 2
summary(lm(data = results, counts ~ gt))
summary(glm(data = results, counts ~ gt, family = quasipoisson()))
maternal_age <- data.frame(aggregate(data_te$maternal_age ~ data_te$case, FUN = mean))
names(maternal_age) <- c("case", "maternal_age")
results<- merge(results, maternal_age, "case")

# fit a Poisson linear model testing for association between number of day 5 embryos submitted and maternal genotype; include maternal age as a covariate
summary(glm(data = results, counts ~ maternal_age + I(maternal_age ^ 2) + gt, family = quasipoisson()))

#################################################

means <- c(mean(results[results$genotype == "GG",]$counts), mean(results[results$genotype == "AG",]$counts), mean(results[results$genotype == "AA",]$counts))

se <- c(std(results[results$genotype == "GG",]$counts), std(results[results$genotype == "AG",]$counts), std(results[results$genotype == "AA",]$counts))

bars <- data.frame(means, se, c("GG", "AG", "AA"))
names(bars) <- c("counts", "se", "genotype")

limits <- aes(ymax = (counts + se), ymin = (counts - se))
c <- ggplot(data = bars, aes(x = genotype, y = means, fill = genotype)) + geom_bar(stat = "identity") + geom_errorbar(limits, alpha = 1, width = 0.25) + ylab("TE samples per mother") + theme(legend.position = "none")

c <- ggplot_gtable(ggplot_build(c))

#################################################

### generate boxplots ###

gt <- read.table("/sandbox/rs2305957_f.gt")
names(gt) <- c("case", "sampleid", "genotype")
pheno <- read.table("/sandbox/all_nonseg_mitotic_blastomere_f.matlab.response", sep=',')
names(pheno) <- c("case", "controls", "cases")
gt <- merge(gt, pheno, "case")

# remove missing genotypes
gt <- gt[gt$genotype!="00",]

# plot proportions stratified by genotype, requiring at least 3 blastomeres per mother
a <- ggplot(data = gt[gt$controls + gt$cases > 2,], aes(fill = factor(genotype), x = genotype, y = (cases / (controls + cases)))) + geom_boxplot() + xlab("Genotype") + ylab("Prop. blastomeres w/ mitotic error") + ylim(0,1) + ggtitle("Discovery") + theme(legend.position="none")
mean(gt[gt$genotype == "GG",]$cases / (gt[gt$genotype == "GG",]$controls + gt[gt$genotype == "GG",]$cases))
mean(gt[gt$genotype == "AG",]$cases / (gt[gt$genotype == "AG",]$controls + gt[gt$genotype == "AG",]$cases))
mean(gt[gt$genotype == "AA",]$cases / (gt[gt$genotype == "AA",]$controls + gt[gt$genotype == "AA",]$cases))

a <- ggplot_gtable(ggplot_build(a))

gt_validate <- read.table("/sandbox/rs2305957_validate_f.gt")
names(gt_validate) <- c("case", "sampleid", "genotype")
pheno_validate <- read.table("/sandbox/validate_nonseg_mitotic_blastomere.results", sep=',')
names(pheno_validate) <- c("case", "controls", "cases")
gt_validate <- merge(gt_validate, pheno_validate, "case")

mean(gt_validate[gt_validate$genotype == "GG",]$cases / (gt_validate[gt_validate$genotype == "GG",]$controls + gt_validate[gt_validate$genotype == "GG",]$cases))
mean(gt_validate[gt_validate$genotype == "AG",]$cases / (gt_validate[gt_validate$genotype == "AG",]$controls + gt_validate[gt_validate$genotype == "AG",]$cases))
mean(gt_validate[gt_validate$genotype == "AA",]$cases / (gt_validate[gt_validate$genotype == "AA",]$controls + gt_validate[gt_validate$genotype == "AA",]$cases))

# remove missing genotypes
gt_validate <- gt_validate[gt_validate$genotype!="00",]

# plot proportions stratified by genotype, requiring at least 3 blastomeres per mother
b <- ggplot(data = gt_validate[gt_validate$controls + gt_validate$cases > 2,], aes(fill = factor(genotype), x = genotype, y = (cases / (controls + cases)))) + geom_boxplot() + xlab("Genotype") + ylab("Prop. blastomeres w/ mitotic error") + ylim(0,1) + ggtitle("Validation") + theme(legend.position = "none")

b <- ggplot_gtable(ggplot_build(b))


#################################################

### plot effect size versus age ###

data <- read.table("~/Desktop/aneuploidy_calls.csv", sep = ",", header=T) # import the data

data_filtered <- filterData(data)

data_blastomere <- selectSampleType(data_filtered, blastomere)

data_blastomere <- callPloidy(data_blastomere)

aneuploid_binom <- aneuploidyByCase(data_blastomere)

gt <- read.table("/sandbox/rs2305957_f.gt")
names(gt) <- c("i", "sample_id", "genotype")
gt$numeric <- 0
gt[gt$genotype == "GG",]$numeric <- 0
gt[gt$genotype == "AG",]$numeric <- 1
gt[gt$genotype == "AA",]$numeric <- 2

aneuploid_binom <- merge(aneuploid_binom, gt, "i")

gt <- read.table("/sandbox/rs2305957_f.gt")
names(gt) <- c("case", "sampleid", "genotype")
pheno <- read.table("/sandbox/all_nonseg_mitotic_blastomere_f.matlab.response", sep=',')
names(pheno) <- c("case", "no_mitotic_error", "mitotic_error")
gt <- merge(gt, pheno, "case")

# remove missing genotypes
gt <- gt[gt$genotype != "00",]
maternal_age <- data.frame(cbind(data[!duplicated(data$case),]$case, data[!duplicated(data$case),]$maternal_age))
names(maternal_age)<-c("case", "maternal_age")
age_gt <- merge(gt, maternal_age, "case")
age_gt <- age_gt[complete.cases(age_gt),]
age_gt$prop <- age_gt$mitotic_error/ (age_gt$no_mitotic_error + age_gt$mitotic_error)

base<-2
mround <- function(x){ 
        2*round(x/2) 
} 

aneuploidyByGenotypeByAge <- function(data, genotype) {
	age_results_frame <- data.frame(matrix(ncol = 4))
	names(age_results_frame) <- c("prop", "se", "maternal_age", "genotype")
	range <- unique(mround(data[data$genotype == genotype,]$maternal_age))[order(unique(mround(data[data$genotype == genotype,]$maternal_age)))]
	for (i in range) {
		age_subset <- data[(mround(data$maternal_age) == i & data$genotype == genotype),]
		prop <- mean(age_subset$prop)
		se <- sqrt((prop * (1 - prop)) / nrow(age_subset))
		age_results_frame <- rbind(age_results_frame, c(prop, se, i, genotype))
	}
	age_results_frame<-age_results_frame[-1,]
	age_results_frame$prop <- as.numeric(as.character(age_results_frame$prop))
	age_results_frame$se <- as.numeric(as.character(age_results_frame$se))
	age_results_frame$maternal_age <- as.numeric(as.character(age_results_frame$maternal_age))
	return(age_results_frame)
}

gg <- aneuploidyByGenotypeByAge(age_gt, "GG")
ag <- aneuploidyByGenotypeByAge(age_gt, "AG")
aa <- aneuploidyByGenotypeByAge(age_gt, "AA")
age_results_frame <- rbind(gg, ag, aa)

limits <- aes(ymax = (prop + se), ymin = (prop - se))
d <- ggplot(data = age_results_frame, aes(x = maternal_age, y = prop, color = factor(genotype))) + geom_point() + geom_line() + geom_errorbar(limits, width = 0.5) + coord_cartesian(ylim = c(-0.05, 1.05)) + ylab("Prop. blastomeres w/ mitotic error") + theme(legend.justification = c(0,0), legend.position = c(.1,0.45)) + xlab('Maternal Age') + labs(fill = "") + scale_color_discrete(name = "Genotype")
d<-ggplot_gtable(ggplot_build(d))

#################################################

aneuploid_binom$prop <- aneuploid_binom$aneuploid_1 / (aneuploid_binom$aneuploid_1 + aneuploid_binom$aneuploid_0)

gg <- aneuploidyByGenotypeByAge(aneuploid_binom, "GG")
ag <- aneuploidyByGenotypeByAge(aneuploid_binom, "AG")
aa <- aneuploidyByGenotypeByAge(aneuploid_binom, "AA")
age_results_frame <- rbind(gg, ag, aa)

e <- ggplot(data = age_results_frame, aes(x = maternal_age, y = prop, color = factor(genotype))) + geom_point() + geom_line() + geom_errorbar(limits, width=0.5) + coord_cartesian(ylim = c(-0.05, 1.05)) + ylab("Prop. blastomeres w/ mitotic error") + theme(legend.justification = c(0,0), legend.position = c(.5,0)) + xlab('Maternal Age') + labs(fill = "") + scale_color_discrete(name = "Genotype")
e <- ggplot_gtable(ggplot_build(e))

#################################################

# put all panels in a single figure

b$widths <- a$widths
c$widths <- a$widths
d$widths <- a$widths
e$widths <- a$widths

grid.arrange(a, b, d, e, c, nrow = 3)
