library(data.table)
devtools::install_github("petrelharp/local_pca/lostruct")
library(lostruct)
char_loc <- read_tped("plink.tped", triallelic=FALSE, phased=FALSE)
eigenstuff <- eigen_windows(char_loc, win=50, k=2)
windist <- pc_dist(eigenstuff, npc=2)
fit2d <- cmdscale(windist, eig=TRUE, k=2)
plot(fit2d$points,xlab="Coordinate 1",ylab="Coordinate 2", col=rainbow(1.2*nrow(windist)))

##subsetting map file by every 50th SNP to match window size in lostruct
MAPPOS <- read.delim2("new_char_Poly.map", header=F, stringsAsFactors = FALSE)
colnames(MAPPOS) <- c("CHR","SNP","CM","BP")
SNP_sub <- window(MAPPOS$SNP, deltat=50)
CHR_sub <- window(MAPPOS$CHR, deltat=50)
BP_sub <- window(MAPPOS$BP, deltat=50)
new_MAPPOS <- as.data.frame(cbind(CHR_sub,SNP_sub,BP_sub))
write.csv(new_MAPPOS,"new_MAPPOS.csv")
##removed last SNP from this MAPPOS file as didn't overlap with fit1d
MAPPOS2 <- read.csv("new_MAPPOS.csv", header=T, stringsAsFactors = FALSE)

#plot PC1 against location in the genome
library(ggplot2)
library(ggman)
library(dplyr)
eigenstuff <- eigen_windows(char_loc, win=50, k=100)
windist <- pc_dist(eigenstuff, npc=100)
fit1d <- cmdscale(windist, eig=TRUE, k=100)
P <- fit1d$points[,1]
PCMAP <-as.data.frame(cbind(MAPPOS2,P))
head(PCMAP)
ggman(PCMAP, chrom="CHR_sub", pvalue="P", snp="SNP_sub", bp="BP_sub", xlabel="Chromosome", ymax=10, pointSize=2) + theme_classic() + theme(axis.text.x=element_text(size=6,angle=45, hjust=1,colour="black"))
