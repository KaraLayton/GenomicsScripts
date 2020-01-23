##can use SNPrelate to prune for LD, calculate per locus Fst, PCA
library(SNPRelate)
library(gdsfmt)
library(SeqArray)
library(MASS)

##converting ped to gds
snpgdsPED2GDS("mappedordered_Poly.ped","mappedordered_Poly.map","all_snps.gds")
genofile<- snpgdsOpen("all_snps.gds")
genofile
get.attr.gdsn(index.gdsn(genofile,"snp.chromosome"))
(g <- read.gdsn(index.gdsn(genofile,"genotype"),start=c(1,1),count=c(5,3)))
get.attr.gdsn(index.gdsn(genofile,"genotype"))
head(read.gdsn(index.gdsn(genofile,"snp.id")))
head(read.gdsn(index.gdsn(genofile,"snp.rs.id")))
head(read.gdsn(index.gdsn(genofile,"snp.allele")))
allele <- read.gdsn(index.gdsn(genofile,"snp.allele"))
write.csv(allele, "allpolysnps_allele.csv")

##LD pruning
snpset<-snpgdsLDpruning(genofile, ld.threshold=0.2, autosome.only=FALSE, verbose=TRUE)
snpset.id<-unlist(snpset)

##PCA with pruned list
pca <- snpgdsPCA(genofile, snp.id = snpset.id, num.thread=2, autosome.only=FALSE)

##check variance
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

##plot PC1 loadings
SnpLoad <- snpgdsPCASNPLoading(pca, genofile)
dim(SnpLoad$snploading)
plot(SnpLoad$snploading[1,], type="h", ylab="PC 1")

##plot PCA
sample.id <-read.gdsn(index.gdsn(genofile, "sample.id"))
head(cbind(sample.id, pop))
tab <- data.frame(sample.id=pca$sample.id, pop = factor(pop)[match(pca$sample.id, sample.id)], EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors= FALSE)
plot(tab$EV1, tab$EV2, col=as.integer(tab$pop), pch=19, xlab="eigenvector 1", ylab="eigenvector 2") 
legend(locator(1),cex=0.7, legend=levels(tab$pop), pch=19, col=1:nlevels(tab$pop))

##plot top four PCs
lbls <- paste("PC",1:4,"\n",format(pc.percent[1:4],digits=2),"%",sep="")
pairs(pca$eigenvect[,1:4],col=tab$pop,labels=lbls)

#parallel coordinates plot
datpop <- factor(pop)[match(pca$sample.id, sample.id)]
parcoord(pca$eigenvect[,1:16],col=datpop)

#correlation between eigenvectors and SNP genotypes
chr <- read.gdsn(index.gdsn(genofile, "snp.position"))
CORR <- snpgdsPCACorr(pca, genofile, eig.which=1:4)
savepar <- par(mfrow=c(3,1),mai=c(0.3,0.55,0.1,0.25))
for(i in 1:3)
{
  plot(abs(CORR$snpcorr[i,]),ylim=c(0,1),xlab="",ylab=paste("PC",i),col=chr,pch="+")
}

##Dendrogram construction
##snpgdsDiss function calculates dissimilarities for each pair of individuals
diss <- snpgdsDiss(genofile,autosome.only=FALSE)
##hierarchical cluster analysis on the dissimilarity matrix
hc <- snpgdsHCluster(diss)
set.seed(100)
rv <- snpgdsCutTree(hc, label.H=FALSE, label.Z=FALSE)
snpgdsDrawTree(rv, main="All SNPs", edgePar=list(col=rgb(0.5, 0.5, 0.5, 0.75), t.col="black"))





