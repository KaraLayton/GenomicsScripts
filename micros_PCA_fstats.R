library(adegenet)
library(ape)
library(pegas)
library(hierfstat)

charpca <- import2genind("clean_char_Feb11_micro.gtx")
charpca

charpcaH <- genind2hierfstat(charpca)
PWS <- pairwise.WCfst(charpcaH, diploid=TRUE)
write.csv(PWS, "PWFst_micros_16Jul.csv")

summary(charpca)
write.csv(charpca, file="charpca.csv")

test_stats <- basic.stats(charpca, diploid=TRUE)
test_stats

##identify population names
pops<-seppop(charpca)
names(pops)
head(pop(charpca),50)
table(pop(charpca))
barplot(table(pop(charpca)), col=funky(17), las=3, cex.names=0.75, xlab="Population", ylab="Sample Size")

temp <- locNames(charpca, withAlleles=TRUE)
head(temp, 10)
temp <- summary(charpca)
barplot((temp$loc.n.all), col=funky(17), las=3, cex.names=0.5, ylim=c(0, 25), xlab="Locus", ylab="Numbers of Alleles")
barplot((temp$pop.n.all), col=funky(17), las=3, cex.names=0.75, ylim=c(0,600), xlab="Population", ylab="Numbers of Alleles")

plot(temp$Hexp, temp$Hobs, pch=20, cex=3, xlim=c(0,0.9), ylim=c(0,1), xlab="Expected heterozygosity", ylab="Observed heterozygosity")
abline(0,1,lty=2)

##PCA
x.char <- scaleGen(charpca, NA.method=c("mean"), scale=FALSE)
x.char
pca.char <- dudi.pca(x.char, center=FALSE, scale=FALSE)
pca.char
s.label(pca.char$li)
s.class(pca.char$li, fac=pop(charpca), col=funky(15))

##examine first axis
s.class(pca.char$li, fac=pop(charpca),
        col=transp(funky(15),.6),
        axesel=FALSE, cstar=0, cpoint=3)
add.scatter.eig(pca.char$eig[1:50],3,1,2, ratio=.3)

##examine second axis
s.class(pca.char$li, fac=pop(charpca),
        xax=2, yax=3, col=transp(funky(15),.6),
        axesel=FALSE, cstar=0, cpoint=3)
add.scatter.eig(pca.char$eig[1:50],3,2,3, ratio=.3)


##checking locus names
locNames(charpop)

##checking allele names for each locus- only first 10 entries
alleles <- locNames(charpop, withAlleles=TRUE)
head(alleles, 10)

##fstats
library(hierfstat)
fstat(charpca)
fstat(charpca, fstonly=TRUE)

charpca.gtest <- gstat.randtest(charpca)
charpca.gtest

charpca.matFst <- pairwise.fst(charpca, res.type="matrix")
write.csv(charpca.matFst, file="charmatrix.csv")
charpca.matFst[1:10,1:4]


