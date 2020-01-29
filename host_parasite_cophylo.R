##PARAFIT ANALYSIS AND COPHYLOGENY
#calculate phylo distances among host and parasite phylogenies
#use distance matrices as input for Parafit analysis
require(ape)
require(MASS)
my.phyloP <- read.tree("parasite_Apr27.tre")
my.phyloH <- read.tree("host_Apr27.tre")

P.dist <- cophenetic(my.phyloP)
H.dist <- cophenetic(my.phyloH)

write.matrix(P.dist, file="parasite phylo Apr 27")
write.matrix(H.dist, file="host phylo Apr 27")

#mark top half of matrix with NA & transform matrix to dataframe
#do this for P.dist too
H.dist[upper.tri(H.dist)] <- NA
cbind(which(!is.na(H.dist),arr.ind = TRUE),na.omit(as.vector(H.dist)))
H.dist.frame <- as.data.frame(as.table(H.dist))
write.csv(H.dist.frame, file="host phylo table")  

#comparing distances among host species that harbour the same parasite
#classified hosts into pair (closely related) and no-pair groups
means <-read.csv("host cophenetic pairs.csv")
means
boxplot(Distance ~ Group, data=means)
wilcox.test(Distance ~ Group, data = means)
t.test(Distance ~ Group, data=means)

#reading host csv file and transforming to matrix
host <- as.matrix(read.csv("host matrix Apr 27.csv", sep=",", header = TRUE, row.names=1))
host

#reading parasite csv file and transforming to matrix
parasite <- as.matrix(read.csv("parasite matrix Apr 27.csv", sep=",", header = TRUE, row.names=1))
parasite

#reading host parasite csv and transforming to matrix
#hosts as row names, parasites as column names, 1/0 for presence/absence
HP <- as.matrix(read.csv("HP Apr 27.csv",sep=",", header = TRUE, row.names=1))
HP

library(ape)
#ran with calliez correction because host matrix with negative eigenvalues
res <- parafit(host, parasite, HP, nperm=9999, correction="cailliez", test.links = TRUE, seed = NULL, silent = FALSE)
res
print(res)

#cophylo plot between host and parasite phylogenies
require(ape)
require(phytools)
host<- read.tree("host tree.tre")
parasite <- read.tree("parasite tree.tre")
assoc <- read.csv("assoc.csv") #2 column dataframe with host, parasite associations
obj <- cophylo(host, parasite, assoc=assoc)
plot(obj, link.type = "curved", link.lty = "solid", lwd= 1, pts= 0.75, fsize=0.5, gap = 10,  space = 30, length.line = 0.25)


