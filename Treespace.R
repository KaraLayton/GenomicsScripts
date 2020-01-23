##thanks to Michelle Kendall for help including ASTRAL trees in this analysis
library(ape)
library(treespace)
teas_boot <- read.tree("Teasdale.ufboot")
all_boot <- read.tree("All_genes.ufboot")
mtDNA_boot <- read.tree("mtDNA_boots.tre")
ASTRAL_boot <- read.tree("Partial_ASTRAL_ten.tre")
teas_ML <- read.tree("Teasdale_ML.treefile")
all_ML <- read.tree("All_genes_ML.treefile")
mtDNA_ML <- read.tree("mtDNA_ML.tre")

#testing ASTRAL tree
ASTRALA <- read.tree("Partial_ASTRAL_tenA.tre")
ASTRALB <- read.tree("Partial_ASTRAL_tenB.tre")
ASTRAL <- c(ASTRALA, ASTRALB)
class(ASTRAL) <- "multiPhylo"
names(ASTRAL)[1:10] <- paste0("ASTRALA",1:10)
names(ASTRAL)[11:20] <- paste0("ASTRALB",11:20)
Dtype <- c(rep("ASTRALA",10),rep("ASTRALB",10))
Dscape <- treespace(ASTRAL, nf=5)
plot(ASTRAL)
ASTRAL_bootR$tip.label
ASTRAL_bootR$edge.length

#rooting trees
teas_bootR <- root(teas_boot, resolve.root=TRUE, outgroup="WAMS92138")
all_bootR <- root(all_boot, resolve.root=TRUE, outgroup="WAMS92138")
mtDNA_bootR <- root(mtDNA_boot, resolve.root=TRUE, outgroup="WAMS92138")
ASTRAL_bootR <- root(ASTRAL_boot, resolve.root=TRUE, outgroup="WAMS92138")
teas_MLR <- root(teas_ML, resolve.root=TRUE, outgroup="WAMS92138")
all_MLR <- root(all_ML, resolve.root=TRUE, outgroup="WAMS92138")
mtDNA_MLR <- root(mtDNA_ML, resolve.root=TRUE, outgroup="WAMS92138")

#collect trees into single multiPhylo object; BEAST trees and bootstrapped and non-bootstrapped versions of NJ and ML
ChromTrees <- c(teas_bootR, all_bootR, mtDNA_bootR, ASTRAL_bootR, teas_MLR, all_MLR, mtDNA_MLR)
class(ChromTrees) <- "multiPhylo"
write.nexus(ChromTrees, file="ChromTrees")

#add tree names
names(ChromTrees)[1:1000] <- paste0("teas_bootR",1:1000)
names(ChromTrees)[1001:2000] <- paste0("all_bootR",1001:2000)
names(ChromTrees)[2001:3000] <- paste0("mtDNA_bootR",2001:3000)
names(ChromTrees)[3001:3010] <- paste0("ASTRAL_bootR",3001:3010)
names(ChromTrees)[[3011]] <- "teas_MLR"
names(ChromTrees)[[3012]] <- "all_MLR"
names(ChromTrees)[[3013]] <- "mtDNA_MLR"

#create vector corresponding to tree inference method
Dtype <- c(rep("teas_bootR",1000),rep("all_bootR",1000),rep("mtDNA_bootR",1000),rep("ASTRAL_bootR",10),"teas_MLR","all_MLR","mtDNA_MLR")
#treescape to find and project differences
Dscape <- treespace(ChromTrees, nf=5)
plotGrovesD3(Dscape$pco, groups=Dtype)

#refining the plot with different colours, symbols, legend
Dcols <- c("darkblue","firebrick", "forestgreen", "cornflowerblue", "coral", "darkseagreen")
Dmethod <- c(rep("Teasdale_boots",1000),rep("AllGenes_boots",1000),rep("mtDNA_boots",1000),"Teasdale_ML", "AllGenes_ML", "mtDNA_ML")
Dbootstraps <- c(rep("replicates",3000), "Teasdale_ML", "AllGenes_ML", "mtDNA_ML")
Dhighlight <- c(rep(1,3000),2,2,2)
plotGrovesD3(Dscape$pco, 
             groups=Dmethod, 
             colors=Dcols,
             col_lab="Dataset",
             size_var=Dhighlight,
             size_range = c(100,300),
             size_lab="",
             symbol_var=Dbootstraps,
             symbol_lab="",
             point_opacity=c(rep(0.4,300),1,1,1), 
             legend_width=200)

#adding tree labels to plot
#here you can identify outlier trees/topologies
plotGrovesD3(Dscape$pco, 
             groups=Dmethod, 
             treeNames = names(ChromTrees), # add the tree names as labels
             colors=Dcols,
             col_lab="Dataset",
             size_var=Dhighlight,
             size_range = c(100,300),
             size_lab="",
             symbol_var=Dbootstraps,
             symbol_lab="",
             point_opacity=c(rep(0.4,300),1,1,1), 
             legend_width=200)

plotGrovesD3(Dscape$pco, 
             groups=Dmethod, 
             tooltip_text = names(ChromTrees), # add the tree names as tooltip text
             colors=Dcols,
             col_lab="Dataset",
             size_var=Dhighlight,
             size_range = c(100,300),
             size_lab="",
             symbol_var=Dbootstraps,
             symbol_lab="",
             point_opacity=c(rep(0.4,300),1,1,1), 
             legend_width=200)

#scree plot
barplot(Dscape$pco$eig, col="navy")

#view plot in 3D
library(rgl)
Dcols3D <- c(rep(Dcols[[1]],1000),rep(Dcols[[2]],1000),rep(Dcols[[3]],1000),Dcols[[3]],Dcols[[4]],Dcols[[5]])
rgl::plot3d(Dscape$pco$li[,1],Dscape$pco$li[,2],Dscape$pco$li[,3],Dscape$pco$li[,4],
            type="s",
            size=c(rep(1.5,300),3,3,3), 
            col=Dcols3D,
            xlab="", ylab="", zlab="")