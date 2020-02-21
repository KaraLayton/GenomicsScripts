##code for mixed-stock fishery analysis in Layton et al. (2020) Evolutionary Applications
library(tidyverse)
library(rubias)
library(stringr)
library(genepopedit)
library(rubias)
library(reshape2)
library(tibble)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(marmap)

######################################BASELINE ANALYSIS##############################################   
##reading in microsat baseline
charaw <- read.csv("Rubias_microsat_baseline.csv", colClasses = "character")
char1_base <- as_tibble(charaw)
saveRDS(char1_base, file="char1.rds", compress="xz")
char2_base <- readRDS("char1.rds")

##reading in SNP baseline
fst_500 <- genepop_flatten("NEW_training_500_panel.txt")
head(RF_500)
for(i in 1:ncol(RF_500)){
  RF_500[, i] = as.character(RF_500[, i])
}
AlleleEx <- max(sapply(RF_500[, 4], FUN = function(x) {
  nchar(as.character(x[!is.na(x)]))
}))

firstAllele <- as.data.frame(sapply(RF_500[, -c(1, 2, 3)], function(x) as.character(substring(x, 1, alleleEx/2))))
secondAllele <- as.data.frame(sapply(RF_500[, -c(1, 2, 3)], function(x) as.character(substring(x, (alleleEx/2) + 1, alleleEx))))

colnames(secondAllele) <- paste0(colnames(secondAllele), "_1")

d <- cbind(firstAllele, secondAllele)
indx <- rbind(names(firstAllele), names(secondAllele))
d <- d[, indx]

repunit <- rep(NA, nrow(d))
sample_type <- rep("reference", nrow(d))

d <- cbind(sample_type, repunit, RF_500[, 1:3],  d)
d$repunit <- d$Population
#ds <- d
##for(i in 1:length(ds)){
  ##ds[ds$Population == ds[i,"popnames"],"repunit"] = ds[i,"group"]#this function didn't work for me
##}
##d<-ds
head(d)
d <- d[ , c(1,2,4,3, 5:ncol(d))]
head(d)
colnames(d)[3] <- "collection"
colnames(d)[4] <- "indiv"
d2 <- as.data.frame(sapply(d[, -c(1:4)], function(x) as.character(as.factor(x))), stringsAsFactors = FALSE)
head(d2)
d3 <- cbind(d[, 1:4], d2)
head(d3)
d3[d3 == "000"] <- NA
rownames(d3[, 1:8]) #use d3 as input for self assignment

######################SELF ASSIGNMENT
micro_SA <- self_assign(char2_base, gen_start_col = 5) ##microsats
SNP_SA <- self_assign(d3), gen_start_col = 6) ##SNPs
#next steps all use micro_SA and char2_base files, but repeat with SNP files too

##grouping by repunit
sa_to_repu <- micro_SA %>% group_by(indiv, collection, repunit, inferred_repunit) %>% summarise(repu_scaled_like = sum(scaled_likelihood))
sa_to_repu$id <- paste(sa_to_repu$indiv, sa_to_repu$repu_scaled_like, sep = "_")
sa_to_repu2 <- sa_to_repu %>% group_by(indiv) %>% summarise(max_repu_scaled_like = max(repu_scaled_like)) %>% ungroup() %>% data.frame()
head(sa_to_repu2)
sa_to_repu2$id <- paste(sa_to_repu2$indiv, sa_to_repu2$max_repu_scaled_like,sep = "_")

test <- sa_to_repu[sa_to_repu$id %in% sa_to_repu2$id, ]
write.csv(x = test,"Baseline_Assign.csv",quote=F,row.names = F)

myresults<-table(test$repunit, test$inferred_repunit)
write.csv(x = myresults, file="self assign matrix.csv")

test$correct <- "No"
test$correct[test$repunit == test$inferred_repunit] = "Yes"
test_out <- data.frame(Pop = rownames(table(test$repunit, test$correct)), No = table(test$repunit, test$correct)[,1], Yes = table(test$repunit, test$correct)[,2])
test_out$Prop_correct <- test_out$Yes / (test_out$No + test_out$Yes)
test_out$Prop_correct<-round(as.numeric(test_out$Prop_correct),digits=3)
test_out
write.csv(test_out,file = "ReportingGroup_efficiency.csv",quote = F) #this gives you efficiency per repunit

new_acc_eff <- read.csv("sa_accur_effic_micros.csv", header=TRUE, stringsAsFactors = FALSE) #accuracy and efficiency for each repunit

comb_accu <- read.csv("sa_accur_effic_combined.csv", header=TRUE, stringsAsFactors = FALSE) ##accuracy and efficiency for each repunit for both micros and SNPs
       
comb_accu$repunit<-factor(comb_accu$repunit,levels=c("KAN","KOM","KOG","NAC","MCC","PAL","STC","NOR","SWA","KIY","R109","PAN","R105","IKATHRFOU","IKI",
                                                             "PUT","KIN","KAM","FRS","ANA","IKLREI","R78",
                                                             "ENG","MBB","PBP")) ##change repunits here depending on micro or SNP datasets

fig4a <- ggplot(data=comb_accu,aes(x=repunit, y=accuracy,fill=factor(dataset)))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values=c("grey51","coral"))+
  geom_point(data = comb_accu, aes(x = repunit, y=efficiency),position=position_dodge(width=1))+
  geom_line(data = comb_accu, aes(x = repunit, y=efficiency, group=dataset,colour=factor(dataset)),position=position_dodge(width=1),size=0.80) +
  scale_colour_manual(values=c("grey31","coral3"))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0,1.05))+
  labs(x="\nReporting Group",y="Accuracy-Efficiency\n")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,color="black"),
        axis.text.x = element_text(size=8, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color="black"),legend.position = "none",
        strip.background = element_rect(colour="black", fill=NA))
fig4a       
        
######################100% SIMULATIONS
six_hundy_scenarios <- lapply(char2_base$repunit, function(x) tibble(repunit = x, ppn = 1.0))
names(six_hundy_scenarios) <- paste("All", char2_base$repunit, sep = "_")
char_sim3 <- assess_reference_loo(reference = char2_base, 
                                           gen_start_col = 5, 
                                           reps = 50, 
                                           mixsize = 100,
                                           alpha_repunit = six_hundy_scenarios,
                                           alpha_collection = 10)

unique(char_sim3$repunit_scenario)
char_sim3<-char_sim3[which(char_sim3$true_pi==1 & char_sim3$post_mean_pi>0.50),]
char_sim3$repunit<-factor(char_sim3$repunit,levels=c("KAN","KOM","KOG","NAC","MCC","PAL","STC","NOR","SWA","KIY","R109","PAN","R105","IKATHRFOU","IKI",
                                                             "PUT","KIN","KAM","FRS","ANA","IKLREI","R78",
                                                             "ENG","MBB","PBP")) ##change repunits here depending on micro or SNP datasets

write.csv(char_sim3,file="charsim3_100Percent.csv",quote=F)

char_sim3 <- read.csv("charsim3_100Percent_SNP_micros.csv", header=TRUE, stringsAsFactors = FALSE) ##combined sim3 results for micros and SNPs

char_sim3$repunit<-factor(char_sim3$repunit,levels=c("KAN","KOM","KOG","NAC","MCC","PAL","STC","NOR","SWA","KIY","R109","PAN","R105","IKA","R104","R103","IKATHRFOU","IKI",
                                                     "PUT","KIN","KAM","FRS","ANA","IKL","REI","IKLREI","R78",
                                                     "ENG","MBB","PBP"))

fig4b<-ggplot(char_sim3,aes(x=repunit,y=mle_pi,colour=Dataset))+
  geom_boxplot(lwd=0.5)+
  scale_colour_manual(values=c("grey51","coral"))+
  #scale_y_continuous(labels=scales::percent,limits=c(0.5,1.0001),minor_breaks = 0.05)+
  geom_hline(yintercept = 0.944,linetype="solid",col="grey51", alpha=0.6, size=1)+
  geom_hline(yintercept = 0.939,linetype="solid",col="coral", alpha=0.6, size=1)+
  labs(x="\nReporting Group",y="Accuracy\n")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,color="black"),
        axis.text.x = element_text(size=8, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color="black"),legend.position = "none",
        strip.background = element_rect(colour="black", fill=NA))
fig4b
            
######################LOO SIMULATION           
char_sim1 <- assess_reference_loo(d3, 6, reps = 100, mixsize = 500)
char_sim1

write.csv(char_sim1,file="charsim1_defaultLOO.csv",quote=F)
d3$collection=d3$repunit
char_sim1$repunit<-factor(char_sim1$repunit,levels=c("ANA","FRS","KIN","KAM","IKL","REI","IKI","PUT","ENG","IKA","P103","P104","P105","KAN","KIY","KOG","KOM","MBB",
                                                     "MCC","NAC","NOR","PAL","PAN","PBP","P109",
                                                     "R78","STC","SWA"))

char_sim1 <- read.csv("charsim1_LOO_SNPs_micros.csv", header=TRUE, stringsAsFactors = FALSE) ##combined sim1 results for micros and SNPs

char_sim1$repunit<-factor(char_sim1$repunit,levels=c("KAN","KOM","KOG","NAC","MCC","PAL","STC","NOR","SWA","KIY","R109","PAN","R105","IKA","R104","R103","IKATHRFOU","IKI",
                                                             "PUT","KIN","KAM","FRS","ANA","IKL","REI","IKLREI","R78",
                                                             "ENG","MBB","PBP"))

figs3<-ggplot(char_sim1, aes(x = true_pi, y = post_mean_pi, colour = dataset)) +
  geom_point(size=2, alpha=0.6) +
  scale_colour_manual(values=c("grey","coral"))+
  labs(x="\nTrue Mixing Proportion",y="Simulated Mixing Proportion\n")+
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~repunit, ncol = 3)+
  scale_y_continuous(breaks = seq(0, 0.20, by = 0.1))+
  scale_x_continuous(breaks=seq(0, 0.20, by = 0.1))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color="black"),legend.position = "none",
        strip.background = element_rect(colour="black", fill=NA))
figs3

######################################MIXTURE ANALYSIS##############################################   
baselinedat4 <- read.csv("Rubias_microsat_baseline_mix.csv", colClasses = "character")
baselinedat4 <- as_tibble(baselinedat4)
saveRDS(baselinedat4, file="baselinedat4.rds", compress="xz")
baselinedat4[baselinedat4=="0"]<-NA

mymix <- read.csv("Rubias_microsat_fisheries_mix.csv", colClasses = "character")
mymix <- as_tibble(mymix)
saveRDS(mymix, file="mymix.rds", compress="xz")
mymix[mymix=="0"]<-NA

test_MIX <- infer_mixture(reference = baselinedat4, mixture = mymix, gen_start_col = 5,reps = 20000,burn_in = 1000)

##mixing proportions
rep_mix_ests <- test_MIX$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi))  # adding mixing proportions over collections in the repunit

write.csv(rep_mix_ests,"rep_mix_ests_MCMC.csv") ##tells you the proportion of individuals in the mixture that assigned to each repunit
       
WGM_mix_EST_IndPosts <- test_MIX$indiv_posteriors
WGM_mix_EST_IndPosts_res <- WGM_mix_EST_IndPosts %>% group_by(indiv) %>% top_n(1, PofZ) %>% ungroup() ##top match
hist(WGM_mix_EST_IndPosts_res$PofZ)
mean(WGM_mix_EST_IndPosts_res$PofZ)
write.csv(WGM_mix_EST_IndPosts_res[,-10], "Mixture_IndividualProbabilities_MCMC.csv",quote=F,row.names = F,append=T)

##number of mixture individuals assigning (with >80% probability) to each repunit
indivs_80 <- read.csv("Mixture_IndividualProbabilities_MCMC_80.csv", header=TRUE, stringsAsFactors = FALSE) ##only individuals that assigned with >80% prob

indivs_80$mixture_collection <-  factor(indivs_80$mixture_collection,
                             levels=c("SKF","NAF","HOP","MKK","PTV","LMV","BTK","CTW","SLW"))

indivs_80$repunit <-  factor(indivs_80$repunit,
                                        levels=c("STC","SWA","KIY","PAN","R109","IKATHRFOU","R105","IKI","PUT",
                                                 "KIN","KAM","FRS","IKLREI","R78","ENG","MBB"))

fig5a<-ggplot(indivs_80, aes(x = mixture_collection, y = PofZ, fill = repunit)) + 
  #facet_grid(~Continent,scales="free_x",space="free_x")+
  geom_bar(stat = "identity") + 
  scale_fill_manual(values=c("olivedrab2","indianred1","red","green3","deeppink3","mediumseagreen","darkorchid1",
                             "yellow","darkorange2","blue4","peru","blue","darkseagreen1","cyan","darkgoldenrod","forestgreen"),name="Reporting Group")+
  #scale_fill_brewer(palette="Rainbow", name="Reporting Group")+
  xlab("Fishery")+
  ylab("Individuals assigned")+ 
  #scale_y_continuous(expand=c(0,0)) + 
  #coord_cartesian(xlim=c(1,16.5),ylim = c(0, 350))+
  #geom_segment(data = anno_seg, aes(x = x1, xend= x2, y = y1, yend= y2), inherit.aes=FALSE, linetype=2, size=0.6) +
  #geom_text(data=anno_seg, aes(x=x1, y=y3), label=c("SKF","NAF","HOP","MKK, PTV","LMV","BTK, CTW, SLW"),inherit.aes=FALSE,angle=90, nudge_x=0.2, size=4.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=14,color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black"),
        axis.text.y = element_text(color="black"), legend.position=c(0.85,0.95), legend.justification=c(0.01,0.95), legend.direction="vertical",
        legend.text=element_text(size=8), legend.background = element_rect(color=NA), legend.key.size=unit(0.4,"cm"), legend.title=element_text(size=8), strip.background = element_blank(),
        strip.text.x = element_blank())
fig5a

##looking at proportion of individuals assigning to a repunit and the repunit's drainage area
drain_mix <- read.csv("drainage_proportion.csv", header=TRUE, stringsAsFactors = FALSE) ##four column dataframe; mixture, repunit, repunit drainage area, mean proportion 
NAF_drain <- subset(drain_mix, mixture=="NAF")
SKF_drain <- subset(drain_mix, mixture=="SKF")
NAF_SKF_drain <- rbind(NAF_drain, SKF_drain)

fit <- lm(NAF_SKF_drain$mean.proportion ~ NAF_SKF_drain$repunit.drainage)
summary(fit)

fig5b<-ggplot(NAF_SKF_drain, aes(x = repunit.drainage, y = mean.proportion)) + 
  geom_point(size=2.5) +
  geom_smooth(method = "lm", se = FALSE, size=0.5) +
  #scale_color_manual(values=c("darkmagenta","yellow","blue4","darkgoldenrod","forestgreen","blue","red","green",
  #"cornflowerblue","indianred1","darkslateblue", "deeppink3","cyan","darkorange2"), name="Reporting Group")+
  xlab("Drainage area (km)")+
  ylab("Mixture proportion")+ 
  scale_y_continuous(breaks = seq(0, 0.62, by = 0.1)) + 
  coord_cartesian(ylim = c(0, 0.62))+
  #scale_x_continuous(breaks = seq(0, 200, by = 50)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=14)) + facet_wrap(~mixture, scales="free")
fig5b

##calculaye least cost distance in marmap                                           
geodata <- read.csv("geosphere_dist.csv", header=TRUE, row.names=1)
bathydata <- getNOAA.bathy(-55,-65,50,60, res=1,keep=T)
summary.bathy(bathydata)
plot(bathydata, image=TRUE)
blues <- colorRampPalette(c("red","purple","blue",
                            "cadetblue1","white"))
plot(bathydata, image = TRUE, bpal = blues(100))
trans1 <- trans.mat(bathydata)
cost <- lc.dist(trans1,geodata,res="dist")

##number of individuals assigned to each repunit and distance between repunit and fishery
ind_dist <- read.csv("rep_inds_dist_80.csv", header=TRUE, stringsAsFactors = FALSE) ##includes geographic distance, proportion assigned to repunit, and number of individuals assigning with >80% probability

ind_dist$mixture_collection <-  factor(ind_dist$mixture,
                                        levels=c("BTK","CTW","HOP","LMV","MKK","NAF","PTV","SKF","SLW"))

fig5c <- ggplot(data=ind_dist, aes(x = geodist, y = individuals, color = mixture)) + geom_point(size=3.5, pch=19) + theme_classic() + xlab("\nGeographic Distance (km)")+
  ylab("Individuals Assigned\n (>80%)\n") + theme(text=element_text(size=14),legend.position=c(1,0.7), legend.justification=c(0.95,0.5))+ scale_color_brewer(palette="Paired", name="Fishery")
fig5c

##proportion of mixture assigning to each repunit with dotted line indicating location of fishery
#calculate credible intervals and remove repunits that have CI proportions encompassing zero
top25 <- rep_mix_ests %>%
  filter(mixture_collection == "SLW") %>% ##swap in different mixture codes here
  arrange(desc(repprop)) %>%
  slice(1:25)

trace_subset <- test_MIX$mix_prop_traces %>%
  filter(mixture_collection == "SLW", sweep > 200) %>%
  group_by(sweep, repunit) %>%
  summarise(repprop = sum(pi)) %>% 
  filter(repunit %in% top25$repunit)

sum_SLW <- trace_subset %>%
  group_by(repunit) %>%
  summarise(loCI = round(quantile(repprop, 0.05),4),
            hiCI = round(quantile(repprop, 0.95),4),
            mean = round(mean(repprop),4))

write.csv(sum_SLW,"SLW_CIs.csv",quote=F)

rep_mix_CIs <- read.csv("Rubias_MixtureProps_withCIs_0removed.csv", header=TRUE, stringsAsFactors = FALSE) ##four column dataframe with repunit, lower CI, upper CI, mean proportion

rep_mix_CIs$repunit<-factor(rep_mix_CIs$repunit,levels=c("KAN","KOM","KOG","NAC","MCC","PAL","STC","NOR","SWA","KIY","R109","PAN","R105","IKATHRFOU","IKI",
                                                             "PUT","KIN","KAM","FRS","ANA","IKLREI","R78",
                                                             "ENG","MBB","PBP"))
rep_mix_CIs$mixture <-  factor(rep_mix_CIs$mixture,
                                 levels=c("SKF","NAF","HOP","MKK","PTV","LMV","BTK","CTW","SLW"))

figs4<-ggplot(rep_mix_CIs, aes(x = repunit, y = mean*2, fill = repunit)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values=c("darkmagenta","yellow","blue4","darkseagreen1","forestgreen","blue","red","green",
                             "peru","mediumseagreen","indianred1","darkorchid1",
                             "darkslateblue","deeppink3","cyan","darkorange2","dimgray","olivedrab2","burlywood3","lightpink1","cornflowerblue","black","grey"
                             ,"salmon","darkgoldenrod"))+ 
  labs(x="Reporting Unit",y=expression("Mixture Proportion"))+
  geom_errorbar(aes(x=repunit,ymin=loCI*2, ymax=hiCI*2), width=.5,color="black")+
  geom_segment(data = anno_seg, aes(x = x1, xend= x2, y = y1, yend= y2), inherit.aes=FALSE, linetype=2, size=0.6) +
  scale_y_continuous(expand=c(0,0),breaks = seq(0, 1, by = 0.2)) + 
  coord_cartesian(ylim = c(0, 1))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=16,color="black"),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust=.4,color = "black"),
        axis.text.y = element_text(color="black"),legend.position = "none")

figs4

figs4 + facet_wrap(~mixture) + expand_limits(x=17)

anno_seg <- data.frame(x1 = c(16.5, 16.5, 14.5, 15.5, 14.5, 10.5, 15.5, 2.5, 16.5), x2 = c(16.5, 16.5, 14.5, 15.5, 14.5, 10.5, 15.5, 2.5, 16.5), 
                       y1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0), y2 = c(1, 1, 1, 1, 1, 1, 1, 1, 1), 
                       mixture = c("BTK", "CTW", "HOP", "LMV", "MKK", "NAF", "PTV", "SKF", "SLW"))


                                         