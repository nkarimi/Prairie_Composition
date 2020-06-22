library(ape)
library()
library(betapart)
library(picante)
library(cowplot)
#Phylo-Ordinations
#per this tutorial: http://kembellab.ca/r-workshop/biodivR/SK_Biodiversity_R.html

#LOAD DATA: Species x plots cover matrix (NAs) for yr 1 and yr 2:
#comm matrix with rows as sites and columns are species

head(df.cover)
#Each cell contains the percent cover of a species in a sample.
#Many multivariate methods are sensitive to the total abundance in a sample, so we should probably convert these absolute abundance estimates to a relative abundance estimate. We can do this with a function from the vegan package.
# Turn percent cover to relative abundance by dividing each value by sample total abundance
df.cover[8:ncol(df.cover)] -> comm
paste(df.cover$plot, df.cover$year, sep="_") -> rownames(comm)
comm <- decostand(comm, method = "total")
#apply(comm, 1, sum) # check total abundance in each sample should be 1

phy <- read.tree("../DATA/phylogeny.analyzed.2016-01-05b.tre") 
phy$tip.label    #140
phy$tip.label[which(phy$tip.label == "Symphyotrichum_novae-angliae")] <-
  "Symphyotrichum_novae.angliae"
tree <- drop.tip(phy, which(!phy$tip.label %in% colnames(comm))) 
#convert phylogenety to a distance matrix
phy.dist <- cophenetic(tree)

##########Part 1 MNTD & MPD beta diversity

#calculate phylogenetic MNTD beta diversity
comm.mntd.dist <- comdistnt(comm, phy.dist, abundance.weighted = TRUE)
# comm.mntd.traits.dist <- comdistnt(comm, trait.dist, abundance.weighted = TRUE)
# calculate Mantel correlation for taxonomic Bray-Curtis vs. phylogenetic MNTD diversity
#mantel(comm.bc.dist, comm.mntd.dist)

comm.mpd.dist <- comdist(comm, phy.dist, abundance.weighted = TRUE)

#Part 2 Ordinations
# NMDS ordination of phylogenetic distances - use monoMDS since we only have among-sample distances
comm.mntd.mds <- monoMDS(comm.mntd.dist, k=3)
comm.mpd.mds <- monoMDS(comm.mpd.dist, k=3)

mds.fig <- ordiplot(comm.mntd.mds, type = "none")
points(mds.fig, "sites", pch = 19, col = "green", select = df.cover$phyD == 
         "L")
points(mds.fig, "sites", pch = 19, col = "blue", select = df.cover$phyD == 
         "M")
points(mds.fig, "sites", pch = 19, col = "red", select = df.cover$phyD == 
         "H")
ordiellipse(comm.mntd.mds, df.cover$phyD, conf = 0.95, label = TRUE)

###Convert for more complex plotting:
as.data.frame(scores(comm.mntd.mds)) -> Sc.mntd.mds
Sc.mntd.mds$year <- df.cover$year 
Sc.mntd.mds$phyD <- as.factor(df.cover$phyD)

#ggplot() + 
#  stat_ellipse(data=Sc.mntd.mds, aes(x=MDS1,y=MDS2, color= year)) +
#  geom_point(data =Sc.mntd.mds, aes(MDS1, MDS2, color = year),size=3) +
#  facet_grid(phyD ~ .)
  #geom_path(data=centroid.Phy.Plug, aes(x = NMDS1, y = NMDS2, group = phyD), arrow = arrow(length = unit(0.15, "cm")), colour="black", size = .85) + 
  #labs(x="nmds1", y="nmds2", color="phyD") + theme_bw() 
#above too messy. subset data:

start.ord <- subset(Sc.mntd.mds, year == "2017")
spring.ord <- subset(Sc.mntd.mds, year=="2019.06")
end.ord <- subset(Sc.mntd.mds, year=="2019.09")
  
mntd.ords <- ggplot() +  
  geom_point(data=start.ord, aes(x = MDS1, y = MDS2, color=year), size = 1) +
  stat_ellipse(type="t", data=start.ord, aes(x=MDS1,y=MDS2, color=year), size=1) +
  
  stat_ellipse(type="t", data=spring.ord, aes(x=MDS1,y=MDS2, color=year), size=1) +
  geom_point(data=spring.ord, aes(x = MDS1, y = MDS2, color = year), size = 1) +
  
  stat_ellipse(type="t", data=end.ord, aes(x=MDS1,y=MDS2, color=year), size=1 ) +
  geom_point(data=end.ord, aes(x = MDS1, y = MDS2, color = year), size = 1) +
  facet_grid(phyD ~ .) +
  theme_bw() +
  labs(title="MNTD") + scale_color_manual(values=c('#44AA99','#332288', '#999933')) +  #"#e66101", "#5e3c99","#33a02c"
  theme(strip.background = element_rect(color="grey", fill="white", size=1.5)) +
  theme(legend.position = "none")
#  "#88CCEE", "#CC6677", "#AA4499", "#882255"))
  
####PLOT MPD
as.data.frame(scores(comm.mpd.mds)) -> Sc.mpd.mds
Sc.mpd.mds$year <- df.cover$year 
Sc.mpd.mds$phyD <- as.factor(df.cover$phyD)

start.ord.mpd <- subset(Sc.mpd.mds, year == "2017")
spring.ord.mpd <- subset(Sc.mpd.mds, year=="2019.06")
end.ord.mpd <- subset(Sc.mpd.mds, year=="2019.09")

mpd.ords <- ggplot() +  
  stat_ellipse(type="norm", data=start.ord.mpd, aes(x=MDS1,y=MDS2, color=year), size=1) +
  geom_point(data=start.ord.mpd, aes(x = MDS1, y = MDS2, color=year), size = 1) +
  
  stat_ellipse(type="norm", data=spring.ord.mpd, aes(x=MDS1,y=MDS2, color=year), size=1) +
  geom_point(data=spring.ord.mpd, aes(x = MDS1, y = MDS2, color = year), size = 1) +
  
  stat_ellipse(type="norm", data=end.ord.mpd, aes(x=MDS1,y=MDS2, color=year), size=1) +
  geom_point(data=end.ord.mpd, aes(x = MDS1, y = MDS2, color = year), size = 1) +
  facet_grid(phyD ~ .) +
  theme_bw() +
  labs(title="MPD") + scale_color_manual(values=c('#44AA99','#332288', '#999933')) +  ## #882255
  theme(strip.background = element_rect(color="grey", fill="white", size=1.5))

cowplot::plot_grid(mntd.ords, mpd.ords, labels = "auto", rel_widths = c(1,1.4))
#ggsave(file="../OUT/phylo_ordinations.pdf", width=6, units ="in")
#ggsave(file="../OUT/phylo_ordinations_sm.pdf", scale=1, width=3, units ="in")

##Add the cover taxonmic order here:
#NMDS.plug.comp
start.NMDS.plug.comp <- subset(NMDS.plug.comp , year == "2017")
end.NMDS.plug.comp <- subset(NMDS.plug.comp , year == "2019.09")
springNMDS.plug.comp <- subset(NMDS.plug.comp , year == "2019.06")
  
taxonomic.ords <- ggplot() +  
  stat_ellipse(data=start.NMDS.plug.comp, aes(x=NMDS1,y=NMDS2, color=year)) +
  geom_point(data=start.NMDS.plug.comp, aes(x = NMDS1, y = NMDS2, color=year), size = 1) +
  
  stat_ellipse(data=springNMDS.plug.comp, aes(x=NMDS1,y=NMDS2, color=year)) +
  geom_point(data=springNMDS.plug.comp, aes(x = NMDS1, y = NMDS2, color = year), size = 1) +
  
  stat_ellipse(data=end.NMDS.plug.comp, aes(x=NMDS1,y=NMDS2, color=year) ) +
  geom_point(data=end.NMDS.plug.comp, aes(x = NMDS1, y = NMDS2, color = year), size = 1) +
  facet_grid(phyD ~ .) +
  theme_bw() +
  labs(title="Taxonomic")



