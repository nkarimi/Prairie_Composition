library(ggplot2)
library(ggrepel)
library(dplyr)
library(reshape2)
library(tidyr)
#Part 1: mono-mixture expectation plot
#Part2: slope plot
#Part3: other slopeplot

###MONOCULTURES::
abund.mono <- read.csv("../DATA/2019_VegetationCover/Mono.2019.09.Plug.csv", as.is = TRUE)
abund.mono.mean <- abund.mono %>%
  group_by(sp) %>%
  summarise(mean_value_mono = mean(cover.new), SD_mono=sd(cover.new)) #, relCov = covTrans/SD) 
names(abund.mono.mean) <- c("Species", "mean_mono", "mono.SD")
  
##MIXTURES::
sppCovMat2019 <- read.csv("../DATA/2019_VegetationCover/COVER.sppCovMat.2019.09.Plug.csv")
#Test maximum cover data for both mono and mixtures (below)
abund.mix <- sppCovMat2019
#abund.spp <- intersect(abund.mono$sp, names(abund.mix[5:ncol(abund.mix)])) #127 checks out
abund.mix.long <- melt(abund.mix, id.vars = c("plot", "year", "type", "Code"), variable.name = "Species")
abund.mix.long2 <- abund.mix.long %>% drop_na()
abund.mix.long2 <- split(abund.mix.long2, abund.mix.long2$type)

abund.mix.long3 <- subset(abund.mix.long2$Plug, select=c(year, Species, value)) #sometimes Species is variable..?
abund.mix.mean <- abund.mix.long3 %>%
  group_by(Species) %>%
  summarise(mean_mix = mean(value), SD_mix=sd(value), CoV= (SD_mix/mean_mix)*100) #, relCov = covTrans/SD) 
#Coefficient of Variation = (Standard Deviation / Mean) * 100

names(abund.mix.mean) <- c("Species", "mix.cover", "SD_mix", "CoV")

#Merge mono and mixture:
abund.dat <- merge(abund.mix.mean, abund.mono.mean, by="Species")
abund.dat <- abund.dat %>%
  group_by(Species) %>%
  mutate(Proportional_abundance_in_mixture_relative_to_expectation = mix.cover / (mean_mono/ 15))

names(abund.dat) <- c("Species","Mixture abundance", "SD_mix","CoV" , "Monoculture abundance", "SD_mono", "Proportional abundance in mixture relative to expectation")

#drop rows with NaN or Inf
abund.dat[!grepl("NaN", abund.dat$`Proportional abundance in mixture relative to expectation`),] -> abund.dat
abund.dat[!grepl("Inf", abund.dat$`Proportional abundance in mixture relative to expectation`),] -> abund.dat
abund.dat[!grepl("Helianthus_strumosus", abund.dat$Species),] -> abund.dat

abund.dat$spEdit <- gsub('_', ' ', abund.dat$Species, fixed = T)
abund.dat$spEdit[abund.dat$`Proportional abundance in mixture relative to expectation` < 5.0 & abund.dat$`Monoculture abundance` > 15] <- ''
abund.dat$spEdit[abund.dat$`Proportional abundance in mixture relative to expectation` > 1.5 & abund.dat$`Monoculture abundance` < 75.00] <- ''
abund.dat$`Relative competitiveness` <- ifelse(abund.dat$`Monoculture abundance` > 50, 'High', 'Low')

p.abund <- ggplot(abund.dat,
            aes(x = `Monoculture abundance`,
                y = `Proportional abundance in mixture relative to expectation`,
                label = spEdit , size = abund.dat$CoV)) +
geom_point() + labs(size = "CoV") +
  scale_size_continuous(range = c(0, 6))

p.abund <- p.abund + geom_smooth(method = "lm", show.legend = FALSE) + theme(legend.position="right")
p.abund <- p.abund + geom_text_repel(aes(fontface='italic', color = `Relative competitiveness`),
                         size = 2.5,
                         segment.size = 0.15, show.legend = FALSE) + theme_bw() 

lmfit <- lm(`Monoculture abundance` ~ `Proportional abundance in mixture relative to expectation`, data=abund.dat)
summary(lmfit)
       
#ggsave(file="../OUT/mono-mix.expectations.pdf")

#####################GET MAX values#####################
Mono19.06 <- read.csv("../DATA/2019_VegetationCover/Mono.2019.06.Plug.csv", as.is = TRUE)
Mono19 <- read.csv("../DATA/2019_VegetationCover/Mono.2019.09.Plug.csv", as.is = TRUE)
abund.mono.max <- rbind(Mono19.06, Mono19) #take sp and mid columns
abund.mono.max <- subset(abund.mono.max, select=c(plot, sp, cover.new, year)) 

mono.max.mean <- abund.mono.max %>%
  group_by(sp, year) %>%
  summarise(mean_mono = mean(cover.new), SD_mono=sd(cover.new)) #, relCov = covTrans/SD) 
#get the max value:
mono.max <- mono.max.mean %>%
  group_by(sp) %>%
  summarize(max.cover = max(mean_mono))
names(mono.max) <- c("Species", "max.cover")
##MIXTURES::
sppCovMat2019 <- read.csv("../DATA/2019_VegetationCover/COVER.sppCovMat.2019.09.Plug.csv")
sppCovMat2019.06 <- read.csv("../DATA/2019_VegetationCover/COVER.sppCovMat.2019.06.Plug.csv")
#Get max for mixutres:
sppCovMat19.Max <- rbind(sppCovMat2019.06, sppCovMat2019) 
abund.dat.max <- sppCovMat19.Max
#abund.spp <- intersect(abund.mono$sp, names(abund.mix[5:ncol(abund.mix)])) #127 checks out
abund.mix.long <- melt(abund.dat.max, id.vars = c("plot", "year", "type", "Code"), variable.name = "Species")

#drop NA columns - remove rows from df.cover.long for which value = NA (but keeping those with cover = 0)
abund.mix.long2 <- abund.mix.long %>% drop_na()
abund.mix.long2 <- split(abund.mix.long2, abund.mix.long2$type)

abund.mix.long3 <- subset(abund.mix.long2$Plug, select=c(year, Species, value))
abund.mix.mean <- abund.mix.long3 %>%
  group_by(Species, year) %>%
  summarise(mean_mix = mean(value), SD_mix=sd(value)) #, relCov = covTrans/SD) 

#take the highest row (max mean):
abund.mix.max <- abund.mix.mean %>%
  group_by(Species) %>%
  summarize(max.cover = max(mean_mix), max.SD = max(SD_mix))

names(abund.mix.max) <- c("Species", "max.cover", "SD_max")

abund.dat.max <- merge(abund.mix.max, mono.max, by="Species")
names(abund.dat.max) <- c("Species", "max.cover", "SD_max", "mono.max")
abund.dat.max <- abund.dat.max %>%
  group_by(Species) %>%
  mutate(Proportional_abundance_in_mixture_relative_to_expectation = max.cover / (mono.max/ 15))
names(abund.dat.max) <- c("Species","Mixture abundance", "SD_mix", "Monoculture abundance", "Proportional abundance in mixture relative to expectation")

abund.dat.max$spEdit <- gsub('_', ' ', abund.dat.max$Species, fixed = T)
abund.dat.max$spEdit[abund.dat.max$`Proportional abundance in mixture relative to expectation` < 5.0 & abund.dat.max$`Monoculture abundance` > 20] <- ''
abund.dat.max$spEdit[abund.dat.max$`Proportional abundance in mixture relative to expectation` > 2 & abund.dat.max$`Monoculture abundance` < 75.00] <- ''
abund.dat.max$`Relative competitiveness` <- ifelse(abund.dat.max$`Monoculture abundance` > 50, 'High', 'Low')

#replace Inf with 0, NaN and outlier
abund.dat.max[!grepl("NaN", abund.dat.max$`Proportional abundance in mixture relative to expectation`),] -> abund.dat.max
abund.dat.max[!grepl("Inf", abund.dat.max$`Proportional abundance in mixture relative to expectation`),] -> abund.dat.max
abund.dat[!grepl("Asclepias_tuberosa", abund.dat$Species),] -> abund.dat.max

p.abundmax <- ggplot(abund.dat.max,
                  aes(x = `Monoculture abundance`,
                      y = `Proportional abundance in mixture relative to expectation`,
                      label = Species, size = abund.dat.max$SD_mix)) +
  geom_point() + labs(size = "SD in mixture cover") +
  scale_size_continuous(range = c(0, 10))

p.abundmax <- p.abundmax + geom_smooth(method = 'lm', show.legend = FALSE) + theme(legend.position="right")
p.abundmax <- p.abundmax + geom_text_repel(aes(fontface='italic', color = `Relative competitiveness`),
                                     size = 2.5,
                                     segment.size = 0.15, show.legend = FALSE) + theme_bw()

#ggsave(file="../OUT/mono-mix.expectations.png")

#########################################
#SLOPEPLOT
####RESCALE FOR SLOPEPLOTS ONLY
        is.na(abund.dat)<-sapply(abund.dat, is.infinite)
        abund.dat[is.na(abund.dat)] <- 0
        
        #MEAN SCALING
        abund.dat$scale.mix <- scale(abund.dat$mean_value_mix)
        abund.dat$scale.mono <- scale(abund.dat$mean_value_mono)
        
        #Subset the data for plotting
        abund.dat.sub <- abund.dat
        abund.dat.sub <- subset(abund.dat, select=c(Species, scale.mix, scale.mono))
        names(abund.dat.sub) <- c("Species", "mix.cover", "mono.cover")
        
        #Make slope chart for comparing mixtures vs mono http://r-statistics.co/Top50-Ggplot2-Visualizations-MasterList-R-Code.html#Correlogram
        left_label <- paste(abund.dat.sub$Species) 
        right_label <- paste(abund.dat.sub$Species)
        abund.dat.sub$class <- ifelse( (abund.dat.sub$mono.cover - abund.dat.sub$mix.cover) < 0, "red", "green")
        
        
        p <- ggplot(abund.dat.sub) + geom_segment(aes(x=1, xend=2, y=mix.cover, yend=mono.cover, col=class), size=.75, show.legend=F) + 
          geom_vline(xintercept=1, linetype="dashed", size=.1) + 
          geom_vline(xintercept=2, linetype="dashed", size=.1) +
          scale_color_manual(labels = c("Up", "Down"), 
                             values = c("green"="#00ba38", "red"="#f8766d")) +  # color of lines
          labs(x="", y="Mean Cover") +  # Axis labels
          xlim(.5, 2.5) + ylim(0,(1.1*(max(abund.dat.sub$mix.cover, abund.dat.sub$mono.cover))))  # X and Y axis limits
        
        # Add texts
        p <- p + geom_text(label=left_label, y=abund.dat.sub$mix.cover, x=rep(1, NROW(abund.dat.sub)), hjust=1.1, size=3.3)
        p <- p + geom_text(label=right_label, y=abund.dat.sub$mono.cover, x=rep(2, NROW(abund.dat.sub)), hjust=-0.1, size=3.3)
        p <- p + geom_text(label="Mixtures", x=1, y=1.1*(max(abund.dat.sub$mix.cover, abund.dat.sub$mono.cover)), hjust=1.2, size=4)  # title
        p <- p + geom_text(label="Monocultures", x=2, y=1.1*(max(abund.dat.sub$mix.cover, abund.dat.sub$mono.cover)), hjust=-0.1, size=4)  # title
        
        # Minify theme
        p <- p + theme(panel.background = element_blank(), 
                  panel.grid = element_blank(),
                  axis.ticks = element_blank(),
                  axis.text.x = element_blank(),
                  panel.border = element_blank(),
                  plot.margin = unit(c(1,2,1,2), "cm"))
        
        p
#        ggsave(filename = "../OUT/scaled.mix-mono.SlopePlot.png", width = 8, height=16)
