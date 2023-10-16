library(readr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggbiplot)
library(ggpubr)
library(rstatix)
library(plotly)
library(NMF)
library(tidyverse)

######### Metadata #########
samplemap <- read_tsv("results/samplemap.tsv")
samplemap$group <-factor(samplemap$group, levels=c("Healthy", "SCCHN", "OMD", "PMD"))

samplemap$treatment <-factor(samplemap$treatment, levels=c("BL", "T1", "T2", "F1", "F2", "F3", "P1"))

########## Isolated DNA ############
p <- ggplot() +
  geom_bar(data=samplemap[samplemap$dropout == "no",], aes(x=sample, 
                                                           y=concentration,
                                                           fill = group),
           stat = "identity") +
  theme_classic()+ 
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

anov <- compare_means(concentration ~ group,
                      data = samplemap[samplemap$dropout == "no",],
                      method = "anova")

pwise <- compare_means(concentration ~ group,
                       data = samplemap[samplemap$dropout == "no",],
                       ref.group = "Healthy",
                       method = "t.test")
p <- ggplot(data=samplemap[samplemap$dropout == "no",],
            aes(x=group, y=concentration)) +
  geom_boxplot(aes(fill = group), outlier.shape = NA) + geom_point() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  stat_compare_means(method = "anova") +
  theme_classic() + 
  stat_pvalue_manual(pwise, label = "p.signif", y.position = 13.5, step.increase = 0.04) +
  ylab("Concentration ng/μl")
p

anov <- compare_means(concentration ~ group,
                      data = samplemap[samplemap$dropout == "no" & samplemap$day=="0",],
                      method = "anova")

pwise <- compare_means(concentration ~ group,
                       data = samplemap[samplemap$dropout == "no" & samplemap$day=="0",],
                       ref.group = "Healthy",
                       method = "t.test")

p <- ggplot(data=samplemap[samplemap$dropout == "no"& samplemap$day=="0",],
            aes(x=group, y=concentration)) +
  geom_boxplot(aes(fill = group), outlier.shape = NA) + geom_point() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  stat_compare_means(method = "anova") +
  theme_classic() + 
  ylab("Concentration ng/μl")
p

p <- ggplot(data=samplemap[samplemap$dropout == "no",],
            aes(x=sample, y=concentration)) +
  geom_bar(stat="identity", aes(fill = group)) +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + 
  ylab("Concentration ng/μl") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_grid(.~group, scales = "free", space = "free_x")
p


################################################################################
######### Tumor fraction data #########
tf <- read_tsv("results/tf_ichorcna.tsv")

#prep for the size-selection analysis
tflong <- melt(tf, 
               id.vars=c("sample"),
               variable.name="size_selection",
               value.name="tumor_fraction")
tflong$size_selection <- as.character(tflong$size_selection)
tflong$size_selection[tflong$size_selection == "tumor_fraction_275-325"] <-"275-325"
tflong$size_selection[tflong$size_selection == "tumor_fraction"] <- "all"
tflong$size_selection <-factor(tflong$size_selection, levels=c("all", "275-325"))

tf <- merge(x = tf, y = samplemap, by = "sample", all.x = TRUE)
tflong <- merge(x = tflong, y = samplemap, by = "sample", all.x = TRUE)

p <- ggplot(tf[tf$dropout=="no",], aes(x=group, y=tumor_fraction, fill=group)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p <- ggplot(tf[tf$dropout=="no" & tf$day=="0",], aes(x=group, y=tumor_fraction)) + 
  geom_boxplot(aes(fill=group), outlier.shape = NA) +
  geom_point() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.title.x=element_blank(),
                          axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent) + ylab("Tumor fraction")
p

p <- ggplot(tf[tf$dropout=="no" & tf$day=="0",], aes(x=group, y=`tumor_fraction_275-325`)) + 
  geom_boxplot(aes(fill=group), outlier.shape = NA) + geom_point() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent)
p

p <- ggplot(tf[tf$dropout=="no",], aes(x=group, y=`tumor_fraction_275-325`)) + 
  geom_boxplot(aes(fill=group), outlier.shape = NA) + geom_point() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p <- ggplot(tf[tf$dropout=="no",], aes(x=group, y=tumor_fraction)) + 
  geom_boxplot(aes(fill=group), outlier.shape = NA) + geom_point() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

ggplotly(p)

p<-ggplot(tf[tf$dropout=="no",], aes(x=sample, y=tumor_fraction, fill=group)) +
  geom_bar(stat="identity")+ 
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.title.x=element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1))
p

p<-ggplot(tf[tf$dropout=="no" & tf$day=="0",], aes(x=sample, y=tumor_fraction, fill=group)) +
  geom_bar(stat="identity")+ 
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent)
p

anov <- compare_means(tumor_fraction ~ group,
                      data = tf[tf$dropout=="no" & tf$day=="0",],
                      method = "anova")

pwise <- compare_means(tumor_fraction ~ group,
                       data = tf[tf$dropout=="no" & tf$day=="0",],
                       ref.group = "Healthy",
                       method = "t.test")

greens <- c("#9CCC65", "#CDDC39", "#388E3C", "#4CAF50", "#00BFA5", "#689F38", "#B8FF12")
yellows <- c("#FFFF00",  "#FFD700", "#F9812A", "#FFFACD", "#FCAA00", "#C49102", "#FF5349")

p<-ggplot(tf[tf$dropout=="no" & (tf$group=="SCCHN" | tf$group=="OMD"),],
          aes(x=treatment, y=tumor_fraction, group=patient)) +
  geom_line(aes(color=patient), linewidth=1.5) + 
  scale_color_manual(values = c(greens, yellows)) +
  geom_hline(yintercept = 0.02, linetype = 'dashed', col = 'red') +
  annotate("text", x=5.5, y=0.0257, label="2%", size = 3) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(.~group, scales = "free", space = "free_x")+
  scale_y_continuous(labels = scales::percent) + ylab("Tumor fraction") +
  xlab("Timepoint")
p

### Compare size selected vs non-size-selected
p<-ggplot(tflong[tflong$dropout=="no",], aes(x=sample, y=tumor_fraction, fill=size_selection)) +
  geom_bar(stat="identity", position="dodge")+ 
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("darkgreen", "darkviolet")) +
  facet_grid(.~group, scales = "free", space = "free_x")
p

p<-ggplot(tflong[tflong$dropout=="no" & tflong$day=="0",], aes(x=sample, y=tumor_fraction, fill=size_selection)) +
  geom_bar(stat="identity", position="dodge")+ 
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("darkgreen", "darkviolet")) +
  facet_grid(.~group, scales = "free", space = "free_x")
p

p<-ggplot(tflong[tflong$dropout=="no",], aes(x=group, y=tumor_fraction, fill=size_selection)) +
  geom_boxplot() +
  scale_fill_manual(values = c("darkgreen", "darkviolet")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p<-ggplot(tflong[tflong$dropout=="no" & tflong$day=="0",], aes(x=group, y=tumor_fraction, fill=size_selection)) +
  geom_boxplot() +
  scale_fill_manual(values = c("darkgreen", "darkviolet")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p<-ggplot(tflong[tflong$dropout=="no",], aes(x=group, y=tumor_fraction, fill=group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(.~size_selection, scales = "free", space = "free_x") +
  scale_y_continuous(labels = scales::percent) + ylab("Tumor fraction")
p

p<-ggplot(tflong[tflong$dropout=="no" & tflong$day=="0",], aes(x=group, y=tumor_fraction, fill=group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(.~size_selection, scales = "free", space = "free_x")
p

tf$patient <- factor(tf$patient, levels=c('PV5', 'PV3', 'PV7', 'PV6', 'PV9', 'PV1', 'PV4'))
p<-ggplot(tf[tf$dropout=="no" & tf$group=="PMD",], aes(x=patient, y=tumor_fraction)) +
  geom_bar(stat="identity", fill="red")+ 
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent)
p
tf$`death after sampling`[tf$patient=="PV4"] <- 730
ggplot(tf[tf$dropout=="no" & tf$group=="PMD",], aes(x=`death after sampling`, y=tumor_fraction)) +
  geom_point(color="red") + theme_classic() + 
  geom_text(label=tf$patient[tf$dropout=="no" & tf$group=="PMD"])+
  scale_y_continuous(labels = scales::percent) + ylab("Tumor fraction") +
  xlab("Days of survival after sampling")
cor(tf$tumor_fraction[tf$dropout=="no" & tf$group=="PMD"], tf$`death after sampling`[tf$dropout=="no" & tf$group=="PMD"])

################################################################################
######## Tumor fraction data with Cristiano #########

cristiano_samplemap <- read_tsv("data/Cristiano_samplemap.tsv")
unique(cristiano_samplemap$`Patient Type`)
cristiano_samplemap$Stage[is.na(cristiano_samplemap$Stage)] <-"NA"

p <- ggplot(cristiano_samplemap[cristiano_samplemap$`Sample Type`== "cfDNA" &
                                  cristiano_samplemap$`Patient Type`!= "Duodenal Cancer" ,],
            aes(x=`Stage`, y=tfx, fill=`Patient Type`)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("green", "pink", "darkblue",
                               "lavender", "grey", "yellow", "darkcyan", "purple")) +
  theme_classic() + facet_grid(.~`Patient Type`, scales = "free_x", space = "free")
p

p <- ggplot(cristiano_samplemap[cristiano_samplemap$`Sample Type`== "cfDNA" &
                                  cristiano_samplemap$`Patient Type`!= "Duodenal Cancer" ,],
            aes(x=`Stage`, y=tfx, fill=Stage)) + 
  geom_violin(scale = "count", adjust=0.5, alpha=0.75) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.005, stackratio = 0.5) + 
  scale_fill_manual(values = c("#e7d87d", "#dd9f40", "#b4451f",
                               "#b01111", "lavender", "grey")) +
  theme_classic() + facet_grid(.~`Patient Type`, scales = "free_x", space = "free") 
p 

p <- ggplot(cristiano_samplemap[cristiano_samplemap$`Sample Type`== "cfDNA" &
                                  cristiano_samplemap$`Patient Type`!= "Duodenal Cancer" ,],
            aes(x=`Stage`, y=`cfDNA Extracted (ng/ml)`, fill=Stage)) + 
  geom_violin(scale = "count", adjust=0.5, alpha=0.75) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 2, stackratio = 0.5) + 
  scale_fill_manual(values = c("#e7d87d", "#dd9f40", "#b4451f",
                               "#b01111", "lavender", "grey")) +
  theme_classic() + facet_grid(.~`Patient Type`, scales = "free_x", space = "free") 
p 


p <- ggplot(cristiano_samplemap[cristiano_samplemap$`Sample Type`== "cfDNA" &
                                  cristiano_samplemap$`Patient Type`!= "Duodenal Cancer" &
                                  cristiano_samplemap$Stage!="X",],
            aes(x=`Stage`, y=tfx, fill=Stage)) + 
  geom_violin(scale = "width", adjust=0.5, alpha=0.75) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.003, stackratio = 0.5) + 
  scale_fill_manual(values = c("#e7d87d", "#dd9f40", "#b4451f",
                               "#b01111", "lavender", "grey")) +
  theme_classic()
p

View(cristiano_samplemap[cristiano_samplemap$`Sample Type`== "cfDNA" &
                           cristiano_samplemap$`Patient Type`!= "Duodenal Cancer" &
                           cristiano_samplemap$Stage!="X",])

+ scale_y_continuous(trans = 'log10')
################################################################################
######## length data #########
count <- read_tsv("results/uniqcounts_100_450.tsv")
count[is.na(count)] <- 0
ratio <- count %>% mutate_at(vars(-length), funs(./sum(.))) 

ratiolong <- melt(ratio, 
                  id.vars=c("length"),
                  variable.name="sample",
                  value.name="ratio")


rat150 <- as.data.frame(colSums(ratio[ratio$length<=150 & ratio$length>=100,][,-1]))
colnames(rat150) <- "ratio150"
rat150$sample <- rownames(rat150)
rat150samp <- merge(x = rat150, y = samplemap, by = "sample", all.x = TRUE)

p<-ggplot(rat150samp[rat150samp$dropout== "no" & rat150samp$day=="0",]
          , aes(x=sample, y=ratio150, fill=group)) +
  geom_bar(stat="identity")+ 
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p <- ggplot(rat150samp[rat150samp$dropout== "no" & rat150samp$day=="0",],
            aes(x=group, y=ratio150, fill=group)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic()
p

p <- ggplot(rat150samp, aes(x=group, y=ratio150, fill=group)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic()
p

p <- ggplot(rat150samp, aes(x=treatment, y=ratio150, fill=group)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + facet_grid(.~group, scales = "free", space = "free")
p

tf_rat150 <- merge(x = rat150samp, y = tf[, c("sample", "tumor_fraction")], by = "sample", all.x = TRUE)

ggplot(tf_rat150[tf_rat150$dropout=="no",], aes(x=tumor_fraction, y=ratio150)) +
  geom_point(aes(color=group)) + theme_classic() +
  scale_color_manual(values = c("grey", "green", "yellow", "red")) +
  scale_x_continuous(labels = scales::percent) + xlab("Tumor fraction") +
  ylab("Frequency of fragments of length (100-150 bp)") +
  geom_vline(xintercept = 0.02, linetype = 'dashed', col = 'red') +
  annotate("text", x=0.026, y=0.02, label="2%", size = 3)
cor(tf_rat150$ratio150[tf_rat150$dropout=="no"], tf_rat150$tumor_fraction[tf_rat150$dropout=="no"])


rat300 <- as.data.frame(colSums(ratio[ratio$length<=325 & ratio$length>=275,][,-1]))
colnames(rat300) <- "ratio300"
rat300$sample <- rownames(rat300)
rat300samp <- merge(x = rat300, y = samplemap, by = "sample", all.x = TRUE)

p<-ggplot(rat300samp[rat300samp$dropout== "no" & rat300samp$day=="0",],
          aes(x=sample, y=ratio300, fill=group)) +
  geom_bar(stat="identity")+ 
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p <- ggplot(rat300samp[rat300samp$dropout== "no" & rat300samp$day=="0",],
            aes(x=group, y=ratio300, fill=group)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic()
p

p<-ggplot(rat300samp[rat300samp$dropout== "no",], aes(x=sample, y=ratio300, fill=group)) +
  geom_bar(stat="identity")+ 
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p <- ggplot(rat300samp[rat300samp$dropout== "no",],
            aes(x=group, y=ratio300, fill=group)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic()
p

rat230 <- as.data.frame(colSums(ratio[ratio$length<=240 & ratio$length>=190,][,-1]))
colnames(rat230) <- "ratio230"
rat230$sample <- rownames(rat230)
rat230samp <- merge(x = rat230, y = samplemap, by = "sample", all.x = TRUE)

p<-ggplot(rat230samp[rat230samp$dropout== "no" & rat230samp$day=="0",],
          aes(x=sample, y=ratio230, fill=group)) +
  geom_bar(stat="identity")+ 
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p <- ggplot(rat230samp[rat230samp$dropout== "no" & rat230samp$day=="0",],
            aes(x=group, y=ratio230, fill=group)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic()
p

p<-ggplot(rat230samp[rat230samp$dropout== "no",], aes(x=sample, y=ratio230, fill=group)) +
  geom_bar(stat="identity")+ 
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p <- ggplot(rat230samp[rat230samp$dropout== "no",],
            aes(x=group, y=ratio230, fill=group)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic()
p

p <- ggplot(rat230samp[rat230samp$dropout== "no",],
            aes(x=treatment, y=ratio230, fill=group)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + facet_grid(. ~ group, scales = "free_x", space = "free_x")
p

rat230samp$grouptreatment <- paste(rat230samp$group, rat230samp$treatment, sep= "_")

anov <- compare_means(ratio230 ~ grouptreatment,
                      data = rat230samp[rat230samp$dropout == "no",],
                      method = "anova")

pwise <- compare_means(ratio230 ~ grouptreatment,
                       data = rat230samp[rat230samp$dropout == "no" & rat230samp$treatment!="P1",],
                       ref.group = "Healthy_BL",
                       method = "t.test")

omdanov <- compare_means(ratio230 ~ grouptreatment,
                         data = rat230samp[rat230samp$dropout == "no" & rat230samp$group=="OMD",],
                         method = "anova")

omdpwise <- compare_means(ratio230 ~ treatment,
                       data = rat230samp[rat230samp$dropout == "no" & rat230samp$group=="OMD",],
                       ref.group = "BL",
                       method = "t.test")



rat230samp$grouptreatment <-factor(rat230samp$grouptreatment, levels=c("Healthy_BL",
                                                                       "SCCHN_BL", "SCCHN_T1", "SCCHN_T2", "SCCHN_F1", "SCCHN_F2", "SCCHN_F3", "SCCHN_P1",
                                                                       "OMD_BL", "OMD_T1", "OMD_T2", "OMD_F1", "OMD_F2", "OMD_F3", "OMD", 
                                                                       "PMD_BL"))

p <- ggplot() +
  geom_boxplot(data=rat230samp[rat230samp$dropout == "no",], aes(x=grouptreatment, 
                                                               y=ratio230,
                                                               fill = group)) +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  stat_compare_means(method = "anova") +
  theme_classic() + 
  stat_pvalue_manual(pwise, label = "p.signif", y.position = 0.26, step.increase = 0.035) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  geom_text(aes(label = paste0("ANOVA, p=", anov$p.adj), x=2.5, y=0.43)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Ratio of fragments 190-240bp")
p

p <- ggplot(rat230samp[rat230samp$dropout== "no",],
            aes(x=treatment, y=ratio230, fill=group)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() 
p

################################################################################
####### Non-negative matrix factorization ########

ratio <- as.data.frame(ratio)
rownames(ratio) <- ratio$length
ratio$length <- NULL

# remove dropouts
ratio <- ratio[, samplemap$sample[samplemap$dropout=="no"]]

# Take only baseline
#ratio <- ratio[, samplemap$sample[samplemap$dropout=="no" & samplemap$day=="0" & samplemap$group!="Healthy"]]

# Deciding on the number of factors

ratiosvd <- svd(ratio)
ratiosvd$d

estim.r <- nmf(ratio, 2:7, nrun=10, seed=123)

opar <- par(mfrow=c(3,3))
consensusmap(estim.r, annCol=NA, labCol=NA, labRow=NA)
par(opar)
plot(estim.r)


## Running NMF

res <- nmf(ratio, 4, nrun=100, seed=123)

fit(res)

V.hat <- fitted(res)
dim(V.hat)

summary(res)

w <- basis(res) #  the values of each signature along the length
dim(w)

h <- coef(res) # the weights of each signature in the samples
dim(h) #  r x p (r = 3 p = 100)

## Plotting ###

w <- as.data.frame(w)
w <- w %>% dplyr::rename(sig1 = V1, sig2 = V2, sig3 = V3, sig4 = V4)
w$length <- rownames(w)

signlong <- melt(w, 
                 id.vars=c("length"),
                 variable.name="signature",
                 value.name="amplitude")

signlong$length <- as.numeric(signlong$length)

p<-ggplot(signlong, aes(x=length, y=amplitude, group=signature)) +
  geom_line(aes(color=signature)) + theme_classic() +
  scale_color_manual(values = c("red", "orange", "black", "blue", "lightslateblue", "grey"))
p
ggplotly(p)

h <- as.data.frame(h)
h <- sweep(h,2,colSums(h),`/`)
h[5,] <- h[2,] + h[3,]

h$signature <- as.vector(c("sig1", "sig2"
                           , "sig3"
                           , "sig4"
                           , "sig2plussig3"
                           #                           , "sig5"
#                           , "sig6"
))
#h[5,-94] <- h[2,-94] + h[3,-94]
#h[5,94] <- "sig2plussig3"

signpersample <- melt(h, 
                      id.vars=c("signature"),
                      variable.name="sample",
                      value.name="amplitude")

signsamplemerged <- merge(x = signpersample, y = samplemap, by = "sample", all.x = TRUE)

p<-ggplot(signsamplemerged, aes(x=sample, y=amplitude, fill=signature)) +
  geom_bar(stat="identity")+ 
  scale_fill_manual(values = c("lightslateblue", "orange", "black", "purple", "blue", "grey")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p<-ggplot(signsamplemerged[signsamplemerged$signature=="sig1",], aes(x=sample, y=amplitude, fill=signature)) +
  geom_bar(stat="identity")+ 
  scale_fill_manual(values = c("lightslateblue")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p<-ggplot(signsamplemerged[signsamplemerged$signature=="sig2",], aes(x=sample, y=amplitude, fill=signature)) +
  geom_bar(stat="identity")+ 
  scale_fill_manual(values = c("orange")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p<-ggplot(signsamplemerged[signsamplemerged$signature=="sig3",], aes(x=sample, y=amplitude, fill=signature)) +
  geom_bar(stat="identity")+ 
  scale_fill_manual(values = c("black")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p<-ggplot(signsamplemerged[signsamplemerged$signature=="sig4",], aes(x=sample, y=amplitude, fill=signature)) +
  geom_bar(stat="identity")+ 
  scale_fill_manual(values = c("purple")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p<-ggplot(signsamplemerged[signsamplemerged$signature=="sig4" & signsamplemerged$group!="P1",],
          aes(x=sample, y=amplitude, fill=signature)) +
  geom_bar(stat="identity")+ 
  scale_fill_manual(values = c("purple")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

#p<-ggplot(signsamplemerged[signsamplemerged$signature=="sig5",], aes(x=sample, y=amplitude, fill=signature)) +
#  geom_bar(stat="identity")+ 
#  scale_fill_manual(values = c("blue")) +
#  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p

#p<-ggplot(signsamplemerged[signsamplemerged$signature=="sig6",], aes(x=sample, y=amplitude, fill=signature)) +
#  geom_bar(stat="identity")+ 
#  scale_fill_manual(values = c("grey")) +
#  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p

p<-ggplot(signsamplemerged, aes(x=treatment, y=amplitude, fill=signature)) +
  geom_boxplot()+ 
  scale_fill_manual(values = c("lightslateblue", "orange", "black", "purple", "blue", "black")) +
  theme_classic() +
  facet_grid(. ~ group, scales = "free_x", space = "free_x")
p

p<-ggplot(signsamplemerged[signsamplemerged$signature=="sig1",], aes(x=treatment, y=amplitude, fill=group)) +
  geom_boxplot()+ 
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p<-ggplot(signsamplemerged[signsamplemerged$signature=="sig2",], aes(x=treatment, y=amplitude, fill=group)) +
  geom_boxplot()+ 
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p<-ggplot(signsamplemerged[signsamplemerged$signature=="sig3",], aes(x=treatment, y=amplitude, fill=group)) +
  geom_boxplot()+  
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p<-ggplot(signsamplemerged[signsamplemerged$signature=="sig4",], aes(x=treatment, y=amplitude, fill=group)) +
  geom_boxplot()+ 
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p<-ggplot(signsamplemerged[signsamplemerged$signature=="sig4" & signsamplemerged$treatment!="P1",],
          aes(x=treatment, y=amplitude, fill=group)) +
  geom_boxplot()+ 
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p<-ggplot(signsamplemerged[signsamplemerged$signature=="sig2plussig3" & signsamplemerged$treatment!="P1",],
          aes(x=treatment, y=amplitude, fill=group)) +
  geom_boxplot()+ 
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

#p<-ggplot(signsamplemerged[signsamplemerged$signature=="sig5",], aes(x=treatment, y=amplitude, fill=group)) +
#  geom_boxplot()+ 
#  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
#  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p

#p<-ggplot(signsamplemerged[signsamplemerged$signature=="sig6",], aes(x=treatment, y=amplitude, fill=group)) +
#  geom_boxplot()+ 
#  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
#  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p


### Plotting the factors with the cfdna length distributions
signlong$amplitude_scaled <- scale(signlong$amplitude)
signlong$amplitude_scaled <- (signlong$amplitude-mean(signlong$amplitude)) /sd(signlong$amplitude)
w_scaled <- as_data_frame(scale(w[, 1:4]))
w_scaled$length <- w$length


sign_z_scaledlong <- melt(w_scaled, 
                          id.vars=c("length"),
                          variable.name="signature",
                          value.name="amplitude")

sign_z_scaledlong$length <- as.numeric(sign_z_scaledlong$length)

p<-ggplot(sign_z_scaledlong, aes(x=length, y=amplitude, group=signature)) +
  geom_line(aes(color=signature)) + theme_classic() +
  scale_color_manual(values = c("lightslateblue", "orange", "black", "purple", "blue", "grey"))
p

ggplot(sign_z_scaledlong, aes(x = length, y = signature , fill = amplitude)) +
  geom_tile() +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") + theme_classic()

meanw <- rowMeans(w[, 1:4])
w_log_scaled <- as_data_frame(scale(w[, 1:4]/meanw))
w_log_scaled$length <- w$length

sign_z_scaledlong <- melt(w_log_scaled, 
                          id.vars=c("length"),
                          variable.name="signature",
                          value.name="amplitude")

sign_z_scaledlong$length <- as.numeric(sign_z_scaledlong$length)

p<-ggplot(sign_z_scaledlong, aes(x=length, y=amplitude, group=signature)) +
  geom_line(aes(color=signature)) + theme_classic() +
  scale_color_manual(values = c("lightslateblue", "orange", "black", "purple", "blue", "grey"))
p

ggplot(sign_z_scaledlong, aes(x = length, y = signature , fill = amplitude)) +
  geom_tile() +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") + theme_classic()

pmds <- (samplemap$sample[samplemap$dropout=="no" & samplemap$group=="PMD"])
pmdmean <- apply(ratio[, (colnames(ratio) %in% pmds)], 1, mean, na.rm=T)
healthies <- (samplemap$sample[samplemap$dropout=="no" & samplemap$group=="Healthy"])
healthymean <- apply(ratio[, (colnames(ratio) %in% healthies)], 1, mean, na.rm=T)
omds <- (samplemap$sample[samplemap$dropout=="no" & samplemap$group=="OMD" &
                            samplemap$treatment=="BL"])
omdmean <- apply(ratio[, (colnames(ratio) %in% omds)], 1, mean, na.rm=T)
hns <- (samplemap$sample[samplemap$dropout=="no" & samplemap$group=="SCCHN" &
                           samplemap$treatment=="BL"])
hnmean <- apply(ratio[, (colnames(ratio) %in% hns)], 1, mean, na.rm=T)

means <- data.frame(Healthy=healthymean, SCCHN=hnmean, OMD=omdmean, PMD=pmdmean,
                    length=rownames(ratio))

meanslong <- melt(means, 
                  id.vars=c("length"),
                  variable.name="sample",
                  value.name="frequency")
meanslong$length <- as.numeric(meanslong$length)
p<-ggplot(meanslong, aes(x=length, y=frequency, group=sample)) +
  geom_line(aes(color=sample)) + theme_classic() +
  scale_color_manual(values = c("grey", "green", "yellow", "red"))
p

selected <- data.frame(length=rownames(ratio), Healthy=healthymean, OMD9=ratio$`OMD9-BL`,
                       PV3=ratio$PV3, PV9=ratio$PV9)

selectedlong <- melt(selected, 
                     id.vars=c("length"),
                     variable.name="sample",
                     value.name="frequency")
selectedlong$length <- as.numeric(selectedlong$length)

p<-ggplot(selectedlong, aes(x=length, y=frequency, group=sample)) +
  geom_line(aes(color=sample)) + theme_classic() +
  scale_color_manual(values = c("grey", "yellow", "red", "orange"))
p

#pmds <- (samplemap$sample[samplemap$dropout=="no" & samplemap$group=="PMD"])
#pmdmean <- apply(ratio[, (colnames(ratio) %in% pmds)], 1, mean, na.rm=T)
#healthies <- (samplemap$sample[samplemap$dropout=="no" & samplemap$group=="Healthy"])
#healthymean <- apply(ratio[, (colnames(ratio) %in% healthies)], 1, mean, na.rm=T)
#omds <- (samplemap$sample[samplemap$dropout=="no" & samplemap$group=="OMD" &
#                            samplemap$treatment=="BL"])
#omdmean <- apply(ratio[, (colnames(ratio) %in% omds)], 1, mean, na.rm=T)
omdt1 <- (samplemap$sample[samplemap$dropout=="no" & samplemap$group=="OMD" &
                             samplemap$treatment=="T1"])
omdt1mean <- apply(ratio[, (colnames(ratio) %in% omdt1)], 1, mean, na.rm=T)
omdt2 <- (samplemap$sample[samplemap$dropout=="no" & samplemap$group=="OMD" &
                             samplemap$treatment=="T2"])
omdt2mean <- apply(ratio[, (colnames(ratio) %in% omdt2)], 1, mean, na.rm=T)

omdf1 <- (samplemap$sample[samplemap$dropout=="no" & samplemap$group=="OMD" &
                             samplemap$treatment=="F1"])
omdf1mean <- apply(ratio[, (colnames(ratio) %in% omdf1)], 1, mean, na.rm=T)

omdf2 <- (samplemap$sample[samplemap$dropout=="no" & samplemap$group=="OMD" &
                             samplemap$treatment=="F2"])
omdf2mean <- apply(ratio[, (colnames(ratio) %in% omdf2)], 1, mean, na.rm=T)

omdf3 <- (samplemap$sample[samplemap$dropout=="no" & samplemap$group=="OMD" &
                             samplemap$treatment=="F3"])
omdf3mean <- apply(ratio[, (colnames(ratio) %in% omdf3)], 1, mean, na.rm=T)

radselected <- data.frame(length=rownames(ratio), Healthy=healthymean, T1=omdt1mean,
                          T2=omdt2mean, F1=omdf1mean, F2=omdf2mean, F3=omdf3mean)

selectedlong <- melt(radselected, 
                     id.vars=c("length"),
                     variable.name="sample",
                     value.name="frequency")
selectedlong$length <- as.numeric(selectedlong$length)

p<-ggplot(selectedlong, aes(x=length, y=frequency, group=sample)) +
  geom_line(aes(color=sample)) + theme_classic() +
  scale_color_manual(values = c("grey", "#FFD400", "orange", "#C89F5D", "#BF8801", "#692C00"))
p

spswide <- as.data.frame(t(h))
spswide$sample <- rownames(spswide)
colnames(spswide) <- spswide["signature",]

spswide$sample <- rownames(spswide)

spswide$sample <- rownames(spswide)

spswide[,c("sample", "sig1")]


tf_sig1 <- merge(x = tf, y = spswide[, c("sample", "sig1")], by = "sample", all.x = TRUE)
tf_sig1$sig1 <- as.numeric(tf_sig1$sig1)

samplelist <- c("PV9", "PV3",
                "OMD9-d000", "OMD9-d002", "OMD9-d005", "OMD9-d090", "OMD9-d180", "OMD9-365")

tf_sig1$plotname <- NA
tf_sig1$plotname[tf_sig1$tumor_fraction>0.02 & tf_sig1$sig1>0.1] <- tf_sig1$sample[tf_sig1$tumor_fraction>0.02 & tf_sig1$sig1>0.1]

ggplot(tf_sig1[tf_sig1$dropout=="no",], aes(x=tumor_fraction, y=sig1)) +
  geom_point(aes(color=group)) + theme_classic() +
  scale_color_manual(values = c("grey", "green", "yellow", "red")) +
  geom_text(aes(x = tumor_fraction+0.002, label = plotname), hjust = 0)+
  scale_x_continuous(labels = scales::percent) + xlab("Tumor fraction") +
  ylab("Signature 1") + geom_vline(xintercept = 0.02, linetype = 'dashed', col = 'red') +
  annotate("text", x=0.026, y=0.02, label="2%", size = 3)


################################################################################
####### Fragment length in amplified region 8q ########

ampcount <- read_tsv("results/uniqcounts_all.tsv")

ampcount[is.na(ampcount)] <- 0
ampratio <- ampcount %>% mutate_at(vars(-length), funs(./sum(.))) 
ampsamplemap <- read_tsv("data/samplemap_8q.tsv")

ampratiolong <- melt(ampratio, 
                  id.vars=c("length"),
                  variable.name="sample",
                  value.name="ratio")

ampratiolong <- merge(x = ampratiolong, y = ampsamplemap, by = "sample", all.x = TRUE)

### Fragment length distribution line plots
p<-ggplot(ampratiolong[ampratiolong$group=="OMD" ,], aes(x=length, y=ratio, group=sample)) +
  geom_line(aes(color=selected)) + xlim(100,150) + ylim(0, 0.0075)  + theme_classic() +
  scale_color_manual(values = c("red", "grey", "yellow"))
p

p<-ggplot(ampratiolong[ampratiolong$group=="OMD" ,], aes(x=length, y=ratio, group=sample)) +
  geom_line(aes(color=selected)) + xlim(200, 450) + ylim(0, 0.0075) + theme_classic() +
  scale_color_manual(values = c("red", "grey", "yellow"))
p

p<-ggplot(ampratiolong, aes(x=length, y=ratio, group=sample)) +
  geom_line(aes(color=selected)) + theme_classic() +
  scale_color_manual(values = c("red", "grey", "yellow"))
p

ggplotly(p)

### Enriched fragment lengths in 8q (amplified in OMD9)

#new_df <- ampratiolong %>%
#  separate(sample, into = c("sample_name", "sample_index"), sep = "_") %>% # separate sample column into samplename and sample_index columns
#  pivot_wider(id_cols = c(sample_name, length, group, day, dropout, treatment), names_from = sample_index, values_from = ratio)

#new_df <- ampratiolong %>%
#  separate(sample, into = c("sample_name", "sample_index"), sep = "_") %>% # separate sample column into samplename and sample_index columns
#  pivot_wider(id_cols = c(sample_name, length, group, day, dropout, treatment), names_from = sample_index, values_from = ratio) %>% # spread samplename to columns
#  mutate(ratio = `8q` / all)# calculate new ratio column

#new_df$log2ratio <- base::log2(new_df$ratio)
rat8q_df <- ampratiolong[ampratiolong$selected=="all",]
df_8q <-ampratiolong[ampratiolong$selected=="8q",]
rat8q_df$rat8q <- df_8q$ratio/rat8q_df$ratio
rat8q_df$log2ratio <- base::log2(rat8q_df$rat8q)

ggplot(rat8q_df[rat8q_df$dropout=="no" & rat8q_df$length>99 & rat8q_df$length<401 & 
                (rat8q_df$patient=="OMD9" | rat8q_df$group=="Healthy"),],
       aes(x = length, y = ID , fill = log2ratio)) +
  geom_tile() +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") + theme_classic()


################################################################################
####### Differential chromatin accessibility using LIQUORICE #######

celltype <- read_tsv("results/celltype_signatures_pruned.tsv")
ct_wide <- spread(celltype, key = celltype, value = signature)
ct_wide[, c(2:9)] <- scale(ct_wide[, c(2:9)])
celltype$signature <- scale(celltype$signature)
celltype <- pivot_longer(ct_wide, colnames(ct_wide[2:9]), 
                         names_to = "cell_type",
                         values_to = "signature")

celltype_merged <- merge(x=celltype , y=samplemap, by="sample", all.x=TRUE)
ct_wide_merged <- merge(x=ct_wide, y=samplemap, by="sample", all.x=TRUE)

heatmap_plot <- ggplot(celltype_merged[celltype_merged$dropout=="no",],
                       aes(x = cell_type, y = sample, fill = signature)) +
  geom_tile() +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
  theme_classic() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        strip.text.y = element_text(size = 6),
        legend.key.size = unit(0.3, 'cm'),
        panel.spacing = unit(0.03, "lines")
        #legend.key.size = unit(0.3, 'cm'),
        #legend.spacing.x = unit(0.03, 'cm')
  ) +
  facet_grid(group~., scales = "free_y", space = "free_y",
             labeller = label_wrap_gen(width = 5))
heatmap_plot


heatmap_plot <- ggplot(celltype_merged[celltype_merged$dropout=="no" &
                                         celltype_merged$treatment=="BL",],
                       aes(x = cell_type, y = sample, fill = signature)) +
  geom_tile() +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
  theme_classic() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        strip.text.y = element_text(size = 6),
        legend.key.size = unit(0.3, 'cm'),
        panel.spacing = unit(0.03, "lines")
        #legend.key.size = unit(0.3, 'cm'),
        #legend.spacing.x = unit(0.03, 'cm')
  ) +
  facet_grid(group~., scales = "free_y", space = "free_y",
             labeller = label_wrap_gen(width = 5))
heatmap_plot

anov <- compare_means(signature ~ group,
                      data = celltype_merged[celltype_merged$dropout=="no" &
                                               celltype_merged$treatment=="BL"& 
                                               celltype_merged$cell_type=="hematopoietic",],
                      method = "anova")

pwise <- compare_means(signature ~ group,
                       data = celltype_merged[celltype_merged$dropout=="no" &
                                                celltype_merged$treatment=="BL"& 
                                                celltype_merged$cell_type=="hematopoietic",],
                       ref.group = "Healthy",
                       method = "t.test")

p<-ggplot(celltype_merged[celltype_merged$dropout=="no" &
                            celltype_merged$treatment=="BL"& 
                            celltype_merged$cell_type=="hematopoietic",],
          aes(x=group, y=signature)) +
  geom_boxplot(aes(fill=group))+  
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  stat_pvalue_manual(pwise, label = "p.signif", y.position = 0.7, step.increase = 0.04) +
  geom_text(aes(label = paste0("ANOVA, p=", anov$p.adj), x=1, y=-1.2)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

anov <- compare_means(signature ~ group,
                      data = celltype_merged[celltype_merged$dropout=="no" &
                                               celltype_merged$treatment=="BL"& 
                                               celltype_merged$cell_type=="hepatocyte",],
                      method = "anova")

pwise <- compare_means(signature ~ group,
                       data = celltype_merged[celltype_merged$dropout=="no" &
                                                celltype_merged$treatment=="BL"& 
                                                celltype_merged$cell_type=="hepatocyte",],
                       ref.group = "Healthy",
                       method = "t.test")

p<-ggplot(celltype_merged[celltype_merged$dropout=="no" &
                            celltype_merged$treatment=="BL"& 
                            celltype_merged$cell_type=="hepatocyte",],
          aes(x=group, y=signature)) +
  geom_boxplot(aes(fill=group))+  
  scale_fill_manual(values = c("grey", "green", "yellow", "red")) +
  stat_pvalue_manual(pwise, label = "p.signif", y.position = 0.6, step.increase = 0.04) +
  geom_text(aes(label = paste0("ANOVA, p=", anov$p.adj), x=1, y=-1.0)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

p<-ggplot(celltype_merged[celltype_merged$dropout=="no" &
                            (celltype_merged$group=="Healthy" |
                              celltype_merged$patient=="OMD9" |
                              celltype_merged$patient=="PV9") & 
                            celltype_merged$cell_type=="prostate",],
          aes(x=sample, y=signature, fill=group)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = c("grey", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p


p<-ggplot(celltype_merged[celltype_merged$dropout=="no" &
                            (celltype_merged$group=="Healthy" |
                               celltype_merged$patient=="OMD9") & 
                            celltype_merged$cell_type=="prostate",],
          aes(x=sample, y=signature, fill=group)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = c("grey", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p


p<-ggplot(celltype_merged[celltype_merged$dropout=="no" &
                            celltype_merged$patient=="OMD9" & 
                            celltype_merged$cell_type=="prostate",],
          aes(x=sample, y=signature, fill=group)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = c("yellow")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p


p<-ggplot(celltype_merged[celltype_merged$dropout=="no" &
                            (celltype_merged$group=="Healthy" |
                               celltype_merged$patient=="OMD9") & 
                            celltype_merged$cell_type=="hematopoietic",],
          aes(x=sample, y=signature, fill=group)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = c("grey", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p


p<-ggplot(celltype_merged[celltype_merged$dropout=="no" &
                            (celltype_merged$group=="Healthy" |
                               celltype_merged$patient=="OMD9" |
                               celltype_merged$patient=="PV9") & 
                            celltype_merged$cell_type=="hematopoietic",],
          aes(x=sample, y=signature, fill=group)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = c("grey", "yellow", "red")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p


# OMD9-tumor load
days <- c(0, 2, 5, 90, 180, 365)
prostate_signature <- celltype_merged$signature[celltype_merged$dropout=="no" &
                                                  celltype_merged$patient=="OMD9" & 
                                                  celltype_merged$cell_type=="prostate"]
tumor_fraction <- tf$tumor_fraction[tf$patient=="OMD9"]
omd9 <- data.frame(days, prostate_signature, tumor_fraction)
omd9$days <- as.factor(omd9$days)
#coeff <- mean(omd9$tumor_fraction)/mean(omd9$prostate_signature)

ggplot(omd9, aes(x=days)) +
  geom_line( aes(y=tumor_fraction, group = 1), lwd=1.1, color="darkgreen") + 
  scale_y_continuous(name = "Tumor fraction")+ theme_classic()+
  scale_y_continuous(labels = scales::percent) + ylab("Tumor fraction") +
  xlab("Days after starting RT")

ggplot(omd9, aes(x=days, y=prostate_signature, group = 1)) +
  geom_line(lwd=1.1, color="purple") +
  scale_y_continuous(name = "Nucleosome accessibility at prostate-specific promoters") + 
  theme_classic() + xlab("Days after starting RT")

################################################################################
####### OMD9 PSA values ###########
psa <- read_tsv("results/OMD9_PSA.tsv")
psa$days <- psa$months*30
psa$months <- as.factor(psa$months)
p <- ggplot() +
  geom_bar(data=psa, aes(x=months, y=PSA, fill="OMD9"),
           stat = "identity") + scale_fill_manual(values = c("yellow")) +
  theme_classic() + ylab("PSA concentration ng/ml")+ xlab("months after RT") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
p
psa$patient <- "OMD9"
p <- ggplot() +
  geom_line(data=psa, aes(x=months, y=PSA, group=patient), color= "darkgrey", lwd=1.1) +
  theme_classic() + ylab("PSA concentration ng/ml")+ xlab("months after RT") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
p

p <- ggplot() +
  geom_bar(data=psa, aes(x=days, y=PSA, fill="OMD9"),
           stat = "identity") + scale_fill_manual(values = c("yellow")) +
  theme_classic() + ylab("PSA concentration ng/ml")+ xlab("days after the start of RT") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
p

## Merge PSA and prostate promoter chromatin accessibility
psaandliquorice <- ct_wide_merged[ct_wide_merged$patient=="OMD9",]
psaandliquorice <- psaandliquorice %>% rename(days = day)
psaandliquorice <- psaandliquorice[,c("prostate", "days")]
psa <- merge(x = psa, y = psaandliquorice, by = "days", all = TRUE)

a <- ggplot(psa, aes(x=days)) +
  geom_line( aes(y=prostate), linewidth=2, color="purple", stat="identity") + 
  geom_line( aes(y=PSA), linewidth=2, color="orange", stat="identity") +
  scale_y_continuous(name = "PSA concentration (ng/ml)",
    sec.axis = sec_axis(~.*0.8 ,name="Nucleosome accessibility at prostate-specific promoters (z-score)"))
a

################################################################################
####### HPV ###########

hpv <- data.frame (day  = c(0, 8, 49, 90, 0, 8, 49, 90, 0, 8, 49, 90, 0, 8, 49, 90),
                   HPV = c(0.303733717, 4.690612866, 0, 0.024438845, 2.112608682, 8.316719668, 1.190032045, 0.027690872, 2.964606464, 4.88602756, 0.087429904, 0, 0, 0.830933895, 0, 0),
                   patient = c("HN001", "HN001", "HN001", "HN001", "HN002", "HN002", "HN002", "HN002", "HN003", "HN003", "HN003", "HN003", "HN004", "HN004", "HN004", "HN004")
)

ggplot(data=hpv, aes(x=day, y=HPV, color=patient)) +
  geom_line(size=1.2) + theme_classic() +
  scale_color_brewer(palette="Dark2")

hpv$day <- as.factor(hpv$day)
ggplot(hpv, aes(x = day, y = HPV)) +
  geom_boxplot(fill = "green") + geom_jitter(shape=16, position=position_jitter(0.2)) + 
  theme_classic()


df
df <- aggregate(hpv$HPV, list(hpv$day), FUN=mean)
df<- df %>% rename(day = Group.1,
                   HPV = x)
df$sd <- aggregate(hpv$HPV, list(hpv$day), FUN=sd)$x
p<- ggplot(df, aes(x=day, y=HPV)) + 
  geom_bar(stat = "identity", fill="green2") +
  geom_errorbar(aes(ymin=HPV, ymax=HPV+sd), width=.2,
                position=position_dodge(.9)) + theme_classic()
p

q<- ggplot(df, aes(x=day, y=HPV)) + 
  geom_bar(stat = "identity", fill="green2") +
  theme_classic()
q

########################################################
###### Viral read counts ######
viral_read <- read_tsv("results/viral_reads.tsv")
viral_read <- merge(x = viral_read, y = samplemap, by = "sample", all.x = TRUE)
greens <- c("#9CCC65", "#CDDC39", "#388E3C", "#4CAF50", "#00BFA5", "#2E7D32", "#689F38")
hpvpatients <- c("HN1", "HN2", "HN3", "HN4")
p <- ggplot() +
  geom_line(data=viral_read[viral_read$patient %in% hpvpatients,],
            aes(x=as.factor(day), y=viral_reads, group=patient, color=patient, lwd=patient)) +
  theme_classic() + scale_color_manual(values = greens) +
  ylab("HPV reads / 100million human") + xlab("Day") + scale_linewidth_manual(values=c(1,1,1,1))+
  geom_hline(yintercept=126, linetype = 'dashed', col = 'grey') +
  annotate("text", x=4.9, y=129.5, label="0.5% tumor fraction at VL=100") + 
  geom_hline(yintercept=12.6, linetype = 'dashed', col = 'grey') +
  annotate("text", x=4.92, y=16, label="0.5% tumor fraction at VL=10")
p

################################################################################
###### Viral fragment lengths ######

viral <- read_tsv("results/viral_length.tsv",
                  col_names="viral")
p <- ggplot(viral, aes(x=viral)) + 
  geom_density() + theme_classic()
p


x = -10:10
y = dnorm(x, mean=0, sd=3)
df.norm = data.frame('x'=x, 'y'=y)

random = data.frame('x'=rnorm(1000, mean = 0, sd = 3))

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

hpvmode <- getmode(viral$viral)
hnmode <- ratiolong$length[ratiolong$sample=="HN3-BL"& ratiolong$ratio==max(ratiolong$ratio[ratiolong$sample=="HN3-BL"])]

hpvsamples <- c("HN001-d008", "HN2-BL", "HN2-d008", "HN2-d049", "HN3-BL", "HN3-d008", "HN4-BL")
p <- ggplot() +
  geom_line(data=ratiolong[ratiolong$sample %in% hpvsamples,], aes(x=length, y=ratio, group=sample), color = "green") +
  xlim(100,400) + theme_classic() + 
  geom_density(data=viral, aes(x=viral), size=0.8, color="#800020") + 
  geom_vline(aes(xintercept = hpvmode))+
  geom_vline(aes(xintercept = hnmode)) +
  annotate(x=hpvmode-6,y=+Inf,label=hpvmode,vjust=1,geom="label") +
  annotate(x=hnmode+6,y=+Inf,label=hnmode,vjust=1,geom="label")
p
ggplotly(p)

################################################################################
###### Viral subtypes ######

subtypes <- read_tsv("results/virus_genomes.tsv")
subtypes$match_rate <- 1-subtypes$substitution_rate

p <- ggplot(subtypes, aes(x=sample, y=`mapped length`)) + 
  geom_bar(stat = "identity", fill="green2") + theme_classic() +
  facet_grid(~genome)
p

p <- ggplot(subtypes, aes(x=sample, y=substitution_rate)) + 
  geom_bar(stat = "identity", fill="green2") + theme_classic() +
  facet_grid(~genome) + scale_y_continuous(labels = scales::percent)
p

p <- ggplot(subtypes, aes(x=sample, y=match_rate)) + 
  geom_bar(stat = "identity", fill="green2") + theme_classic() +
  facet_grid(~genome) + scale_y_continuous(labels = scales::percent)
p

p <- ggplot(subtypes, aes(x=sample, y=match_rate)) + 
  geom_point(color="green2", size=3) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~genome) + scale_y_continuous(labels = scales::percent)
p

################################################################################
###### Swimmer plot #######
library(swimplot)
cohort <- read_tsv("data/cohort_data.tsv")
#Converting dates (character) to days (01.01.1970-based)
cohort[,3:14] <- sapply(cohort[,3:14], function(x) as.Date(x, "%d.%m.%Y"))
#Setting the earliest day to day 0 for each patient
cohort[,3:14] <- t(apply(cohort[,3:14], 1, function(x) x-min(x, na.rm=T)))
cohort$group <-factor(cohort$group, levels=c("OMD", "SCCHN"))
#Split dataframe
##Creating arm dataframe
arm_df <- melt(cohort[colnames(cohort) %in% (c("patient", "group", "start_of_RT", "end_of_RT",
                             "end_of_study"))], 
                            id.vars=c("patient", "group"),
                            variable.name="arm",
                            value.name="day")

#After the melt, this became a factor, so we set it to charachter
arm_df$arm <- as.character(arm_df$arm)
#Renaming the arms to have during RT and non-RT phases for each cohort arm
arm_df$arm[arm_df$group=="SCCHN" & arm_df$arm!="end_of_RT"] <- "SCCHN"
arm_df$arm[arm_df$group=="SCCHN" & arm_df$arm=="end_of_RT"] <- "SCCHN_RT"
arm_df$arm[arm_df$group=="OMD" & arm_df$arm!="end_of_RT"] <- "OMD"
arm_df$arm[arm_df$group=="OMD" & arm_df$arm=="end_of_RT"] <- "OMD_RT"

# plot arms
arm_plot <- swimmer_plot(df=arm_df,id="patient",end="day",name_fill='arm',
                         id_order="patient", col="black",alpha=0.75,width=.8) + 
  scale_fill_manual(name="arm",values=c("OMD" ="#FFE699", "OMD_RT"="#BF9000",
                                        "SCCHN"="#A9D18E", "SCCHN_RT"="#548235"))
arm_plot

##Creating adverse event dataframe
ae_df <- melt(cohort[colnames(cohort) %in% (c("patient", "local_progression", 
                                              "oligoprogression", "polyprogression",
                                              "date_of_death"))], 
               id.vars=c("patient"),
               variable.name="event",
               value.name="day")

ae_df <- drop_na(ae_df)
ae_df$event <- as.character(ae_df$event)
ae_df$event[ae_df$event=="date_of_death"] <- "Death"
ae_df$event[ae_df$event=="local_progression"] <- "local progression"
ae_df$event[ae_df$event=="oligoprogression"] <- "oligoprogression"
ae_df$event[ae_df$event=="polyprogression"] <- "polyprogression"
ae_df$event_type <- ae_df$event
ae_df$event_type <- "Progression"
ae_df$event_type[ae_df$event=="Death"] <- "Death"

AE_plot <- arm_plot + swimmer_points(df_points=ae_df, id="patient", time="day",
                                     name_shape ="event",size=2.5, fill='white',
                                     col = "black") + 
  scale_shape_manual(name="event",values=c("local progression"=0,
                                         "oligoprogression"=5,
                                         "polyprogression"=6,
                                         "Death"=19))
AE_plot


colnames(cohort)

##Creating response dataframe
resp_df <- cohort[colnames(cohort) %in% (c("patient", "PR_start", 
                                                "PR_end", "CR_start",
                                                "CR_end", "end_of_study"))]
resp_df$resp_type <- "CR"
resp_df$resp_type[resp_df$PR_start>0] <- "PR"
resp_df$resp_start <- resp_df$CR_start
resp_df$resp_end <- resp_df$CR_end
resp_df$resp_start[resp_df$resp_type=="PR"] <- resp_df$PR_start[resp_df$resp_type=="PR"]
resp_df$resp_end[resp_df$resp_type=="PR"] <- resp_df$PR_end[resp_df$resp_type=="PR"]
resp_df$resp_cont <- 1
resp_df$resp_cont[resp_df$end_of_study>resp_df$resp_end] <- 0

resp_df[is.na(resp_df)] <- 0

Response_plot <- arm_plot + swimmer_lines(df_lines=resp_df, id="patient",
                                          start = "resp_start", end="resp_end",
                                          name_col="resp_type", size=1, position_n)
Response_plot

resp_plot <- arm_plot + swimmer_arrows(df_arrows=resp_df, id="patient",
                                      arrow_start='end_of_study', cont = "resp_cont",
                                      name_col='resp_type', type ="open", show.legend=FALSE, cex=1)
resp_plot

swimmer_arrows(df_arrows=resp_df, id='patient', arrow_start='end_of_study',
               cont = 'resp_cont',name_col='resp_type',type = "open",cex=1)

###########################

