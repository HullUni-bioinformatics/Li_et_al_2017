####START####



####Set working directory under your home directory ('R'is my home directory)
setwd("~/R")

###Filtration experiment###

###packages require in this script###
##The package(reshape) should after package(dplyr), otherwise the rename funtion in package(reshape) will not work

#plyr
library(plyr)
#libraries pakage,including ggplot2
library(tidyverse)
#cowplot
library(cowplot)
#scales
library(scales)
#reshape
library(reshape)
library(reshape2)
#gtable
library(gtable)
#data.table
library(dtplyr)
#gridExtra
library(gridExtra)
#geom_text_repel
library(ggrepel)
#Chi square indepandent
library(rcompanion)
#Dunn test for multiple comparisons after Kruskal-Wallis
library(FSA)
#PCA plotting
library(factoextra)
#Species Prediction And Diversity Estimation
library(SpadeR)

library(car)


FT_DC_original <- read.csv(file= "July2016_12s_onestep_filtration/Reanalysis_May_2017/AppendixS3_Filtration_July_2016_FT_DC.csv", header = TRUE)


FT_DC_original$Pond <-mapvalues(FT_DC_original$Pond, c("E1","E2","E3","E4"), 
                                           c("A","B","C","D"))



FT_DC <- FT_DC_original [which(FT_DC_original$SampleID != "T3-1-3" &
                              FT_DC_original$SampleID != "T4-1-3" &
                              FT_DC_original$SampleID != "T2-2-3"),]

#summaries FT and DC------------------------

FT_DC_average <- ddply(FT_DC, .(Pond,Treatment), summarize,
                  FT_mean = round(mean(Time), 2),
                  FT_sd = round(sd(Time), 2),
                  DC_mean = round(mean(DNA), 2),
                  DC_sd = round(sd(DNA), 2))

ddply(FT_DC, .(Treatment), summarize,
               FT_mean = round(mean(Time), 2),
               FT_sd = round(sd(Time), 2),
               DC_mean = round(mean(DNA), 2),
               DC_sd = round(sd(DNA), 2))

#write.csv(FT_DC_average,file = "July2016_12s_onestep_filtration/Reanalysis_May_2017/FT_DC_average.csv",row.names = FALSE)



# Fig2_filtration time ----------------------------------------------------



#The Statistic results are from one-way ANOVA as below:
###ANOVA_time###
FT_DC_E1 <- FT_DC [which(FT_DC$Pond=='A'),]
FT_DC_E2 <- FT_DC [which(FT_DC$Pond=='B'),]
FT_DC_E3 <- FT_DC [which(FT_DC$Pond=='C'),]
FT_DC_E4 <- FT_DC [which(FT_DC$Pond=='D'),]


# Shapiro-Wilk test of normality
shapiro.test(FT_DC_E1$Time)
qqnorm(FT_DC_E1$Time)

shapiro.test(FT_DC_E2$Time)

shapiro.test(FT_DC_E3$Time)

shapiro.test(FT_DC_E4$Time)
qqnorm(FT_DC_E4$Time)

#Not normal distribution, use Kruskal-Wallis one-way analysis of variance

#kruskal.test for non-parametric test https://rcompanion.org/rcompanion/d_06.html
kruskal.test(Time ~ Treatment, data = FT_DC_E1)

### Dunn test

DT_FT_E1 <- dunnTest(Time ~ Treatment, data = FT_DC_E1, method="bh") 


DT_FT_E1 <- DT_FT_E1$res

DT_FT_E1


ST_FT_E1 <-cldList(comparison = DT_FT_E1$Comparison, p.value = DT_FT_E1$P.unadj, threshold= 0.05)

ST_FT_E1$Pond <- "A"



kruskal.test(Time ~ Treatment, data = FT_DC_E2)

### Dunn test

DT_FT_E2 <- dunnTest(Time ~ Treatment, data = FT_DC_E2, method="bh") 


DT_FT_E2 <- DT_FT_E2$res

DT_FT_E2


ST_FT_E2 <-cldList(comparison = DT_FT_E2$Comparison, p.value = DT_FT_E2$P.unadj, threshold= 0.05)

ST_FT_E2$Pond <- "B"



kruskal.test(Time ~ Treatment, data = FT_DC_E3)

### Dunn test

DT_FT_E3 <- dunnTest(Time ~ Treatment, data = FT_DC_E3, method="bh") 


DT_FT_E3 <- DT_FT_E3$res

DT_FT_E3


ST_FT_E3 <-cldList(comparison = DT_FT_E3$Comparison, p.value = DT_FT_E3$P.unadj, threshold= 0.05)

ST_FT_E3$Pond <- "C"


kruskal.test(Time ~ Treatment, data = FT_DC_E4)

### Dunn test

DT_FT_E4 <- dunnTest(Time ~ Treatment, data = FT_DC_E4, method="bh") 


DT_FT_E4 <- DT_FT_E4$res

DT_FT_E4


ST_FT_E4 <-cldList(comparison = DT_FT_E4$Comparison, p.value = DT_FT_E4$P.unadj, threshold= 0.05)

ST_FT_E4$Pond <- "D"

FT_Statistic <- rbind(ST_FT_E1,ST_FT_E2,ST_FT_E3,ST_FT_E4)

FT_Statistic <- rename (FT_Statistic,c("Group"="Treatment"))

FT_Statistic$Treatment <-mapvalues(FT_Statistic$Treatment, c("T1_.45","T2_.8","T3_1.2", "T4_Sterivex", "T5_PF_.45", "TP5_PF"), 
                                c("T1_0.45","T2_0.8","T3_1.2","T4_Sterivex","T5_PF_0.45","TP5_PF"))




#Use FT_DC_original, then distinguish outlier as white

FT_boxplot <- ggplot(data=FT_DC_original,aes(x=Treatment,y=Time,fill=Outlier))+
                geom_boxplot(outlier.shape = NA,fill=NA)+#not fill the bars
                scale_fill_manual(values=c("black", "white"))+
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1,width=1, position = "dodge")+
                stat_summary(data=FT_DC, fun.y=mean, mapping = aes(group =Outlier),geom="point", shape=5, size=5)+#The mean calculate exculing outilers
                facet_wrap(~Pond)+
                labs(x="Treatment", y= "Filtration time (min)")+theme_bw()+
                theme(text=element_text(size=20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none")



#The Statistic results are from kruskal-Wallis one-way analysis of variance and  Dunn test as above:



FT_Statistic$y <- 135

#c(rep(50,6),rep(25,6),rep(50,6),rep(20,6))

FT_boxplot+geom_text(data=FT_Statistic,aes(Treatment,y, label=Letter, fill=NA),hjust = 0.5, size=7, parse = TRUE)



ggsave("Fig2_FT_boxplot.jpeg",path = "July2016_12s_onestep_filtration/Reanalysis_May_2017/Figures/", width = 12, height = 10, units = "in", dpi=500)

# FigS5_FT_treat --------------------------------------------------


FT_DC_T1 <- FT_DC [which(FT_DC$Treatment=='T1_0.45'),]
FT_DC_T2 <- FT_DC [which(FT_DC$Treatment=='T2_0.8'),]
FT_DC_T3 <- FT_DC [which(FT_DC$Treatment=='T3_1.2'),]
FT_DC_T4 <- FT_DC [which(FT_DC$Treatment=='T4_Sterivex'),]
FT_DC_T5 <- FT_DC [which(FT_DC$Treatment=='T5_PF_0.45'),]
FT_DC_TP5 <- FT_DC [which(FT_DC$Treatment=='TP5_PF'),]


kruskal.test(Time ~ Pond, data = FT_DC_T1)


DT_FT_T1 <- dunnTest(Time ~ Pond, data = FT_DC_T1, method="bh") 

DT_FT_T1 <- DT_FT_T1$res

DT_FT_T1


ST_FT_T1 <-cldList(comparison = DT_FT_T1$Comparison, p.value = DT_FT_T1$P.unadj, threshold= 0.05)

ST_FT_T1$Treatment <- "T1_0.45"

kruskal.test(Time ~ Pond, data = FT_DC_T2)

DT_FT_T2 <- dunnTest(Time ~ Pond, data = FT_DC_T2, method="bh") 

DT_FT_T2 <- DT_FT_T2$res

DT_FT_T2


ST_FT_T2 <-cldList(comparison = DT_FT_T2$Comparison, p.value = DT_FT_T2$P.unadj, threshold= 0.05)

ST_FT_T2$Treatment <- "T2_0.8"

kruskal.test(Time ~ Pond, data = FT_DC_T3)

DT_FT_T3 <- dunnTest(Time ~ Pond, data = FT_DC_T3, method="bh") 

DT_FT_T3 <- DT_FT_T3$res

DT_FT_T3


ST_FT_T3 <-cldList(comparison = DT_FT_T3$Comparison, p.value = DT_FT_T3$P.unadj, threshold= 0.05)

ST_FT_T3$Treatment <- "T3_1.2"

kruskal.test(Time ~ Pond, data = FT_DC_T4)
DT_FT_T4 <- dunnTest(Time ~ Pond, data = FT_DC_T4, method="bh") 

DT_FT_T4 <- DT_FT_T4$res

DT_FT_T4


ST_FT_T4 <-cldList(comparison = DT_FT_T4$Comparison, p.value = DT_FT_T4$P.unadj, threshold= 0.05)

ST_FT_T4$Treatment <- "T4_Sterivex"

kruskal.test(Time ~ Pond, data = FT_DC_T5)
DT_FT_T5 <- dunnTest(Time ~ Pond, data = FT_DC_T5, method="bh") 

DT_FT_T5 <- DT_FT_T5$res

DT_FT_T5


ST_FT_T5 <-cldList(comparison = DT_FT_T5$Comparison, p.value = DT_FT_T5$P.unadj, threshold= 0.05)

ST_FT_T5$Treatment <- "T5_PF_0.45"



kruskal.test(Time ~ Pond, data = FT_DC_TP5)

DT_FT_TP5 <- dunnTest(Time ~ Pond, data = FT_DC_TP5, method="bh") 


DT_FT_TP5 <- DT_FT_TP5$res

DT_FT_TP5


ST_FT_TP5 <-cldList(comparison = DT_FT_TP5$Comparison, p.value = DT_FT_TP5$P.unadj, threshold= 0.05)

ST_FT_TP5$Treatment<- "TP5_PF"


FT_Treament_ST <- rbind(ST_FT_T1,ST_FT_T2,ST_FT_T3,ST_FT_T4,ST_FT_T5,ST_FT_TP5)

FT_Treament_ST <- rename (FT_Treament_ST,c("Group"="Pond"))




FT_Treatment_boxplot <- ggplot(data=FT_DC_original,aes(x=Pond,y=Time,fill=Outlier))+
  geom_boxplot(outlier.shape = NA,fill=NA)+#not fill the bars
  scale_fill_manual(values=c("black", "white"))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1,width=1, position = "dodge")+
  stat_summary(data=FT_DC, fun.y=mean, mapping = aes(group =Outlier),geom="point", shape=5, size=5)+#The mean calculate exculing outilers
  facet_wrap(~Treatment)+
  labs(x="Pond", y= "Filtration time (min)")+theme_bw()+
  theme(text=element_text(size=20),legend.position = "none")


FT_Treament_ST$y <- 135
  
#c(rep(80,4),rep(55,4),rep(20,4),rep(110,4),rep(20,4),rep(4,4))

FT_Treatment_boxplot+geom_text(data=FT_Treament_ST,aes(Pond,y, label=Letter, fill=NA),hjust = 0.5, size=7, parse = TRUE)

ggsave("FigS5_FT_treat.jpeg",path = "July2016_12s_onestep_filtration/Reanalysis_May_2017/Figures/", width = 12, height = 9, units = "in", dpi=500)


# Fig3_DNA concentration --------------------------------------------------




#The Statistic results are from one-way ANOVA as below:

#ANOVA_DNA#

shapiro.test(FT_DC_E1$DNA)
shapiro.test(FT_DC_E2$DNA)
shapiro.test(FT_DC_E3$DNA)
shapiro.test(FT_DC_E4$DNA)


library(agricolae)


av1_DNA <- aov (DNA~Treatment,data=FT_DC_E1)
summary (av1_DNA)
(HSD.test(av1_DNA, "Treatment")) # outer parentheses print result

ST_DNA_E1 <-(HSD.test(av1_DNA, "Treatment"))$groups

ST_DNA_E1$Treatment <-row.names(ST_DNA_E1)

ST_DNA_E1$Pond <- "A"

row.names(ST_DNA_E1)<- NULL



av2_DNA <- aov (DNA~Treatment,data=FT_DC_E2)
summary (av2_DNA)
(HSD.test(av2_DNA, "Treatment")) # outer parentheses print result

ST_DNA_E2 <-(HSD.test(av2_DNA, "Treatment"))$groups

ST_DNA_E2$Treatment <-row.names(ST_DNA_E2)

ST_DNA_E2$Pond <- "B"

row.names(ST_DNA_E2)<- NULL



av3_DNA <- aov (DNA~Treatment,data=FT_DC_E3)
summary (av3_DNA)
(HSD.test(av3_DNA, "Treatment")) # outer parentheses print result

ST_DNA_E3 <-(HSD.test(av3_DNA, "Treatment"))$groups

ST_DNA_E3$Treatment <-row.names(ST_DNA_E3)

ST_DNA_E3$Pond <- "C"

row.names(ST_DNA_E3)<- NULL



av4_DNA <- aov (DNA~Treatment,data=FT_DC_E4)
summary (av4_DNA)
(HSD.test(av4_DNA, "Treatment")) # outer parentheses print result

ST_DNA_E4 <-(HSD.test(av4_DNA, "Treatment"))$groups

ST_DNA_E4$Treatment <-row.names(ST_DNA_E4)

ST_DNA_E4$Pond <- "D"

row.names(ST_DNA_E4)<- NULL

DNA_Statistic <- rbind(ST_DNA_E1,ST_DNA_E2,ST_DNA_E3,ST_DNA_E4)



DC_boxplot <- ggplot(data=FT_DC_original,aes(x=Treatment,y=DNA,fill=Outlier))+
              geom_boxplot(outlier.shape = NA,fill=NA)+#not fill the bars
              scale_fill_manual(values=c("black", "white"))+
              geom_dotplot(binaxis='y', stackdir='center', dotsize=1,width=1, position = "dodge")+
              stat_summary(data=FT_DC, fun.y=mean, mapping = aes(group =Outlier),geom="point", shape=5, size=5)+#The mean calculate exculing outilers
              facet_wrap(~Pond)+
              labs(x="Treatment", y= "DNA concentration (ng/µL)")+theme_bw()+
              theme(text=element_text(size=20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none")



DNA_Statistic$y <- 135


DC_boxplot+geom_text(data=DNA_Statistic,aes(Treatment,y, label=groups, fill=NA),hjust = 0.5, size=7, parse = TRUE)



ggsave("Fig3_DC_boxplot.jpeg",path = "July2016_12s_onestep_filtration/Reanalysis_May_2017/Figures/", width = 12, height = 10, units = "in", dpi=500)





###metaBEAT command###

#metaBEAT.py \
#-Q Querymap.txt -R REFmap.txt --cluster --clust_match 1 --clust_cov 3 --blast --min_ident 1 \
#-m 12S -n 5 -E -v -o 12S-trim_30-merged-nonchimeras-cluster_1c3-blast-min_ident_1.0 &> log


#Avoid false positive approaches#
#1.Low-frequency noise threshold (0.001)
#2.The taxon were found in the negative controls were excluded



#1.Low-frequency noise threshold (0.001)


Filtration <- read.csv(file = "July2016_12s_onestep_filtration/Reanalysis_May_2017/AppendixS2_Filtration_12S_OTU_Nov_2017.csv")


Filtration$Pond <-mapvalues(Filtration$Pond, c("E1","E2","E3","E4"), 
                                              c("A","B","C","D"))

#add the reads assign to Leuciscus_idus into Leuciscus_leuciscus
for (i in 1:nrow(Filtration)){Filtration$Leuciscus_leuciscus[i]=Filtration$Leuciscus_leuciscus[i]+Filtration$Leuciscus_idus[i]; 
Filtration$Leuciscus_idus[i]<-0}
Filtration$Leuciscus_idus <- NULL

#Summary the the total read counts for each sample
Filtration$SUM <- rowSums(Filtration[6:18])

#Reshape the dataset


#Using the 1:5, and 19 variables as ID

Filtration_rs <- melt(Filtration,id=c(1:5,19))

#Rename the variables

Filtration_rs <- rename (Filtration_rs,c("variable"="Species","value"="Reads"))



#Calculate each species reads percentage for each sample

Filtration_rs$Ratio <- Filtration_rs$Reads/Filtration_rs$SUM

Threshold <- 0.001 #based on the positive reads in samples percentage

for (a in 1:nrow(Filtration_rs)) {
  if (!is.na(Filtration_rs$Ratio[a])) {
    if(Filtration_rs$Ratio[a] < Threshold){Filtration_rs$Reads[a] <- 0} 
  }
}



###Samples###
#2.The taxon were found in the negative controls were excluded

Filtration_rs_sample <- Filtration_rs[which(Filtration_rs$Treatment != 'Ctrl' 
                                      & Filtration_rs$Species != 'Astatotilapia_calliptera' 
                                      & Filtration_rs$Species != 'Salmo_trutta' 
                                      & Filtration_rs$Species != 'Alburnus_alburnus'
                                      & Filtration_rs$Species != 'Gobio_gobio'
                                      & Filtration_rs$Species != 'unassigned'),]



Filtration_rs_sample <- Filtration_rs_sample[which(Filtration_rs_sample$Pond !="A" |
                                                     Filtration_rs_sample$Species != "Scardinius_erythrophthalmus"),]

Filtration_rs_sample <- Filtration_rs_sample[which(Filtration_rs_sample$Pond !="A" |
                                                     Filtration_rs_sample$Species != "Leuciscus_leuciscus"),]

Filtration_rs_sample <- Filtration_rs_sample[which(Filtration_rs_sample$Pond !="B" |
                                                     Filtration_rs_sample$Species != "Squalius_cephalus"),]

Filtration_rs_sample <- Filtration_rs_sample[which(Filtration_rs_sample$Pond !="B" |
                                                     Filtration_rs_sample$Species != "Tinca_tinca"),]

Filtration_rs_sample <- Filtration_rs_sample[which(Filtration_rs_sample$Pond !="B" |
                                                     Filtration_rs_sample$Species != "Scardinius_erythrophthalmus"),]


Filtration_rs_sample <- Filtration_rs_sample[which(Filtration_rs_sample$Pond !="C" |
                                                     Filtration_rs_sample$Species != "Barbus_barbus"),]

Filtration_rs_sample <- Filtration_rs_sample[which(Filtration_rs_sample$Pond !="C" |
                                                     Filtration_rs_sample$Species != "Leuciscus_leuciscus"),]

Filtration_rs_sample <- Filtration_rs_sample[which(Filtration_rs_sample$Pond !="C" |
                                                     Filtration_rs_sample$Species != "Scardinius_erythrophthalmus"),]


Filtration_rs_sample <- Filtration_rs_sample[which(Filtration_rs_sample$Pond !="C" |
                                                     Filtration_rs_sample$Species != "Abramis_brama"),]

Filtration_rs_sample <- Filtration_rs_sample[which(Filtration_rs_sample$Pond !="D" |
                                                     Filtration_rs_sample$Species != "Squalius_cephalus"),]

Filtration_rs_sample <- Filtration_rs_sample[which(Filtration_rs_sample$Pond !="D" |
                                                     Filtration_rs_sample$Species != "Rutilus_rutilus"),]



Filtration_rs_sample <- ddply(Filtration_rs_sample, 'SampleID', mutate, Percent_reads = Reads/sum(Reads))

Filtration_rs_sample$Treatment <-mapvalues(Filtration_rs_sample$Treatment, c("T1","T2","T3","T4","T5","P5"), 
                                       c("T1_0.45","T2_0.8","T3_1.2","T4_Sterivex","T5_PF_0.45","TP5_PF"))
#reorder

Filtration_rs_sample$Treatment <- factor(Filtration_rs_sample$Treatment , levels = (c("T1_0.45","T2_0.8","T3_1.2","T4_Sterivex","T5_PF_0.45","TP5_PF")))

Filtration_rs_sample$Species <-mapvalues(Filtration_rs_sample$Species, c("Abramis_brama","Barbus_barbus","Carassius_carassius","Squalius_cephalus","Leuciscus_leuciscus",
                                                                 "Rutilus_rutilus","Scardinius_erythrophthalmus","Tinca_tinca"), 
                                                                  c("BRE","BAR","CAR","CHU","DAC","ROA","RUD","TEN"))


Filtration_rs_sample$Species <- factor(Filtration_rs_sample$Species, levels = (c("BRE","BAR","CAR","CHU","DAC","ROA","RUD","TEN")))







# FigS4_abundance_rep ---------------------------------------------------------------


#identify one more outliers T4-1-3; with two were identified before T3-1-3; T2-2-3

Filtration_rs_sample$Pond <-mapvalues(Filtration_rs_sample$Pond, c("A","B","C","D"), 
                                                             c("E1","E2","E3","E4"))

ddply(Filtration_rs_sample, .(Pond,Treatment,Replicate), summarize,
                       read_sum = round(sum(Reads), 2))


ggplot(Filtration_rs_sample,aes(x=Replicate,y=Percent_reads,fill=Species))+
geom_bar(stat="identity",position="stack",width = 0.8)+facet_wrap(~Treatment+Pond,ncol=4,labeller=labeller(.multi_line = FALSE))+
scale_fill_manual(values=c("firebrick3","green4","blue","orange2","tan4","purple3","gray10","hotpink"))+
scale_y_continuous(labels=percent)+
labs(x="Replicate", y= "Read counts composition")+theme_bw()+
theme(text=element_text(size=15))



ggsave("FigS4_abundance_rep.jpeg",path = "July2016_12s_onestep_filtration/Reanalysis_May_2017/Figures/", width = 10, height = 12, units = "in", dpi=500)






# FigS6_heterogeneity_rep_boxplot_outliers -------------------------------------------------------------------

Filtration_rs_sample$Pond <-mapvalues(Filtration_rs_sample$Pond, c("E1","E2","E3","E4"), 
                                                                    c("A","B","C","D"))


#exclude outliers T2-2-3; T3-1-3; T4-1-3

Filtration_rs_sample_outlier <- Filtration_rs_sample [which(Filtration_rs_sample$SampleID != "T2-2-3"&
                                                            Filtration_rs_sample$SampleID != "T3-1-3" &
                                                            Filtration_rs_sample$SampleID != "T4-1-3"),]

ggplot(Filtration_rs_sample_outlier,aes(x=Treatment,y=Percent_reads))+
  geom_boxplot(aes(colour=Species),outlier.shape = 10, outlier.size = 3)+
  geom_point(size=2, position = position_jitterdodge(jitter.width =0.7,dodge.width = 0.9),aes(colour=Species))+
  facet_wrap(~Pond)+
  scale_colour_manual(values=c("firebrick3","green4","blue","orange2","tan4","purple3","gray10","hotpink"))+
  scale_y_continuous(labels=percent)+
  labs(x="Treatment", y= "Read counts composition")+theme_bw()+
  theme(text=element_text(size=20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.title = element_text(size=15),legend.text = element_text(size=15))



ggsave("FigS6_heterogeneity_rep_boxplot_outliers.jpeg",path = "July2016_12s_onestep_filtration/Reanalysis_May_2017/Figures/", width = 12, height = 10, units = "in", dpi=500)








 
####SpadeR####
 
 
names(Filtration_rs_sample)
 

 
Sample_back <-Filtration_rs_sample[c("SampleID","Pond","Species","Reads")]
 
Sample_back_E1 <- Sample_back[which(Sample_back$Pond == 'A'),]

 
Sample_back_E1_wide <- reshape(Sample_back_E1, idvar= c("Pond","Species"), timevar= "SampleID", direction = "wide")
 
Sample_back_E1_wide <- data.frame(Sample_back_E1_wide[,-1:-2])
 
#excluding identified outliers  T3-1-3; T4-1-3
Sample_back_E1_wide$Reads.T3.1.3 <-NULL

 
Sample_back_E1_wide$Reads.T4.1.3 <-NULL


 
#SimilarityMult(Sample_back_E1_wide, datatype = c("abundance"), q = 2, nboot = 200, goal="relative") 
#the code can not be run, Error in rmultinom(1, ni[k], p[, k]) : NA in probability vector
#write this dataframe then run in the online version of SpadeR; https://chao.shinyapps.io/SpadeR/
#then reformat (txt) the result to csv file


 
write.csv(Sample_back_E1_wide,file = "July2016_12s_onestep_filtration/Reanalysis_May_2017/spadeR/Sample_back_E1_wide.csv",row.names = FALSE)
 
 
 
Sample_back_E2 <- Sample_back[which(Sample_back$Pond == 'B'&
                                      Sample_back$Species != "BRE"),]
 
 
Sample_back_E2_wide <- reshape(Sample_back_E2, idvar= c("Pond","Species"), timevar= "SampleID", direction = "wide")
 
Sample_back_E2_wide <- data.frame(Sample_back_E2_wide[,-1:-2])

#excluding identified outliers T2-2-3
Sample_back_E2_wide$Reads.T2.2.3 <-NULL


 
##SimilarityMult(Sample_back_E2_wide, datatype = c("abundance"), q = 1, nboot = 200, goal="relative")
#the code can not be run, Error in rmultinom(1, ni[k], p[, k]) : NA in probability vector
#write this dataframe then run in the online version of SpadeR; https://chao.shinyapps.io/SpadeR/
#then reformat (txt) the result to csv file
 
write.csv(Sample_back_E2_wide,file = "July2016_12s_onestep_filtration/Reanalysis_May_2017/spadeR/Sample_back_E2_wide.csv",row.names = FALSE)
 
 
 
Sample_back_E3 <- Sample_back[which(Sample_back$Pond == 'C'),]
 
 
Sample_back_E3_wide <- reshape(Sample_back_E3, idvar= c("Pond","Species"), timevar= "SampleID", direction = "wide")
 
Sample_back_E3_wide <- data.frame(Sample_back_E3_wide[,-1:-2])




 
#SimilarityMult(Sample_back_E3_wide, datatype = c("abundance"), q = 1, nboot = 200, goal="relative")
#write this dataframe then run in the online version of SpadeR; https://chao.shinyapps.io/SpadeR/
#then reformat (txt) the result to csv file
 
write.csv(Sample_back_E3_wide,file = "July2016_12s_onestep_filtration/Reanalysis_May_2017/spadeR/Sample_back_E3_wide.csv",row.names = FALSE)
 
 
 
Sample_back_E4 <- Sample_back[which(Sample_back$Pond == 'D'),]
 
 
Sample_back_E4_wide <- reshape(Sample_back_E4, idvar= c("Pond","Species"), timevar= "SampleID", direction = "wide")
 
Sample_back_E4_wide <- data.frame(Sample_back_E4_wide[,-1:-2])
 
 
 
#the code can be run
#SimilarityMult(Sample_back_E4_wide, datatype = c("abundance"), q = 1, nboot = 200, goal="relative")
#write this dataframe then run in the online version of SpadeR; https://chao.shinyapps.io/SpadeR/
#then reformat (txt) the result to csv file
 

write.csv(Sample_back_E4_wide,file = "July2016_12s_onestep_filtration/Reanalysis_May_2017/spadeR/Sample_back_E4_wide.csv",row.names = FALSE)





E1_spade <- read.csv(file = "July2016_12s_onestep_filtration/Reanalysis_May_2017/spadeR/E1_q1.csv",header = TRUE)




for (c in 1:nrow(E1_spade)) {
 if(E1_spade$X[c] <= 5 & E1_spade$Y[c] <=5 ){E1_spade$Group[c] <- 'Within';E1_spade$Treatment[c] <- 'T1'} 
 else if (E1_spade$X[c] > 5 & E1_spade$X[c] <= 10 & E1_spade$Y[c] >5 & E1_spade$Y[c] <=10 ){E1_spade$Group[c] <- 'Within';E1_spade$Treatment[c] <- 'T2'} 
 else if (E1_spade$X[c] > 10 & E1_spade$X[c] <= 14 & E1_spade$Y[c] >10 & E1_spade$Y[c] <=14 ){E1_spade$Group[c] <- 'Within';E1_spade$Treatment[c] <- 'T3'} 
 else if (E1_spade$X[c] > 14 & E1_spade$X[c] <= 18 & E1_spade$Y[c] >14 & E1_spade$Y[c] <=18 ){E1_spade$Group[c] <- 'Within';E1_spade$Treatment[c] <- 'T4'} 
 else if (E1_spade$X[c] > 18 & E1_spade$X[c] <= 23 & E1_spade$Y[c] >18 & E1_spade$Y[c] <=23 ){E1_spade$Group[c] <- 'Within';E1_spade$Treatment[c] <- 'T5'}
 else if (E1_spade$X[c] > 23 & E1_spade$X[c] <= 28 & E1_spade$Y[c] >23 & E1_spade$Y[c] <=28 ){E1_spade$Group[c] <- 'Within';E1_spade$Treatment[c] <- 'TP5'} 
 else if (E1_spade$X[c] <= 5 & E1_spade$Y[c] >5 & E1_spade$Y[c] <=10 ){E1_spade$Group[c] <- 'Among';E1_spade$Treatment[c] <- 'T1T2'}
 else if (E1_spade$X[c] <= 5 & E1_spade$Y[c] >10 & E1_spade$Y[c] <=14 ){E1_spade$Group[c] <- 'Among';E1_spade$Treatment[c] <- 'T1T3'}
 else if (E1_spade$X[c] <= 5 & E1_spade$Y[c] >14 & E1_spade$Y[c] <=18 ){E1_spade$Group[c] <- 'Among';E1_spade$Treatment[c] <- 'T1T4'}
 else if (E1_spade$X[c] <= 5 & E1_spade$Y[c] >18 & E1_spade$Y[c] <=23 ){E1_spade$Group[c] <- 'Among';E1_spade$Treatment[c] <- 'T1T5'}
 else if (E1_spade$X[c] <= 5 & E1_spade$Y[c] >23 & E1_spade$Y[c] <=28 ){E1_spade$Group[c] <- 'Among';E1_spade$Treatment[c] <- 'T1TP5'}
 else if (E1_spade$X[c] > 5 & E1_spade$X[c] <= 10 & E1_spade$Y[c] >10 & E1_spade$Y[c] <=14 ){E1_spade$Group[c] <- 'Among';E1_spade$Treatment[c] <- 'T2T3'} 
 else if (E1_spade$X[c] > 5 & E1_spade$X[c] <= 10 & E1_spade$Y[c] >14 & E1_spade$Y[c] <=18 ){E1_spade$Group[c] <- 'Among';E1_spade$Treatment[c] <- 'T2T4'}
 else if (E1_spade$X[c] > 5 & E1_spade$X[c] <= 10 & E1_spade$Y[c] >18 & E1_spade$Y[c] <=23 ){E1_spade$Group[c] <- 'Among';E1_spade$Treatment[c] <- 'T2T5'}
 else if (E1_spade$X[c] > 5 & E1_spade$X[c] <= 10 & E1_spade$Y[c] >23 & E1_spade$Y[c] <=28 ){E1_spade$Group[c] <- 'Among';E1_spade$Treatment[c] <- 'T2TP5'}
 else if (E1_spade$X[c] > 10 & E1_spade$X[c] <= 14 & E1_spade$Y[c] >14 & E1_spade$Y[c] <=18){E1_spade$Group[c] <- 'Among';E1_spade$Treatment[c] <- 'T3T4'}
 else if (E1_spade$X[c] > 10 & E1_spade$X[c] <= 14 & E1_spade$Y[c] >18 & E1_spade$Y[c] <=23 ){E1_spade$Group[c] <- 'Among';E1_spade$Treatment[c] <- 'T3T5'}
 else if (E1_spade$X[c] > 10 & E1_spade$X[c] <= 14 & E1_spade$Y[c] >23 & E1_spade$Y[c] <=28 ){E1_spade$Group[c] <- 'Among';E1_spade$Treatment[c] <- 'T3TP5'}
 else if (E1_spade$X[c] > 14 & E1_spade$X[c] <= 18 & E1_spade$Y[c] >18 & E1_spade$Y[c] <=23 ){E1_spade$Group[c] <- 'Among';E1_spade$Treatment[c] <- 'T4T5'}
 else if (E1_spade$X[c] > 14 & E1_spade$X[c] <= 18 & E1_spade$Y[c] >23 & E1_spade$Y[c] <=28 ){E1_spade$Group[c] <- 'Among';E1_spade$Treatment[c] <- 'T4TP5'}
 else if (E1_spade$X[c] > 18 & E1_spade$X[c] <= 23 & E1_spade$Y[c] >23 & E1_spade$Y[c] <=28 ){E1_spade$Group[c] <- 'Among';E1_spade$Treatment[c] <- 'T5TP5'}
 else {E1_spade$Group[c] <- 'NA';E1_spade$Treatment[c] <- 'NA'} 
 
}

E1_spade$Pond <- 'A1'

for (d in 1:nrow(E1_spade)) {
 if (E1_spade$Group[d] == 'Among'){E1_spade$Treatment[d] <- 'Among'}
}


E2_spade <- read.csv(file = "July2016_12s_onestep_filtration/Reanalysis_May_2017/spadeR/E2_q1.csv",header = TRUE)




for (c in 1:nrow(E2_spade)) {
  if(E2_spade$X[c] <= 5 & E2_spade$Y[c] <=5 ){E2_spade$Group[c] <- 'Within';E2_spade$Treatment[c] <- 'T1'} 
  else if (E2_spade$X[c] > 5 & E2_spade$X[c] <= 9 & E2_spade$Y[c] >5 & E2_spade$Y[c] <=9 ){E2_spade$Group[c] <- 'Within';E2_spade$Treatment[c] <- 'T2'} 
  else if (E2_spade$X[c] > 9 & E2_spade$X[c] <= 14 & E2_spade$Y[c] >9 & E2_spade$Y[c] <=14 ){E2_spade$Group[c] <- 'Within';E2_spade$Treatment[c] <- 'T3'} 
  else if (E2_spade$X[c] > 14 & E2_spade$X[c] <= 19 & E2_spade$Y[c] >14 & E2_spade$Y[c] <=19 ){E2_spade$Group[c] <- 'Within';E2_spade$Treatment[c] <- 'T4'} 
  else if (E2_spade$X[c] > 19 & E2_spade$X[c] <= 24 & E2_spade$Y[c] >19 & E2_spade$Y[c] <=24 ){E2_spade$Group[c] <- 'Within';E2_spade$Treatment[c] <- 'T5'}
  else if (E2_spade$X[c] > 24 & E2_spade$X[c] <= 29 & E2_spade$Y[c] >24 & E2_spade$Y[c] <=29 ){E2_spade$Group[c] <- 'Within';E2_spade$Treatment[c] <- 'TP5'} 
  else if (E2_spade$X[c] <= 5 & E2_spade$Y[c] >5 & E2_spade$Y[c] <=9 ){E2_spade$Group[c] <- 'Among';E2_spade$Treatment[c] <- 'T1T2'}
  else if (E2_spade$X[c] <= 5 & E2_spade$Y[c] >9 & E2_spade$Y[c] <=14 ){E2_spade$Group[c] <- 'Among';E2_spade$Treatment[c] <- 'T1T3'}
  else if (E2_spade$X[c] <= 5 & E2_spade$Y[c] >14 & E2_spade$Y[c] <=19 ){E2_spade$Group[c] <- 'Among';E2_spade$Treatment[c] <- 'T1T4'}
  else if (E2_spade$X[c] <= 5 & E2_spade$Y[c] >19 & E2_spade$Y[c] <=24 ){E2_spade$Group[c] <- 'Among';E2_spade$Treatment[c] <- 'T1T5'}
  else if (E2_spade$X[c] <= 5 & E2_spade$Y[c] >24 & E2_spade$Y[c] <=29 ){E2_spade$Group[c] <- 'Among';E2_spade$Treatment[c] <- 'T1TP5'}
  else if (E2_spade$X[c] > 5 & E2_spade$X[c] <= 9 & E2_spade$Y[c] >9 & E2_spade$Y[c] <=14 ){E2_spade$Group[c] <- 'Among';E2_spade$Treatment[c] <- 'T2T3'} 
  else if (E2_spade$X[c] > 5 & E2_spade$X[c] <= 9 & E2_spade$Y[c] >14 & E2_spade$Y[c] <=19 ){E2_spade$Group[c] <- 'Among';E2_spade$Treatment[c] <- 'T2T4'}
  else if (E2_spade$X[c] > 5 & E2_spade$X[c] <= 9 & E2_spade$Y[c] >19 & E2_spade$Y[c] <=24 ){E2_spade$Group[c] <- 'Among';E2_spade$Treatment[c] <- 'T2T5'}
  else if (E2_spade$X[c] > 5 & E2_spade$X[c] <= 9 & E2_spade$Y[c] >24 & E2_spade$Y[c] <=29 ){E2_spade$Group[c] <- 'Among';E2_spade$Treatment[c] <- 'T2TP5'}
  else if (E2_spade$X[c] > 9 & E2_spade$X[c] <= 14 & E2_spade$Y[c] >14 & E2_spade$Y[c] <=19){E2_spade$Group[c] <- 'Among';E2_spade$Treatment[c] <- 'T3T4'}
  else if (E2_spade$X[c] > 9 & E2_spade$X[c] <= 14 & E2_spade$Y[c] >19 & E2_spade$Y[c] <=24 ){E2_spade$Group[c] <- 'Among';E2_spade$Treatment[c] <- 'T3T5'}
  else if (E2_spade$X[c] > 9 & E2_spade$X[c] <= 14 & E2_spade$Y[c] >24 & E2_spade$Y[c] <=29 ){E2_spade$Group[c] <- 'Among';E2_spade$Treatment[c] <- 'T3TP5'}
  else if (E2_spade$X[c] > 14 & E2_spade$X[c] <= 19 & E2_spade$Y[c] >19 & E2_spade$Y[c] <=24 ){E2_spade$Group[c] <- 'Among';E2_spade$Treatment[c] <- 'T4T5'}
  else if (E2_spade$X[c] > 14 & E2_spade$X[c] <= 19 & E2_spade$Y[c] >24 & E2_spade$Y[c] <=29 ){E2_spade$Group[c] <- 'Among';E2_spade$Treatment[c] <- 'T4TP5'}
  else if (E2_spade$X[c] > 19 & E2_spade$X[c] <= 24 & E2_spade$Y[c] >24 & E2_spade$Y[c] <=29 ){E2_spade$Group[c] <- 'Among';E2_spade$Treatment[c] <- 'T5TP5'}
  else {E2_spade$Group[c] <- 'NA';E2_spade$Treatment[c] <- 'NA'} 
  
}

E2_spade$Pond <- 'B1'

for (d in 1:nrow(E2_spade)) {
  if (E2_spade$Group[d] == 'Among'){E2_spade$Treatment[d] <- 'Among'}
}



E3_spade <- read.csv(file = "July2016_12s_onestep_filtration/Reanalysis_May_2017/spadeR/E3_q1.csv",header = TRUE)





for (c in 1:nrow(E3_spade)) {
  if(E3_spade$X[c] <= 5 & E3_spade$Y[c] <=5 ){E3_spade$Group[c] <- 'Within';E3_spade$Treatment[c] <- 'T1'} 
  else if (E3_spade$X[c] > 5 & E3_spade$X[c] <= 10 & E3_spade$Y[c] >5 & E3_spade$Y[c] <=10 ){E3_spade$Group[c] <- 'Within';E3_spade$Treatment[c] <- 'T2'} 
  else if (E3_spade$X[c] > 10 & E3_spade$X[c] <= 15 & E3_spade$Y[c] >10 & E3_spade$Y[c] <=15 ){E3_spade$Group[c] <- 'Within';E3_spade$Treatment[c] <- 'T3'} 
  else if (E3_spade$X[c] > 15 & E3_spade$X[c] <= 20 & E3_spade$Y[c] >15 & E3_spade$Y[c] <=20 ){E3_spade$Group[c] <- 'Within';E3_spade$Treatment[c] <- 'T4'} 
  else if (E3_spade$X[c] > 20 & E3_spade$X[c] <= 25 & E3_spade$Y[c] >20 & E3_spade$Y[c] <=25 ){E3_spade$Group[c] <- 'Within';E3_spade$Treatment[c] <- 'T5'}
  else if (E3_spade$X[c] > 25 & E3_spade$X[c] <= 30 & E3_spade$Y[c] >25 & E3_spade$Y[c] <=30 ){E3_spade$Group[c] <- 'Within';E3_spade$Treatment[c] <- 'TP5'} 
  else if (E3_spade$X[c] <= 5 & E3_spade$Y[c] >5 & E3_spade$Y[c] <=10 ){E3_spade$Group[c] <- 'Among';E3_spade$Treatment[c] <- 'T1T2'}
  else if (E3_spade$X[c] <= 5 & E3_spade$Y[c] >10 & E3_spade$Y[c] <=15 ){E3_spade$Group[c] <- 'Among';E3_spade$Treatment[c] <- 'T1T3'}
  else if (E3_spade$X[c] <= 5 & E3_spade$Y[c] >15 & E3_spade$Y[c] <=20 ){E3_spade$Group[c] <- 'Among';E3_spade$Treatment[c] <- 'T1T4'}
  else if (E3_spade$X[c] <= 5 & E3_spade$Y[c] >20 & E3_spade$Y[c] <=25 ){E3_spade$Group[c] <- 'Among';E3_spade$Treatment[c] <- 'T1T5'}
  else if (E3_spade$X[c] <= 5 & E3_spade$Y[c] >25 & E3_spade$Y[c] <=30 ){E3_spade$Group[c] <- 'Among';E3_spade$Treatment[c] <- 'T1TP5'}
  else if (E3_spade$X[c] > 5 & E3_spade$X[c] <= 10 & E3_spade$Y[c] >10 & E3_spade$Y[c] <=15 ){E3_spade$Group[c] <- 'Among';E3_spade$Treatment[c] <- 'T2T3'} 
  else if (E3_spade$X[c] > 5 & E3_spade$X[c] <= 10 & E3_spade$Y[c] >15 & E3_spade$Y[c] <=20 ){E3_spade$Group[c] <- 'Among';E3_spade$Treatment[c] <- 'T2T4'}
  else if (E3_spade$X[c] > 5 & E3_spade$X[c] <= 10 & E3_spade$Y[c] >20 & E3_spade$Y[c] <=25 ){E3_spade$Group[c] <- 'Among';E3_spade$Treatment[c] <- 'T2T5'}
  else if (E3_spade$X[c] > 5 & E3_spade$X[c] <= 10 & E3_spade$Y[c] >25 & E3_spade$Y[c] <=30 ){E3_spade$Group[c] <- 'Among';E3_spade$Treatment[c] <- 'T2TP5'}
  else if (E3_spade$X[c] > 10 & E3_spade$X[c] <= 15 & E3_spade$Y[c] >15 & E3_spade$Y[c] <=20){E3_spade$Group[c] <- 'Among';E3_spade$Treatment[c] <- 'T3T4'}
  else if (E3_spade$X[c] > 10 & E3_spade$X[c] <= 15 & E3_spade$Y[c] >20 & E3_spade$Y[c] <=25 ){E3_spade$Group[c] <- 'Among';E3_spade$Treatment[c] <- 'T3T5'}
  else if (E3_spade$X[c] > 10 & E3_spade$X[c] <= 15 & E3_spade$Y[c] >25 & E3_spade$Y[c] <=30 ){E3_spade$Group[c] <- 'Among';E3_spade$Treatment[c] <- 'T3TP5'}
  else if (E3_spade$X[c] > 15 & E3_spade$X[c] <= 20 & E3_spade$Y[c] >20 & E3_spade$Y[c] <=25 ){E3_spade$Group[c] <- 'Among';E3_spade$Treatment[c] <- 'T4T5'}
  else if (E3_spade$X[c] > 15 & E3_spade$X[c] <= 20 & E3_spade$Y[c] >25 & E3_spade$Y[c] <=30 ){E3_spade$Group[c] <- 'Among';E3_spade$Treatment[c] <- 'T4TP5'}
  else if (E3_spade$X[c] > 20 & E3_spade$X[c] <= 25 & E3_spade$Y[c] >25 & E3_spade$Y[c] <=30 ){E3_spade$Group[c] <- 'Among';E3_spade$Treatment[c] <- 'T5TP5'}
  else {E3_spade$Group[c] <- 'NA';E3_spade$Treatment[c] <- 'NA'} 
  
}

E3_spade$Pond <- 'C1'

for (d in 1:nrow(E3_spade)) {
  if (E3_spade$Group[d] == 'Among'){E3_spade$Treatment[d] <- 'Among'}
}




E4_spade <- read.csv(file = "July2016_12s_onestep_filtration/Reanalysis_May_2017/spadeR/E4_q1.csv",header = TRUE)





for (c in 1:nrow(E4_spade)) {
  if(E4_spade$X[c] <= 5 & E4_spade$Y[c] <=5 ){E4_spade$Group[c] <- 'Within';E4_spade$Treatment[c] <- 'T1'} 
  else if (E4_spade$X[c] > 5 & E4_spade$X[c] <= 10 & E4_spade$Y[c] >5 & E4_spade$Y[c] <=10 ){E4_spade$Group[c] <- 'Within';E4_spade$Treatment[c] <- 'T2'} 
  else if (E4_spade$X[c] > 10 & E4_spade$X[c] <= 15 & E4_spade$Y[c] >10 & E4_spade$Y[c] <=15 ){E4_spade$Group[c] <- 'Within';E4_spade$Treatment[c] <- 'T3'} 
  else if (E4_spade$X[c] > 15 & E4_spade$X[c] <= 20 & E4_spade$Y[c] >15 & E4_spade$Y[c] <=20 ){E4_spade$Group[c] <- 'Within';E4_spade$Treatment[c] <- 'T4'} 
  else if (E4_spade$X[c] > 20 & E4_spade$X[c] <= 25 & E4_spade$Y[c] >20 & E4_spade$Y[c] <=25 ){E4_spade$Group[c] <- 'Within';E4_spade$Treatment[c] <- 'T5'}
  else if (E4_spade$X[c] > 25 & E4_spade$X[c] <= 30 & E4_spade$Y[c] >25 & E4_spade$Y[c] <=30 ){E4_spade$Group[c] <- 'Within';E4_spade$Treatment[c] <- 'TP5'} 
  else if (E4_spade$X[c] <= 5 & E4_spade$Y[c] >5 & E4_spade$Y[c] <=10 ){E4_spade$Group[c] <- 'Among';E4_spade$Treatment[c] <- 'T1T2'}
  else if (E4_spade$X[c] <= 5 & E4_spade$Y[c] >10 & E4_spade$Y[c] <=15 ){E4_spade$Group[c] <- 'Among';E4_spade$Treatment[c] <- 'T1T3'}
  else if (E4_spade$X[c] <= 5 & E4_spade$Y[c] >15 & E4_spade$Y[c] <=20 ){E4_spade$Group[c] <- 'Among';E4_spade$Treatment[c] <- 'T1T4'}
  else if (E4_spade$X[c] <= 5 & E4_spade$Y[c] >20 & E4_spade$Y[c] <=25 ){E4_spade$Group[c] <- 'Among';E4_spade$Treatment[c] <- 'T1T5'}
  else if (E4_spade$X[c] <= 5 & E4_spade$Y[c] >25 & E4_spade$Y[c] <=30 ){E4_spade$Group[c] <- 'Among';E4_spade$Treatment[c] <- 'T1TP5'}
  else if (E4_spade$X[c] > 5 & E4_spade$X[c] <= 10 & E4_spade$Y[c] >10 & E4_spade$Y[c] <=15 ){E4_spade$Group[c] <- 'Among';E4_spade$Treatment[c] <- 'T2T3'} 
  else if (E4_spade$X[c] > 5 & E4_spade$X[c] <= 10 & E4_spade$Y[c] >15 & E4_spade$Y[c] <=20 ){E4_spade$Group[c] <- 'Among';E4_spade$Treatment[c] <- 'T2T4'}
  else if (E4_spade$X[c] > 5 & E4_spade$X[c] <= 10 & E4_spade$Y[c] >20 & E4_spade$Y[c] <=25 ){E4_spade$Group[c] <- 'Among';E4_spade$Treatment[c] <- 'T2T5'}
  else if (E4_spade$X[c] > 5 & E4_spade$X[c] <= 10 & E4_spade$Y[c] >25 & E4_spade$Y[c] <=30 ){E4_spade$Group[c] <- 'Among';E4_spade$Treatment[c] <- 'T2TP5'}
  else if (E4_spade$X[c] > 10 & E4_spade$X[c] <= 15 & E4_spade$Y[c] >15 & E4_spade$Y[c] <=20){E4_spade$Group[c] <- 'Among';E4_spade$Treatment[c] <- 'T3T4'}
  else if (E4_spade$X[c] > 10 & E4_spade$X[c] <= 15 & E4_spade$Y[c] >20 & E4_spade$Y[c] <=25 ){E4_spade$Group[c] <- 'Among';E4_spade$Treatment[c] <- 'T3T5'}
  else if (E4_spade$X[c] > 10 & E4_spade$X[c] <= 15 & E4_spade$Y[c] >25 & E4_spade$Y[c] <=30 ){E4_spade$Group[c] <- 'Among';E4_spade$Treatment[c] <- 'T3TP5'}
  else if (E4_spade$X[c] > 15 & E4_spade$X[c] <= 20 & E4_spade$Y[c] >20 & E4_spade$Y[c] <=25 ){E4_spade$Group[c] <- 'Among';E4_spade$Treatment[c] <- 'T4T5'}
  else if (E4_spade$X[c] > 15 & E4_spade$X[c] <= 20 & E4_spade$Y[c] >25 & E4_spade$Y[c] <=30 ){E4_spade$Group[c] <- 'Among';E4_spade$Treatment[c] <- 'T4TP5'}
  else if (E4_spade$X[c] > 20 & E4_spade$X[c] <= 25 & E4_spade$Y[c] >25 & E4_spade$Y[c] <=30 ){E4_spade$Group[c] <- 'Among';E4_spade$Treatment[c] <- 'T5TP5'}
  else {E4_spade$Group[c] <- 'NA';E4_spade$Treatment[c] <- 'NA'} 
  
}

E4_spade$Pond <- 'D1'


for (d in 1:nrow(E4_spade)) {
  if (E4_spade$Group[d] == 'Among'){E4_spade$Treatment[d] <- 'Among'}
}



spadeR_all <- rbind(E1_spade,E2_spade,E3_spade,E4_spade)




spadeR_all$Treatment <-mapvalues(spadeR_all$Treatment, c("T1","T2","T3","T4","T5","TP5"), 
                                   c("T1_0.45","T2_0.8","T3_1.2","T4_Sterivex","T5_PF_0.45","TP5_PF"))





### Dunn test


E1_spade$Treatment <- mapvalues(E1_spade$Treatment, c("Among"), 
                                c("Umong"))

DT_spadeR_E1 <- dunnTest(Estimate ~ Treatment, data = E1_spade, method="bh") 


DT_spadeR_E1 <- DT_spadeR_E1$res

DT_spadeR_E1 


ST_spadeR_E1 <-cldList(comparison = DT_spadeR_E1$Comparison, p.value = DT_spadeR_E1$P.unadj, threshold= 0.05)

ST_spadeR_E1$Pond <- "A1"




### Dunn test

E2_spade$Treatment <- mapvalues(E2_spade$Treatment, c("Among"), 
                                c("Umong"))


DT_spadeR_E2 <- dunnTest(Estimate ~ Treatment, data = E2_spade, method="bh") 


DT_spadeR_E2 <- DT_spadeR_E2$res

DT_spadeR_E2 


ST_spadeR_E2 <-cldList(comparison = DT_spadeR_E2$Comparison, p.value = DT_spadeR_E2$P.unadj, threshold= 0.05)

ST_spadeR_E2$Pond <- "B1"





E3_spade$Treatment <- mapvalues(E3_spade$Treatment, c("Among"), 
                                c("Umong"))


DT_spadeR_E3 <- dunnTest(Estimate ~ Treatment, data = E3_spade, method="bh") 


DT_spadeR_E3 <- DT_spadeR_E3$res

DT_spadeR_E3 


ST_spadeR_E3 <-cldList(comparison = DT_spadeR_E3$Comparison, p.value = DT_spadeR_E3$P.unadj, threshold= 0.05)

ST_spadeR_E3$Pond <- "C1"





E4_spade$Treatment <- mapvalues(E4_spade$Treatment, c("Among"), 
                                c("Umong"))


DT_spadeR_E4 <- dunnTest(Estimate ~ Treatment, data = E4_spade, method="bh") 


DT_spadeR_E4 <- DT_spadeR_E4$res

DT_spadeR_E4 


ST_spadeR_E4 <-cldList(comparison = DT_spadeR_E4$Comparison, p.value = DT_spadeR_E4$P.unadj, threshold= 0.05)

ST_spadeR_E4$Pond <- "D1"





spadeR_Statistic <- rbind(ST_spadeR_E1,ST_spadeR_E2,ST_spadeR_E3,ST_spadeR_E4)

spadeR_Statistic <- rename (spadeR_Statistic,c("Group"="Treatment"))

spadeR_Statistic$Treatment <-mapvalues(spadeR_Statistic$Treatment, c("T1","T2","T3", "T4", "T5", "TP5","Umong"), 
                                   c("T1_0.45","T2_0.8","T3_1.2","T4_Sterivex","T5_PF_0.45","TP5_PF","Among"))





#summarise Similarity Horn index ------------------

Horn_average <- ddply(spadeR_all, .(Pond,Treatment), summarize,
                Estimate_mean = round(mean(Estimate), 2),
                Estimate_sd = round(sd(Estimate), 2))

ddply(spadeR_all, .(Treatment), summarize,
      Estimate_mean = round(mean(Estimate), 2),
      Estimate_sd = round(sd(Estimate), 2))

#write.csv(Horn_average,file = "July2016_12s_onestep_filtration/Reanalysis_May_2017/Horn_average.csv",row.names = FALSE)

# Fig5_Horn-------------------------------------------------------------------



spadeR_all_boxplot <- ggplot(spadeR_all,aes(x=Treatment,y=Estimate))+
                      geom_boxplot(outlier.shape = NA)+
                      geom_jitter(size=2,width =0.2)+
                      facet_wrap(~Pond)+
                      scale_y_continuous(breaks=c(0.00,0.25,0.50,0.75,1.00))+
                      scale_x_discrete(limits=c("T1_0.45","T2_0.8","T3_1.2","T4_Sterivex","T5_PF_0.45","TP5_PF","Among"))+
                      labs(x="Treatment",y= "Horn similarity index")+theme_bw()+
                      theme(text=element_text(size=20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


spadeR_Statistic$y <- 1.1



spadeR_all_boxplot2 <- spadeR_all_boxplot+geom_text(data=spadeR_Statistic,aes(Treatment,y, label=Letter),hjust = 0.5, size=7, parse = TRUE)




#ggsave("Fig5_Similarity_Pond.jpeg",path = "July2016_12s_onestep_filtration/Reanalysis_May_2017/Figures/", width = 12, height = 9, units = "in", dpi=500)





# Fig5_NMDS----------------------------------------------------------------



library(vegan)


#NMDS_E1-----------------------------------------------------------------------------
NMDS_E1 <- Filtration_rs_sample[which(Filtration_rs_sample$Pond == 'A'),]



#excluding identified outliers  T3-1-3; T4-1-3

NMDS_E1 <- NMDS_E1[which(NMDS_E1$SampleID != 'T3-1-3' &
                         NMDS_E1$SampleID != 'T4-1-3'),]

#excluding some variables TP, Pond, Replicate, SUM, Percent_reads, Ratio
NMDS_E1 <- NMDS_E1 [c(-2,-4:-6,-9:-10)]



NMDS_E1_wide <- reshape(NMDS_E1, idvar= c("SampleID","Treatment"), timevar= "Species", direction = "wide")

NMDS_E1_wide <- rename (NMDS_E1_wide,c("Reads.BRE"="BRE","Reads.BAR"="BAR","Reads.CAR"="CAR",
                                     "Reads.ROA"="ROA","Reads.CHU"="CHU","Reads.TEN"="TEN"))


NMDS_E1_wide.rownames <- data.frame(NMDS_E1_wide[,-1], row.names=NMDS_E1_wide[,1])


NMDS_E1_wide.rownames <- na.omit(NMDS_E1_wide.rownames) 



NMDS_E1_wide.active <- NMDS_E1_wide.rownames[, 2:7]

NMDS_E1_wide.active <- na.omit(NMDS_E1_wide.active) 




E1_NMDS<-metaMDS(NMDS_E1_wide.active, distance="bray", k=2, trymax=1000)

#check stressplot
stressplot(E1_NMDS)



#Using ggplot for the NMDS plot https://chrischizinski.github.io/rstats/vegan-ggplot2/

#creat the treatment inforamation as group information

E1_treat<-c(rep("T1_0.45",5),rep("T2_0.8",5),rep("T3_1.2",4),rep("T4_Sterivex",4),rep("T5_PF_0.45",5),rep("TP5_PF",5))

#Using the scores function from vegan to extract the site scores and convert to a data.frame
E1_data.scores <- as.data.frame(scores(E1_NMDS))

# create a column of replicate names, from the rownames of data.scores
E1_data.scores$Rep <- rownames(E1_data.scores)

#  add the Treatment variable created earlier
E1_data.scores$Treatment <- E1_treat

E1_data.scores$Pond <- "A2"

#look at the data
head(E1_data.scores)  

E1_data.scores

#Using the scores function from vegan to extract the species scores and convert to a data.frame
E1_species.scores <- as.data.frame(scores(E1_NMDS, "species"))

# create a column of species, from the rownames of species.scores
E1_species.scores$species <- rownames(E1_species.scores)

E1_species.scores$Pond <- "A2"

#look at the data
head(E1_species.scores)


#Using the chull function from vegan to extract the hull values and convert to a data.frame

E1_grp.T1 <- E1_data.scores[E1_data.scores$Treatment == "T1_0.45", ][chull(E1_data.scores[E1_data.scores$Treatment == 
                                                                   "T1_0.45", c("NMDS1", "NMDS2")]), ]

E1_grp.T2 <- E1_data.scores[E1_data.scores$Treatment == "T2_0.8", ][chull(E1_data.scores[E1_data.scores$Treatment == 
                                                                                      "T2_0.8", c("NMDS1", "NMDS2")]), ]

E1_grp.T3 <- E1_data.scores[E1_data.scores$Treatment == "T3_1.2", ][chull(E1_data.scores[E1_data.scores$Treatment == 
                                                                                      "T3_1.2", c("NMDS1", "NMDS2")]), ]
E1_grp.T4 <- E1_data.scores[E1_data.scores$Treatment == "T4_Sterivex", ][chull(E1_data.scores[E1_data.scores$Treatment == 
                                                                                      "T4_Sterivex", c("NMDS1", "NMDS2")]), ]
E1_grp.T5 <- E1_data.scores[E1_data.scores$Treatment == "T5_PF_0.45", ][chull(E1_data.scores[E1_data.scores$Treatment == 
                                                                                      "T5_PF_0.45", c("NMDS1", "NMDS2")]), ]
E1_grp.TP5 <- E1_data.scores[E1_data.scores$Treatment == "TP5_PF", ][chull(E1_data.scores[E1_data.scores$Treatment == 
                                                                                      "TP5_PF", c("NMDS1", "NMDS2")]), ]

#combine Treatment data

E1_hull.data <- rbind(E1_grp.T1, E1_grp.T2, E1_grp.T3, E1_grp.T4, E1_grp.T5, E1_grp.TP5)  
E1_hull.data$Pond <- "A2"



#NMDS_E2-----------------------------------------------------------------------------


NMDS_E2 <- Filtration_rs_sample[which(Filtration_rs_sample$Pond == 'B'),]



#excluding identified outliers  T2-2-3; T5-2-1

NMDS_E2 <- NMDS_E2[which(NMDS_E2$SampleID != 'T2-2-3'),]



#excluding some variables TP, Pond, Replicate, SUM, Percent_reads, Ratio
NMDS_E2 <- NMDS_E2 [c(-2,-4:-6,-9:-10)]



NMDS_E2_wide <- reshape(NMDS_E2, idvar= c("SampleID","Treatment"), timevar= "Species", direction = "wide")

NMDS_E2_wide <- rename (NMDS_E2_wide,c("Reads.BRE"="BRE","Reads.BAR"="BAR","Reads.CAR"="CAR",
                                       "Reads.ROA"="ROA","Reads.DAC"="DAC"))





NMDS_E2_wide.rownames <- data.frame(NMDS_E2_wide[,-1], row.names=NMDS_E2_wide[,1])


NMDS_E2_wide.rownames <- na.omit(NMDS_E2_wide.rownames) 



NMDS_E2_wide.active <- NMDS_E2_wide.rownames[, 3:6]

NMDS_E2_wide.active <- na.omit(NMDS_E2_wide.active) 



E2_NMDS=metaMDS(NMDS_E2_wide.active, distance="bray", k=2, trymax=1000)

#check stressplot
stressplot(E2_NMDS)



#Using ggplot for the NMDS plot https://chrischizinski.github.io/rstats/vegan-ggplot2/

#creat the treatment inforamation as group information

E2_treat=c(rep("T1_0.45",5),rep("T2_0.8",4),rep("T3_1.2",5),rep("T4_Sterivex",5),rep("T5_PF_0.45",5),rep("TP5_PF",5))

#Using the scores function from vegan to extract the site scores and convert to a data.frame
E2_data.scores <- as.data.frame(scores(E2_NMDS))

# create a column of replicate names, from the rownames of data.scores
E2_data.scores$Rep <- rownames(E2_data.scores)

#  add the Treatment variable created earlier
E2_data.scores$Treatment <- E2_treat  

E2_data.scores$Pond <- "B2"

#look at the data
head(E2_data.scores)  

E2_data.scores

#Using the scores function from vegan to extract the species scores and convert to a data.frame
E2_species.scores <- as.data.frame(scores(E2_NMDS, "species"))

# create a column of species, from the rownames of species.scores
E2_species.scores$species <- rownames(E2_species.scores)

#look at the data

E2_species.scores$Pond <- "B2"
head(E2_species.scores)


#Using the chull function from vegan to extract the hull values and convert to a data.frame

E2_grp.T1 <- E2_data.scores[E2_data.scores$Treatment == "T1_0.45", ][chull(E2_data.scores[E2_data.scores$Treatment == 
                                                                                            "T1_0.45", c("NMDS1", "NMDS2")]), ]

E2_grp.T2 <- E2_data.scores[E2_data.scores$Treatment == "T2_0.8", ][chull(E2_data.scores[E2_data.scores$Treatment == 
                                                                                           "T2_0.8", c("NMDS1", "NMDS2")]), ]

E2_grp.T3 <- E2_data.scores[E2_data.scores$Treatment == "T3_1.2", ][chull(E2_data.scores[E2_data.scores$Treatment == 
                                                                                           "T3_1.2", c("NMDS1", "NMDS2")]), ]
E2_grp.T4 <- E2_data.scores[E2_data.scores$Treatment == "T4_Sterivex", ][chull(E2_data.scores[E2_data.scores$Treatment == 
                                                                                                "T4_Sterivex", c("NMDS1", "NMDS2")]), ]
E2_grp.T5 <- E2_data.scores[E2_data.scores$Treatment == "T5_PF_0.45", ][chull(E2_data.scores[E2_data.scores$Treatment == 
                                                                                               "T5_PF_0.45", c("NMDS1", "NMDS2")]), ]
E2_grp.TP5 <- E2_data.scores[E2_data.scores$Treatment == "TP5_PF", ][chull(E2_data.scores[E2_data.scores$Treatment == 
                                                                                            "TP5_PF", c("NMDS1", "NMDS2")]), ]

#combine Treatment data

E2_hull.data <- rbind(E2_grp.T1, E2_grp.T2, E2_grp.T3, E2_grp.T4, E2_grp.T5, E2_grp.TP5)  
E2_hull.data$Pond <- "B2"


#NMDS_E3-------------------------------------------------------------------------


NMDS_E3 <- Filtration_rs_sample[which(Filtration_rs_sample$Pond == 'C'),]



#excluding identified outliers  T4-3-2
NMDS_E3 <- NMDS_E3[which(NMDS_E3$SampleID != 'T4-3-2'),]




#excluding some variables TP, Pond, Replicate, SUM, Percent_reads, Ratio
NMDS_E3 <- NMDS_E3 [c(-2,-4:-6,-9:-10)]



NMDS_E3_wide <- reshape(NMDS_E3, idvar= c("SampleID","Treatment"), timevar= "Species", direction = "wide")

NMDS_E3_wide <- rename (NMDS_E3_wide,c("Reads.CHU"="CHU","Reads.TEN"="TEN","Reads.CAR"="CAR",
                                       "Reads.ROA"="ROA"))



NMDS_E3_wide.rownames <- data.frame(NMDS_E3_wide[,-1], row.names=NMDS_E3_wide[,1])


NMDS_E3_wide.rownames <- na.omit(NMDS_E3_wide.rownames) 



NMDS_E3_wide.active <- NMDS_E3_wide.rownames[, 2:5]

NMDS_E3_wide.active <- na.omit(NMDS_E3_wide.active) 



E3_NMDS<-metaMDS(NMDS_E3_wide.active, distance="bray", k=2, trymax=1000)

#check stressplot
stressplot(E3_NMDS)



#Using ggplot for the NMDS plot https://chrischizinski.github.io/rstats/vegan-ggplot2/

#creat the treatment inforamation as group information

E3_treat=c(rep("T1_0.45",5),rep("T2_0.8",5),rep("T3_1.2",5),rep("T4_Sterivex",4),rep("T5_PF_0.45",5),rep("TP5_PF",5))

#Using the scores function from vegan to extract the site scores and convert to a data.frame
E3_data.scores <- as.data.frame(scores(E3_NMDS))

# create a column of replicate names, from the rownames of data.scores
E3_data.scores$Rep <- rownames(E3_data.scores)

#  add the Treatment variable created earlier
E3_data.scores$Treatment <- E3_treat  

E3_data.scores$Pond <- "C2"

#look at the data
head(E3_data.scores)  

E3_data.scores

#Using the scores function from vegan to extract the species scores and convert to a data.frame
E3_species.scores <- as.data.frame(scores(E3_NMDS, "species"))

# create a column of species, from the rownames of species.scores
E3_species.scores$species <- rownames(E3_species.scores)

#look at the data
E3_species.scores$Pond <- "C2"
head(E3_species.scores)


#Using the chull function from vegan to extract the hull values and convert to a data.frame

E3_grp.T1 <- E3_data.scores[E3_data.scores$Treatment == "T1_0.45", ][chull(E3_data.scores[E3_data.scores$Treatment == 
                                                                                            "T1_0.45", c("NMDS1", "NMDS2")]), ]

E3_grp.T2 <- E3_data.scores[E3_data.scores$Treatment == "T2_0.8", ][chull(E3_data.scores[E3_data.scores$Treatment == 
                                                                                           "T2_0.8", c("NMDS1", "NMDS2")]), ]

E3_grp.T3 <- E3_data.scores[E3_data.scores$Treatment == "T3_1.2", ][chull(E3_data.scores[E3_data.scores$Treatment == 
                                                                                           "T3_1.2", c("NMDS1", "NMDS2")]), ]
E3_grp.T4 <- E3_data.scores[E3_data.scores$Treatment == "T4_Sterivex", ][chull(E3_data.scores[E3_data.scores$Treatment == 
                                                                                                "T4_Sterivex", c("NMDS1", "NMDS2")]), ]
E3_grp.T5 <- E3_data.scores[E3_data.scores$Treatment == "T5_PF_0.45", ][chull(E3_data.scores[E3_data.scores$Treatment == 
                                                                                               "T5_PF_0.45", c("NMDS1", "NMDS2")]), ]
E3_grp.TP5 <- E3_data.scores[E3_data.scores$Treatment == "TP5_PF", ][chull(E3_data.scores[E3_data.scores$Treatment == 
                                                                                            "TP5_PF", c("NMDS1", "NMDS2")]), ]

#combine Treatment data

E3_hull.data <- rbind(E3_grp.T1, E3_grp.T2, E3_grp.T3, E3_grp.T4, E3_grp.T5, E3_grp.TP5)  
E3_hull.data$Pond <- "C2"



#NMDS_E4-------------------------------------------------------------------------


NMDS_E4 <- Filtration_rs_sample[which(Filtration_rs_sample$Pond == 'D'),]




#excluding some variables TP, Pond, Replicate, SUM, Percent_reads, Ratio
NMDS_E4 <- NMDS_E4 [c(-2,-4:-6,-9:-10)]



NMDS_E4_wide <- reshape(NMDS_E4, idvar= c("SampleID","Treatment"), timevar= "Species", direction = "wide")

NMDS_E4_wide <- rename (NMDS_E4_wide,c("Reads.BRE"="BRE","Reads.BAR"="BAR","Reads.RUD"="RUD","Reads.TEN"="TEN","Reads.CAR"="CAR",
                                       "Reads.DAC"="DAC"))



NMDS_E4_wide.rownames <- data.frame(NMDS_E4_wide[,-1], row.names=NMDS_E4_wide[,1])


NMDS_E4_wide.rownames <- na.omit(NMDS_E4_wide.rownames) 



NMDS_E4_wide.active <- NMDS_E4_wide.rownames[, 2:7]

NMDS_E4_wide.active <- na.omit(NMDS_E4_wide.active) 



E4_NMDS=metaMDS(NMDS_E4_wide.active, distance="bray", k=2, trymax=1000)

#check stressplot
stressplot(E4_NMDS)



#Using ggplot for the NMDS plot https://chrischizinski.github.io/rstats/vegan-ggplot2/

#creat the treatment inforamation as group information

E4_treat=c(rep("T1_0.45",5),rep("T2_0.8",5),rep("T3_1.2",5),rep("T4_Sterivex",5),rep("T5_PF_0.45",5),rep("TP5_PF",5))

#Using the scores function from vegan to extract the site scores and convert to a data.frame
E4_data.scores <- as.data.frame(scores(E4_NMDS))

# create a column of replicate names, from the rownames of data.scores
E4_data.scores$Rep <- rownames(E4_data.scores)

#  add the Treatment variable created earlier
E4_data.scores$Treatment <- E4_treat  
E4_data.scores$Pond <- "D2"

#look at the data
head(E4_data.scores)  

E4_data.scores

#Using the scores function from vegan to extract the species scores and convert to a data.frame
E4_species.scores <- as.data.frame(scores(E4_NMDS, "species"))

# create a column of species, from the rownames of species.scores
E4_species.scores$species <- rownames(E4_species.scores)

#look at the data
E4_species.scores$Pond <- "D2"
head(E4_species.scores)


#Using the chull function from vegan to extract the hull values and convert to a data.frame

E4_grp.T1 <- E4_data.scores[E4_data.scores$Treatment == "T1_0.45", ][chull(E4_data.scores[E4_data.scores$Treatment == 
                                                                                            "T1_0.45", c("NMDS1", "NMDS2")]), ]

E4_grp.T2 <- E4_data.scores[E4_data.scores$Treatment == "T2_0.8", ][chull(E4_data.scores[E4_data.scores$Treatment == 
                                                                                           "T2_0.8", c("NMDS1", "NMDS2")]), ]

E4_grp.T3 <- E4_data.scores[E4_data.scores$Treatment == "T3_1.2", ][chull(E4_data.scores[E4_data.scores$Treatment == 
                                                                                           "T3_1.2", c("NMDS1", "NMDS2")]), ]
E4_grp.T4 <- E4_data.scores[E4_data.scores$Treatment == "T4_Sterivex", ][chull(E4_data.scores[E4_data.scores$Treatment == 
                                                                                                "T4_Sterivex", c("NMDS1", "NMDS2")]), ]
E4_grp.T5 <- E4_data.scores[E4_data.scores$Treatment == "T5_PF_0.45", ][chull(E4_data.scores[E4_data.scores$Treatment == 
                                                                                               "T5_PF_0.45", c("NMDS1", "NMDS2")]), ]
E4_grp.TP5 <- E4_data.scores[E4_data.scores$Treatment == "TP5_PF", ][chull(E4_data.scores[E4_data.scores$Treatment == 
                                                                                            "TP5_PF", c("NMDS1", "NMDS2")]), ]

#combine Treatment data

E4_hull.data <- rbind(E4_grp.T1, E4_grp.T2, E4_grp.T3, E4_grp.T4, E4_grp.T5, E4_grp.TP5)  
E4_hull.data$Pond <- "D2"




NMDS_data.scores_pond <- rbind(E1_data.scores,E2_data.scores,E3_data.scores,E4_data.scores)
NMDS_species.scores_pond <-  rbind(E1_species.scores,E2_species.scores,E3_species.scores,E4_species.scores)

NMDS_hull.data_pond <- rbind(E1_hull.data,E2_hull.data,E3_hull.data,E4_hull.data)


Pond_NMDS_0.5 <- ggplot() + 
                stat_ellipse(data=NMDS_hull.data_pond,aes(x=NMDS1,y=NMDS2,colour=Treatment,group=Treatment),type = "norm",linetype = "dashed", level = 0.5, size=0.8) + # add the convex hulls
                geom_text(data=NMDS_species.scores_pond,aes(x=NMDS1,y=NMDS2,label=species),colour="magenta2",size=4) +  # add the species labels
                geom_point(data=NMDS_data.scores_pond,aes(x=NMDS1,y=NMDS2,shape=Treatment,colour=Treatment),size=4) + # add the point markers
                facet_wrap(~Pond)+
                scale_colour_manual(values=c("firebrick1","green","blue","cyan2","gold","gray10")) +
                theme_bw()+theme(text=element_text(size=20),legend.position="bottom",legend.title = element_text(size=15),legend.text = element_text(size=15))




Horn_NMDS <- grid.arrange(spadeR_all_boxplot2, Pond_NMDS_0.5,ncol = 1, nrow = 2,
            heights = c(3.5,3.5))



ggsave("Fig5_Horn_NMDS_0.5.jpeg",Horn_NMDS,path = "July2016_12s_onestep_filtration/Reanalysis_May_2017/Figures/", width = 12, height = 18, units = "in", dpi=500)



####Analysis of similarities (ANOSIM)####

#E1.anosim -------------------------------------------------------------------


#Global test
NMDS_E1.dist <- vegdist(NMDS_E1_wide.active,method="bray")


E1.env <- droplevels(NMDS_E1_wide$Treatment)

E1.env <- as.data.frame(E1.env)

E1.env<- rename (E1.env,c("E1.env"="Treatment"))



E1.anosim<- anosim(NMDS_E1.dist, E1.env$Treatment)
summary(E1.anosim)
plot(E1.anosim)



NMDS_E1 <- Filtration_rs_sample[which(Filtration_rs_sample$Pond == 'A'),]


#excluding identified outliers  T3-1-3; T4-1-3

NMDS_E1 <- NMDS_E1[which(NMDS_E1$SampleID != 'T3-1-3' &
                           NMDS_E1$SampleID != 'T4-1-3'),]

#Subset the dataset then pairwise compared

NMDS_E1_1 <- NMDS_E1[which(NMDS_E1$Treatment == "T1_0.45"|
                             NMDS_E1$Treatment == "T2_0.8"),]

NMDS_E1_1 <- NMDS_E1[which(NMDS_E1$Treatment == "T1_0.45"|
                             NMDS_E1$Treatment == "T3_1.2"),]

NMDS_E1_1 <- NMDS_E1[which(NMDS_E1$Treatment == "T2_0.8"|
                             NMDS_E1$Treatment == "T3_1.2"),]

NMDS_E1_1 <- NMDS_E1[which(NMDS_E1$Treatment == "T1_0.45"|
                             NMDS_E1$Treatment == "T4_Sterivex"),]

NMDS_E1_1 <- NMDS_E1[which(NMDS_E1$Treatment == "T1_0.45"|
                             NMDS_E1$Treatment == "T5_PF_0.45"),]

NMDS_E1_1 <- NMDS_E1[which(NMDS_E1$Treatment == "TP5_PF"|
                             NMDS_E1$Treatment == "T5_PF_0.45"),]

NMDS_E1_1 <- NMDS_E1[which(NMDS_E1$Treatment == "T1_0.45"|
                             NMDS_E1$Treatment == "TP5_PF"),]

####Run the above codes one by one, then run the anosim as below

NMDS_E1_1 <- NMDS_E1_1 [c(-2,-4:-6,-9:-10)]

NMDS_E1_wide_1 <- reshape(NMDS_E1_1, idvar= c("SampleID","Treatment"), timevar= "Species", direction = "wide")

NMDS_E1_wide_1 <- rename (NMDS_E1_wide_1,c("Reads.BRE"="BRE","Reads.BAR"="BAR","Reads.CAR"="CAR",
                                       "Reads.ROA"="ROA","Reads.CHU"="CHU","Reads.TEN"="TEN"))


NMDS_E1_wide.rownames_1 <- data.frame(NMDS_E1_wide_1[,-1], row.names=NMDS_E1_wide_1[,1])


NMDS_E1_wide.rownames_1 <- na.omit(NMDS_E1_wide.rownames_1) 



NMDS_E1_wide.active_1 <- NMDS_E1_wide.rownames_1[, 2:7]

NMDS_E1_wide.active_1 <- na.omit(NMDS_E1_wide.active_1) 




  
NMDS_E1.dist_1 <- vegdist(NMDS_E1_wide.active_1,method="bray")


E1.env_1 <- droplevels(NMDS_E1_wide_1$Treatment)

E1.env_1 <- as.data.frame(E1.env_1)

E1.env_1<- rename (E1.env_1,c("E1.env_1"="Treatment"))



E1.anosim_1<- anosim(NMDS_E1.dist_1, E1.env_1$Treatment)
summary(E1.anosim_1)
plot(E1.anosim_1)


#E2.anosim -------------------------------------------------------------------

#Global test

NMDS_E2.dist <- vegdist(NMDS_E2_wide.active,method="bray")


E2.env <- NMDS_E2_wide$Treatment

E2.env <- as.data.frame(E2.env)

E2.env<- rename (E2.env,c("E2.env"="Treatment"))


E2.anosim <- anosim(NMDS_E2.dist, E2.env$Treatment)
summary(E2.anosim)
plot(E2.anosim)




NMDS_E2 <- Filtration_rs_sample[which(Filtration_rs_sample$Pond == 'B'),]



#excluding identified outliers  T2-2-3; T5-2-1

NMDS_E2 <- NMDS_E2[which(NMDS_E2$SampleID != 'T2-2-3'),]


#Subset the dataset then pairwise compared

NMDS_E2_1 <- NMDS_E2[which(NMDS_E2$Treatment == "T1_0.45"|
                             NMDS_E2$Treatment == "T2_0.8"),]

NMDS_E2_1 <- NMDS_E2[which(NMDS_E2$Treatment == "T1_0.45"|
                             NMDS_E2$Treatment == "T3_1.2"),]

NMDS_E2_1 <- NMDS_E2[which(NMDS_E2$Treatment == "T2_0.8"|
                             NMDS_E2$Treatment == "T3_1.2"),]

NMDS_E2_1 <- NMDS_E2[which(NMDS_E2$Treatment == "T1_0.45"|
                             NMDS_E2$Treatment == "T4_Sterivex"),]

NMDS_E2_1 <- NMDS_E2[which(NMDS_E2$Treatment == "T1_0.45"|
                             NMDS_E2$Treatment == "T5_PF_0.45"),]

NMDS_E2_1 <- NMDS_E2[which(NMDS_E2$Treatment == "TP5_PF"|
                             NMDS_E2$Treatment == "T5_PF_0.45"),]

NMDS_E2_1 <- NMDS_E2[which(NMDS_E2$Treatment == "T1_0.45"|
                             NMDS_E2$Treatment == "TP5_PF"),]

####Run the above codes one by one, then run the anosim as below



NMDS_E2_1 <- NMDS_E2_1 [c(-2,-4:-6,-9:-10)]

NMDS_E2_wide_1 <- reshape(NMDS_E2_1, idvar= c("SampleID","Treatment"), timevar= "Species", direction = "wide")

NMDS_E2_wide_1 <- rename (NMDS_E2_wide_1,c("Reads.BRE"="BRE","Reads.BAR"="BAR","Reads.CAR"="CAR",
                                           "Reads.ROA"="ROA","Reads.DAC"="DAC"))


NMDS_E2_wide.rownames_1 <- data.frame(NMDS_E2_wide_1[,-1], row.names=NMDS_E2_wide_1[,1])


NMDS_E2_wide.rownames_1 <- na.omit(NMDS_E2_wide.rownames_1) 



NMDS_E2_wide.active_1 <- NMDS_E2_wide.rownames_1[, 3:6]

NMDS_E2_wide.active_1 <- na.omit(NMDS_E2_wide.active_1) 



NMDS_E2.dist_1 <- vegdist(NMDS_E2_wide.active_1,method="bray")


E2.env_1 <- droplevels(NMDS_E2_wide_1$Treatment)

E2.env_1 <- as.data.frame(E2.env_1)

E2.env_1<- rename (E2.env_1,c("E2.env_1"="Treatment"))



E2.anosim_1<- anosim(NMDS_E2.dist_1, E2.env_1$Treatment)
summary(E2.anosim_1)
plot(E2.anosim_1)




#E3.anosim ---------------------------------------------------

#Global test
NMDS_E3.dist <- vegdist(NMDS_E3_wide.active,method="bray")


E3.env <- NMDS_E3_wide$Treatment

E3.env <- as.data.frame(E3.env)

E3.env<- rename (E3.env,c("E3.env"="Treatment"))


E3.anosim <- anosim(NMDS_E3.dist, E3.env$Treatment)
summary(E3.anosim)
plot(E3.anosim)



NMDS_E3 <- Filtration_rs_sample[which(Filtration_rs_sample$Pond == 'C'),]

#excluding identified outliers  T4-3-2
NMDS_E3 <- NMDS_E3[which(NMDS_E3$SampleID != 'T4-3-2'),]


#Subset the dataset then pairwise compared

NMDS_E3_1 <- NMDS_E3[which(NMDS_E3$Treatment == "T1_0.45"|
                             NMDS_E3$Treatment == "T2_0.8"),]

NMDS_E3_1 <- NMDS_E3[which(NMDS_E3$Treatment == "T1_0.45"|
                             NMDS_E3$Treatment == "T3_1.2"),]

NMDS_E3_1 <- NMDS_E3[which(NMDS_E3$Treatment == "T2_0.8"|
                             NMDS_E3$Treatment == "T3_1.2"),]

NMDS_E3_1 <- NMDS_E3[which(NMDS_E3$Treatment == "T1_0.45"|
                             NMDS_E3$Treatment == "T4_Sterivex"),]

NMDS_E3_1 <- NMDS_E3[which(NMDS_E3$Treatment == "T1_0.45"|
                             NMDS_E3$Treatment == "T5_PF_0.45"),]

NMDS_E3_1 <- NMDS_E3[which(NMDS_E3$Treatment == "TP5_PF"|
                             NMDS_E3$Treatment == "T5_PF_0.45"),]

NMDS_E3_1 <- NMDS_E3[which(NMDS_E3$Treatment == "T1_0.45"|
                             NMDS_E3$Treatment == "TP5_PF"),]

####Run the above codes one by one, then run the anosim as below




NMDS_E3_1 <- NMDS_E3_1 [c(-2,-4:-6,-9:-10)]

NMDS_E3_wide_1 <- reshape(NMDS_E3_1, idvar= c("SampleID","Treatment"), timevar= "Species", direction = "wide")

NMDS_E3_wide_1 <- rename (NMDS_E3_wide_1,c("Reads.CHU"="CHU","Reads.TEN"="TEN","Reads.CAR"="CAR",
                                           "Reads.ROA"="ROA"))



NMDS_E3_wide.rownames_1 <- data.frame(NMDS_E3_wide_1[,-1], row.names=NMDS_E3_wide_1[,1])


NMDS_E3_wide.rownames_1 <- na.omit(NMDS_E3_wide.rownames_1) 



NMDS_E3_wide.active_1 <- NMDS_E3_wide.rownames_1[, 2:5]

NMDS_E3_wide.active_1 <- na.omit(NMDS_E3_wide.active_1) 





NMDS_E3.dist_1 <- vegdist(NMDS_E3_wide.active_1,method="bray")


E3.env_1 <- droplevels(NMDS_E3_wide_1$Treatment)

E3.env_1 <- as.data.frame(E3.env_1)

E3.env_1<- rename (E3.env_1,c("E3.env_1"="Treatment"))



E3.anosim_1<- anosim(NMDS_E3.dist_1, E3.env_1$Treatment)
summary(E3.anosim_1)
plot(E3.anosim_1)




#E4.anosim ----------------------------------------------------

#Global test

NMDS_E4.dist <- vegdist(NMDS_E4_wide.active,method="bray")


E4.env <- NMDS_E4_wide$Treatment

E4.env <- as.data.frame(E4.env)

E4.env<- rename (E4.env,c("E4.env"="Treatment"))


E4.anosim <- anosim(NMDS_E4.dist, E4.env$Treatment)
summary(E4.anosim)
plot(E4.anosim)





NMDS_E4 <- Filtration_rs_sample[which(Filtration_rs_sample$Pond == 'D'),]



#Subset the dataset then pairwise compared

NMDS_E4_1 <- NMDS_E4[which(NMDS_E4$Treatment == "T1_0.45"|
                             NMDS_E4$Treatment == "T2_0.8"),]

NMDS_E4_1 <- NMDS_E4[which(NMDS_E4$Treatment == "T1_0.45"|
                             NMDS_E4$Treatment == "T3_1.2"),]

NMDS_E4_1 <- NMDS_E4[which(NMDS_E4$Treatment == "T2_0.8"|
                             NMDS_E4$Treatment == "T3_1.2"),]

NMDS_E4_1 <- NMDS_E4[which(NMDS_E4$Treatment == "T1_0.45"|
                             NMDS_E4$Treatment == "T4_Sterivex"),]

NMDS_E4_1 <- NMDS_E4[which(NMDS_E4$Treatment == "T1_0.45"|
                             NMDS_E4$Treatment == "T5_PF_0.45"),]

NMDS_E4_1 <- NMDS_E4[which(NMDS_E4$Treatment == "TP5_PF"|
                             NMDS_E4$Treatment == "T5_PF_0.45"),]

NMDS_E4_1 <- NMDS_E4[which(NMDS_E4$Treatment == "T1_0.45"|
                             NMDS_E4$Treatment == "TP5_PF"),]

####Run the above codes one by one, then run the anosim as below




NMDS_E4_1 <- NMDS_E4_1 [c(-2,-4:-6,-9:-10)]

NMDS_E4_wide_1 <- reshape(NMDS_E4_1, idvar= c("SampleID","Treatment"), timevar= "Species", direction = "wide")

NMDS_E4_wide_1 <- rename (NMDS_E4_wide_1,c("Reads.BRE"="BRE","Reads.BAR"="BAR","Reads.RUD"="RUD","Reads.TEN"="TEN","Reads.CAR"="CAR",
                                           "Reads.DAC"="DAC"))



NMDS_E4_wide.rownames_1 <- data.frame(NMDS_E4_wide_1[,-1], row.names=NMDS_E4_wide_1[,1])


NMDS_E4_wide.rownames_1 <- na.omit(NMDS_E4_wide.rownames_1) 



NMDS_E4_wide.active_1 <- NMDS_E4_wide.rownames_1[, 2:7]

NMDS_E4_wide.active_1 <- na.omit(NMDS_E4_wide.active_1) 





NMDS_E4.dist_1 <- vegdist(NMDS_E4_wide.active_1,method="bray")


E4.env_1 <- droplevels(NMDS_E4_wide_1$Treatment)

E4.env_1 <- as.data.frame(E4.env_1)

E4.env_1<- rename (E4.env_1,c("E4.env_1"="Treatment"))



E4.anosim_1<- anosim(NMDS_E4.dist_1, E4.env_1$Treatment)
summary(E4.anosim_1)
plot(E4.anosim_1)











Filtration_sample_RN <- Filtration_rs_sample[which(Filtration_rs_sample$Reads >0),]

#exclude outliers T3-1-3; T4-1-3; T2-2-3

Filtration_sample_RN <- Filtration_sample_RN [which(Filtration_sample_RN$SampleID != "T3-1-3" &
                                                    Filtration_sample_RN$SampleID != "T4-1-3" &
                                                    Filtration_sample_RN$SampleID != "T2-2-3" ),]

Filtration_sample_RN2 = data.frame(TP=factor(Filtration_sample_RN$TP),Treatment=factor(Filtration_sample_RN$Treatment),
                            Pond=factor(Filtration_sample_RN$Pond),Species=factor(Filtration_sample_RN$Species))



Filtration_rs_dat =dcast(Filtration_sample_RN2, TP+Treatment+Pond~Species, fun.aggregate = length)


Filtration_rs_dat.melt =melt(Filtration_rs_dat,id=c("TP","Treatment","Pond"))

Filtration_rs_dat.melt <- rename (Filtration_rs_dat.melt,c("variable"="Species","value"="RN"))



Filtration_rs_sum <- aggregate(Reads~TP+Treatment+Pond+Species, data = Filtration_sample_RN, FUN=sum)

Filtration_rs_sum_dat <-dcast(Filtration_rs_sum, TP+Treatment+Pond~Species)

Filtration_rs_sum_dat.melt = melt(Filtration_rs_sum_dat,id= c("TP","Treatment","Pond"))

Filtration_rs_sum_dat.melt <- rename(Filtration_rs_sum_dat.melt,c("variable"="Species","value"="Total_Reads"))

Merge_Filtration <- merge(Filtration_rs_dat.melt, Filtration_rs_sum_dat.melt, by=c("TP","Treatment","Pond","Species"))



#calculate the average reads



for (b in 1:nrow(Merge_Filtration)) {
    if(Merge_Filtration$TP[b] == "T3E1"){Merge_Filtration$Reads[b] <- Merge_Filtration$Total_Reads[b]/4} 
    else if(Merge_Filtration$TP[b] == "T4E1"){Merge_Filtration$Reads[b] <- Merge_Filtration$Total_Reads[b]/4} 
    else if(Merge_Filtration$TP[b] == "T2E2"){Merge_Filtration$Reads[b] <- Merge_Filtration$Total_Reads[b]/4}
    else {Merge_Filtration$Reads[b] <- Merge_Filtration$Total_Reads[b]/5}
}




Merge_Filtration <- Merge_Filtration[which(Merge_Filtration$Pond !="A" |
                                             Merge_Filtration$Species != "RUD"),]

Merge_Filtration <- Merge_Filtration[which(Merge_Filtration$Pond !="A" |
                                             Merge_Filtration$Species != "DAC"),]

Merge_Filtration <- Merge_Filtration[which(Merge_Filtration$Pond !="B" |
                                             Merge_Filtration$Species != "CHU"),]

Merge_Filtration <- Merge_Filtration[which(Merge_Filtration$Pond !="B" |
                                             Merge_Filtration$Species != "TEN"),]

Merge_Filtration <- Merge_Filtration[which(Merge_Filtration$Pond !="B" |
                                             Merge_Filtration$Species != "RUD"),]


Merge_Filtration <- Merge_Filtration[which(Merge_Filtration$Pond !="C" |
                                             Merge_Filtration$Species != "BAR"),]

Merge_Filtration <- Merge_Filtration[which(Merge_Filtration$Pond !="C" |
                                             Merge_Filtration$Species != "DAC"),]

Merge_Filtration <- Merge_Filtration[which(Merge_Filtration$Pond !="C" |
                                             Merge_Filtration$Species != "RUD"),]


Merge_Filtration <- Merge_Filtration[which(Merge_Filtration$Pond !="C" |
                                             Merge_Filtration$Species != "BRE"),]

Merge_Filtration <- Merge_Filtration[which(Merge_Filtration$Pond !="D" |
                                             Merge_Filtration$Species != "CHU"),]

Merge_Filtration <- Merge_Filtration[which(Merge_Filtration$Pond !="D" |
                                             Merge_Filtration$Species != "ROA"),]

Merge_Filtration [is.na(Merge_Filtration)] <- 0






#change the read counts to percentage
Merge_Filtration <- ddply(Merge_Filtration, 'TP', mutate, Percent_reads = Reads/sum(Reads))




Merge_Filtration <-Merge_Filtration[order(Merge_Filtration$Species),]




ggplot(Merge_Filtration,aes(x=Treatment,y=Percent_reads,fill=Species))+
            geom_bar(stat="identity",position="stack",width = 0.8)+facet_wrap(~Pond)+
            scale_fill_manual(values=c("firebrick3","green4","blue","orange2","tan4","purple3","gray10","hotpink"))+
            scale_y_continuous(labels=percent)+
            scale_x_discrete(limits=c("T1_0.45","T2_0.8","T3_1.2","T4_Sterivex","T5_PF_0.45","TP5_PF"))+
            labs(x="Treatment", y= "Averaged read counts percentage")+theme_bw()+
            theme(text=element_text(size=20),axis.text.x = element_text(angle = -15,vjust = 0.7))





####Pond fish summary####


Pond_summary <- read.csv(file = "July2016_12s_onestep_filtration/Reanalysis_May_2017/Calverton stock information summary 2015.csv")

Pond_summary$Species <-mapvalues(Pond_summary$Species, c("Abramis_brama","Barbus_barbus","Carassius_carassius","Squalius_cephalus","Leuciscus_leuciscus",
                                                                 "Rutilus_rutilus","Scardinius_erythrophthalmus","Tinca_tinca"), 
                                                             c("BRE","BAR","CAR","CHU","DAC","ROA","RUD","TEN"))

Pond_summary$Pond <-mapvalues(Pond_summary$Pond, c("E1","E2","E3","E4"), 
                                                  c("A","B","C","D"))

# FigS1 AND S2 ------------------------------------------------------------



ggplot(Pond_summary, aes(x= Date, y= Number,colour=Species,group=Species))+
  geom_point(stat="identity",size=3.5) + geom_line(size=0.8,linetype = 3)+
  facet_wrap(~Pond,ncol=2)+
  scale_colour_manual(values=c("firebrick3","green4","blue","orange2","tan4","purple3","gray10","hotpink"))+
  scale_x_discrete(limits=c("Jun-15","Jul-15","Aug-15",'Sep-15','Oct-15','Nov-15','Dec-15',"Jan-16"))+
  labs(x="Month", y="Fish abundance")+theme_bw()+
  theme(text=element_text(size=20),axis.text.x = element_text(angle = -25,vjust = 0.7),legend.title = element_text(size=15),legend.text = element_text(size=15))

ggsave("FigS1_number.jpeg",path = "July2016_12s_onestep_filtration/Reanalysis_May_2017/Figures/", width = 12, height = 9, units = "in", dpi=500)




ggplot(Pond_summary, aes(x= Date, y= Biomass,colour=Species,group=Species))+
  geom_point(stat="identity",size=3.5) + geom_line(size=0.8,linetype = 3)+
  facet_wrap(~Pond,ncol=2)+
  scale_colour_manual(values=c("firebrick3","green4","blue","orange2","tan4","purple3","gray10","hotpink"))+
  scale_x_discrete(limits=c("Jun-15","Jul-15","Aug-15",'Sep-15','Oct-15','Nov-15','Dec-15',"Jan-16"))+
  labs(x="Month", y="Fish biomass (kg)")+theme_bw()+
  theme(text=element_text(size=20),axis.text.x = element_text(angle = -25,vjust = 0.7),legend.title = element_text(size=15),legend.text = element_text(size=15))

ggsave("FigS2_biomass.jpeg",path = "July2016_12s_onestep_filtration/Reanalysis_May_2017/Figures/", width = 12, height = 9, units = "in", dpi=500)


Pond_summary_Aug <- Pond_summary[which(Pond_summary$Date == 'Aug-15'),]

Pond_summary_Aug <- ddply(Pond_summary_Aug, 'Pond', mutate, Percent_number = Number/sum(Number))

Pond_summary_Aug <- ddply(Pond_summary_Aug, 'Pond', mutate, Percent_biomass = Biomass/sum(Biomass))




# Fig4_relative abundance_bio_num --------------------------------------------------------------------


Merge_Filtration_sub <- Merge_Filtration[c(2:4,8)]

Pond_summary_Aug_sub <- Pond_summary_Aug[c(1:2,8)]

Pond_summary_Aug_sub$Treatment <- 'Bio'


Pond_summary_Aug_sub <- rename (Pond_summary_Aug_sub,c("Percent_biomass"="Percent_reads"))


Pond_summary_Aug_sub_num <- Pond_summary_Aug[c(1:2,7)]

Pond_summary_Aug_sub_num$Treatment <- 'Abu'

Pond_summary_Aug_sub_num <- rename (Pond_summary_Aug_sub_num,c("Percent_number"="Percent_reads"))

Merge_Filtration_sub_bio_num<-rbind(Merge_Filtration_sub,Pond_summary_Aug_sub,Pond_summary_Aug_sub_num)


ggplot(Merge_Filtration_sub_bio_num,aes(x=Treatment,y=Percent_reads,fill=Species))+
  geom_bar(stat="identity",position="stack",width = 0.8)+facet_wrap(~Pond)+
  scale_fill_manual(values=c("firebrick3","green4","blue","orange2","tan4","purple3","gray10","hotpink"))+
  scale_y_continuous(labels=percent)+
  scale_x_discrete(limits=c("T1_0.45","T2_0.8","T3_1.2","T4_Sterivex","T5_PF_0.45","TP5_PF","Bio","Abu"))+
  labs(x="Treatment", y= "Species composition")+theme_bw()+
  theme(text=element_text(size=20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.title = element_text(size=15),legend.text = element_text(size=15))



ggsave("Fig4_relative abundance_bio_num.jpeg",path = "July2016_12s_onestep_filtration/Reanalysis_May_2017/Figures/", width = 12, height = 10, units = "in", dpi=500)




#####Relationships#####

#Combine the sequencing reads and fish summary infromation

corfun<-function(x, y) {
  corr=(cor.test(x, y,
                 alternative="two.sided", method="pearson"))}


Filtration_summary <-merge(Merge_Filtration, Pond_summary_Aug, by=c("Pond","Species"))






Filtration_summary$Pond <-mapvalues(Filtration_summary$Pond, c("A","B","C","D"), 
                                                    c("E1","E2","E3","E4"))




for (c in 1:nrow(Filtration_summary)) {
  if(Filtration_summary$Reads[c] == 0){Filtration_summary$Reads[c] <- 1} 
}



Filtration_summary$lgReads <- log10(Filtration_summary$Reads)
Filtration_summary$lgNum <- log10(Filtration_summary$Number)
Filtration_summary$lgBio <- log10(Filtration_summary$Biomass)

names(Filtration_summary)


Filtration_summary_sub <- Filtration_summary[c("Pond","Species","Treatment","lgReads","lgNum","lgBio")]

Filtration_summary_sub_E1 <- Filtration_summary_sub[which(Filtration_summary_sub$Pond == "E1"),]
Filtration_summary_sub_E2 <- Filtration_summary_sub[which(Filtration_summary_sub$Pond == "E2"),]
Filtration_summary_sub_E3 <- Filtration_summary_sub[which(Filtration_summary_sub$Pond == "E3"),]
Filtration_summary_sub_E4 <- Filtration_summary_sub[which(Filtration_summary_sub$Pond == "E4"),]





#Analysis of Covariance
#http://rcompanion.org/rcompanion/e_04.html

#gReads ~ lgNum + Treatment + lgNum:Treatment

options(contrasts = c("contr.treatment", "contr.poly"))

model.1 = lm (lgReads ~ lgNum + Treatment + lgNum:Treatment,
              data = Filtration_summary_sub_E1)

model.2 = lm (lgReads ~ lgNum + Treatment + lgNum:Treatment,
              data = Filtration_summary_sub_E2)

model.3 = lm (lgReads ~ lgNum + Treatment + lgNum:Treatment,
              data = Filtration_summary_sub_E3)

model.4 = lm (lgReads ~ lgNum + Treatment + lgNum:Treatment,
              data = Filtration_summary_sub_E4)


library(car)

Anova(model.1, type="II")


Anova(model.2, type="II")

Anova(model.3, type="II")

Anova(model.4, type="II")

#lgReads ~ lgBio + Treatment + lgBio:Treatment

model.1 = lm (lgReads ~ lgBio + Treatment + lgBio:Treatment,
              data = Filtration_summary_sub_E1)

model.2 = lm (lgReads ~ lgBio + Treatment + lgBio:Treatment,
              data = Filtration_summary_sub_E2)

model.3 = lm (lgReads ~ lgBio + Treatment + lgBio:Treatment,
              data = Filtration_summary_sub_E3)

model.4 = lm (lgReads ~ lgBio + Treatment + lgBio:Treatment,
              data = Filtration_summary_sub_E4)
library(car)

Anova(model.1, type="II")


Anova(model.2, type="II")

Anova(model.3, type="II")

Anova(model.4, type="II")

# Fig6 -----------------------------------------------------------------
#lgreads and lgnum
Filtration3 <- ggplot(Filtration_summary, aes(x= Number, y=Reads, shape=Pond))+
                geom_point(stat="identity",size=3)+
                facet_wrap(~Pond+Treatment,ncol=6,labeller=labeller(.multi_line = FALSE))+
                scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                              labels = trans_format("log10", math_format(10^.x)))+
                scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                              labels = trans_format("log10", math_format(10^.x)))+
                geom_smooth(formula= y ~ x,method="lm",se=FALSE, size=0.5, na.rm=TRUE,colour='black')+
                labs(x="Fish abundance", y="Averaged read counts")+theme_bw()+
                geom_text_repel(aes(label=Species),size=3)+
                theme(text=element_text(size=15))





Filtration3_log_num_reads<- ddply(Filtration_summary, .(Pond,Treatment), summarise, s=corfun(lgNum,lgReads)$statistic,
                                  cor.est=corfun(lgNum,lgReads)$estimate,
                                  pval=corfun(lgNum,lgReads)$p.value,
                                  alt=corfun(lgNum,lgReads)$alternative,
                                  df=corfun(lgNum,lgReads)$parameter)
#Table2 summary

Corr_fishnum <- ddply(Filtration_summary, .(Treatment), summarise, s=corfun(lgNum,lgReads)$statistic,
                cor.est=corfun(lgNum,lgReads)$estimate,
                pval=corfun(lgNum,lgReads)$p.value,
                alt=corfun(lgNum,lgReads)$alternative,
                df=corfun(lgNum,lgReads)$parameter)

#write.csv(Corr_fishnum,file = "July2016_12s_onestep_filtration/Reanalysis_May_2017/Corr_fishnum.csv",row.names = FALSE)


Filtration3_len <- length(levels(Filtration_summary$Pond))*length(levels(Filtration_summary$Treatment))

Filtration3_vars <- data.frame(expand.grid(levels(Filtration_summary$Treatment),levels(Filtration_summary$Pond)))

colnames(Filtration3_vars) <- c("Treatment","Pond")

Filtration3_dat1 <- data.frame(x = rep(1e3, Filtration3_len), y = rep(1.1e1, Filtration3_len), 
                               Filtration3_vars, labs=round(Filtration3_log_num_reads$cor.est,2))

Filtration3_dat2 <- data.frame(x = rep(1e3, Filtration3_len), y = rep(0.4e1, Filtration3_len), 
                               Filtration3_vars, labs=round(Filtration3_log_num_reads$pval,2))
Filtration3_dat3 <- data.frame(x = rep(1e2, Filtration3_len), y = rep(5e3, Filtration3_len), 
                               Filtration3_vars, labs=paste0(LETTERS[1:24]))



for (c in 1:nrow(Filtration3_dat2)) {
  if (!is.na(Filtration3_dat2$labs[c])) {
    if(Filtration3_dat2$labs[c]< 1e-02){Filtration3_dat2$pvalue[c] <- "P< 0.01"} 
    else if (Filtration3_dat2$labs[c]< 5e-02){Filtration3_dat2$pvalue[c] <- "P< 0.05"} 
    else {Filtration3_dat2$pvalue[c] <- paste("P==", Filtration3_dat2$labs[c]) }
  }
}


Filtration3 + geom_text(data=Filtration3_dat1, aes(x, y, label=paste("Cor==", labs)),hjust = 0, size=4, parse = TRUE)+
  geom_text(data=Filtration3_dat2, aes(x, y, label=pvalue),hjust = 0, size=4, parse = TRUE)+
  geom_text(data=Filtration3_dat3, aes(x, y, label=labs),hjust = 0, size=4, parse = TRUE)


ggsave("Fig6_Lg_Num_reads.jpeg",path ="July2016_12s_onestep_filtration/Reanalysis_May_2017/Figures/", width = 12, height = 9, units = "in", dpi=500)



# FigS7 -----------------------------------------------------------------
#lgreads and lgnum

Filtration4 <- ggplot(Filtration_summary, aes(x= Biomass, y=Reads, shape=Pond))+
  geom_point(stat="identity",size=3)+
  facet_wrap(~Pond+Treatment,ncol=6,labeller=labeller(.multi_line = FALSE))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  geom_smooth(formula= y ~ x,method="lm",se=FALSE, size=0.5, na.rm=TRUE,colour='black')+
  labs(x="Fish biomass (kg)", y="Averaged read counts")+theme_bw()+
  geom_text_repel(aes(label=Species),size=3)+
  theme(text=element_text(size=15))






Filtration4_log_bio_reads<- ddply(Filtration_summary, .(Pond,Treatment), summarise, s=corfun(lgBio,lgReads)$statistic,
                                  cor.est=corfun(lgBio,lgReads)$estimate,
                                  pval=corfun(lgBio,lgReads)$p.value,
                                  alt=corfun(lgBio,lgReads)$alternative,
                                  df=corfun(lgBio,lgReads)$parameter)

#Table2 summary

Corr_fishbio <-ddply(Filtration_summary, .(Treatment), summarise, s=corfun(lgBio,lgReads)$statistic,
                cor.est=corfun(lgBio,lgReads)$estimate,
                pval=corfun(lgBio,lgReads)$p.value,
                alt=corfun(lgBio,lgReads)$alternative,
                df=corfun(lgBio,lgReads)$parameter)


#write.csv(Corr_fishbio,file = "July2016_12s_onestep_filtration/Reanalysis_May_2017/Corr_fishbio.csv",row.names = FALSE)


Filtration4_len <- length(levels(Filtration_summary$Pond))*length(levels(Filtration_summary$Treatment))

Filtration4_vars <- data.frame(expand.grid(levels(Filtration_summary$Treatment),levels(Filtration_summary$Pond)))

colnames(Filtration4_vars) <- c("Treatment","Pond")

Filtration4_dat1 <- data.frame(x = rep(0.2e0, Filtration4_len), y = rep(0.3e2, Filtration4_len), 
                               Filtration4_vars, labs=round(Filtration4_log_bio_reads$cor.est,2))

Filtration4_dat2 <- data.frame(x = rep(0.2e0, Filtration4_len), y = rep(1e1, Filtration4_len), 
                               Filtration4_vars, labs=round(Filtration4_log_bio_reads$pval,2))
Filtration4_dat3 <- data.frame(x = rep(0.2e0, Filtration4_len), y = rep(3e3, Filtration4_len), 
                               Filtration4_vars, labs=paste0(LETTERS[1:24]))

for (c in 1:nrow(Filtration4_dat2)) {
  if (!is.na(Filtration4_dat2$labs[c])) {
    if(Filtration4_dat2$labs[c]< 1e-02){Filtration4_dat2$pvalue[c] <- "P< 0.01"} 
    else if (Filtration4_dat2$labs[c]< 5e-02){Filtration4_dat2$pvalue[c] <- "P< 0.05"} 
    else {Filtration4_dat2$pvalue[c] <- paste("P==", Filtration4_dat2$labs[c]) }
  }
}


Filtration4 + geom_text(data=Filtration4_dat1, aes(x, y, label=paste("Cor==", labs)),hjust = 0, size=4, parse = TRUE)+
  geom_text(data=Filtration4_dat2, aes(x, y, label=pvalue),hjust = 0, size=4, parse = TRUE)+
  geom_text(data=Filtration4_dat3, aes(x, y, label=labs),hjust = 0, size=4, parse = TRUE)

ggsave("FigS7_Lg_Bio_reads.jpeg",path ="July2016_12s_onestep_filtration/Reanalysis_May_2017/Figures/", width = 12, height = 9, units = "in", dpi=500)



#####two way ANOVA####

#http://www.biostathandbook.com/twowayanova.html
#https://rcompanion.org/rcompanion/d_08.html

#FT_DC

shapiro.test(transformTukey(FT_DC$Time))
qqnorm(transformTukey(FT_DC$Time))

shapiro.test(sqrt(FT_DC$DNA))
qqnorm(transformTukey(FT_DC$DNA))


if(!require(FSA)){install.packages("FSA")}
if(!require(multcompView)){install.packages("multcompView")}
if(!require(lsmeans)){install.packages("lsmeans")}
if(!require(grid)){install.packages("grid")}
if(!require(nlme)){install.packages("nlme")}
if(!require(lme4)){install.packages("lme4")}
if(!require(lmerTest)){install.packages("lmerTest")} 
if(!require(Rmisc)){install.packages("Rmisc")} 

library(Rmisc)##need this one

####FT####

model_FT = lm(transformTukey(Time) ~ Treatment + Pond + Pond:Treatment,
                 data=FT_DC)

Anova(model_FT, type="II") 

FT_DC_P0 <- FT_DC[which(FT_DC$Treatment !='TP5_PF'),]


model_FT_0 = lm(transformTukey(Time) ~ Treatment + Pond + Pond:Treatment,
              data=FT_DC_P0)

Anova(model_FT_0, type="II") 



FT_DC_P1 <- FT_DC[which(FT_DC$Treatment =='T1_0.45'|
                        FT_DC$Treatment =='T2_0.8'|
                        FT_DC$Treatment =='T3_1.2'),]

model_FT_P1 = lm(transformTukey(Time) ~ Treatment + Pond + Pond:Treatment,
              data=FT_DC_P1)

Anova(model_FT_P1, type="II") 




FT_DC_P2 <- FT_DC[which(FT_DC$Treatment =='T1_0.45'|
                          FT_DC$Treatment =='T4_Sterivex'),]


model_FT_P2 = lm(transformTukey(Time) ~ Treatment + Pond + Pond:Treatment,
                 data=FT_DC_P2)

Anova(model_FT_P2, type="II") 




FT_DC_P3 <- FT_DC[which(FT_DC$Treatment =='T1_0.45'|
                          FT_DC$Treatment =='T5_PF_0.45'),]

model_FT_P3 = lm(transformTukey(Time) ~ Treatment + Pond + Pond:Treatment,
                 data=FT_DC_P3)

Anova(model_FT_P3, type="II") 


####DNA####

model_DNA = lm(sqrt(DNA) ~ Treatment + Pond + Pond:Treatment,
              data=FT_DC)

Anova(model_DNA, type="II") 


FT_DC_P0 <- FT_DC[which(FT_DC$Treatment !='TP5_PF'),]


model_DNA_0 = lm(sqrt(DNA) ~ Treatment + Pond + Pond:Treatment,
                data=FT_DC_P0)

Anova(model_DNA_0, type="II") 



FT_DC_P1 <- FT_DC[which(FT_DC$Treatment =='T1_0.45'|
                          FT_DC$Treatment =='T2_0.8'|
                          FT_DC$Treatment =='T3_1.2'),]

model_DNA_P1 = lm(sqrt(DNA) ~ Treatment + Pond + Pond:Treatment,
                 data=FT_DC_P1)

Anova(model_DNA_P1, type="II") 




FT_DC_P2 <- FT_DC[which(FT_DC$Treatment =='T1_0.45'|
                          FT_DC$Treatment =='T4_Sterivex'),]


model_DNA_P2 = lm(sqrt(DNA) ~ Treatment + Pond + Pond:Treatment,
                 data=FT_DC_P2)

Anova(model_DNA_P2, type="II") 




FT_DC_P3 <- FT_DC[which(FT_DC$Treatment =='T1_0.45'|
                          FT_DC$Treatment =='T5_PF_0.45'),]

model_DNA_P3 = lm(sqrt(DNA) ~ Treatment + Pond + Pond:Treatment,
                 data=FT_DC_P3)

Anova(model_DNA_P3, type="II") 






####Species detection####

Total_sepcies <- summarySE(Filtration_rs_sample, measurevar="Reads",groupvars=c("SampleID"))

Total_sepcies <- Total_sepcies [c('SampleID','N')]
Total_sepcies <- rename(Total_sepcies, c('N'='TotalN'))



Filtration_rs_sample_sepcies <- Filtration_rs_sample[which(Filtration_rs_sample$Reads>0),]

Num_sp <- summarySE(Filtration_rs_sample_sepcies, measurevar="Reads",groupvars=c("SampleID"))


Num_sp <- Num_sp [c('SampleID','N')]

Num_sp <- rbind(Num_sp, c('T4-1-3',0))




FT_DC <- merge(FT_DC,Num_sp, by=c("SampleID"))

FT_DC <- merge(FT_DC,Total_sepcies, by=c("SampleID"))

FT_DC$N <- as.numeric(FT_DC$N)

FT_DC$Rate <- round(FT_DC$N/FT_DC$TotalN,2)


shapiro.test(transformTukey(FT_DC$Rate))




model_Rate = lm(transformTukey(Rate) ~ Treatment + Pond + Pond:Treatment,
               data=FT_DC)

Anova(model_Rate, type="II") 


FT_DC_P0 <- FT_DC[which(FT_DC$Treatment !='TP5_PF'),]


model_Rate_0 = lm(transformTukey(Rate) ~ Treatment + Pond + Pond:Treatment,
                 data=FT_DC_P0)

Anova(model_Rate_0, type="II") 



FT_DC_P1 <- FT_DC[which(FT_DC$Treatment =='T1_0.45'|
                          FT_DC$Treatment =='T2_0.8'|
                          FT_DC$Treatment =='T3_1.2'),]

model_Rate_P1 = lm(transformTukey(Rate) ~ Treatment + Pond + Pond:Treatment,
                  data=FT_DC_P1)

Anova(model_Rate_P1, type="II") 




FT_DC_P2 <- FT_DC[which(FT_DC$Treatment =='T1_0.45'|
                          FT_DC$Treatment =='T4_Sterivex'),]


model_Rate_P2 = lm(transformTukey(Rate) ~ Treatment + Pond + Pond:Treatment,
                  data=FT_DC_P2)

Anova(model_Rate_P2, type="II") 




FT_DC_P3 <- FT_DC[which(FT_DC$Treatment =='T1_0.45'|
                          FT_DC$Treatment =='T5_PF_0.45'),]

model_Rate_P3 = lm(transformTukey(Rate) ~ Treatment + Pond + Pond:Treatment,
                  data=FT_DC_P3)

Anova(model_Rate_P3, type="II") 





ddply(FT_DC, .(Pond,Treatment), summarize,
                       Rate_mean = round(mean(Rate), 2),
                       Rate_sd = round(sd(Rate), 2))


ddply(FT_DC, .(Treatment), summarize,
      Rate_mean = round(mean(Rate), 2),
      Rate_sd = round(sd(Rate), 2))

###spadeR#####
spadeR_all_with <- spadeR_all[which(spadeR_all$Group == 'Within'),]



qqnorm(transformTukey(spadeR_all_with$Estimate,plotit=FALSE))

model_SpadeR = lm(transformTukey(Estimate) ~ Treatment + Pond + Pond:Treatment,
                  data=spadeR_all_with)

Anova(model_SpadeR, type="II") 


spadeR_all_with_P0 <- spadeR_all_with[which(spadeR_all_with$Treatment !='TP5_PF'),]


model_SpadeR_0 = lm(transformTukey(Estimate) ~ Treatment + Pond + Pond:Treatment,
                    data=spadeR_all_with_P0)

Anova(model_SpadeR_0, type="II") 



spadeR_all_with_P1 <- spadeR_all_with[which(spadeR_all_with$Treatment =='T1_0.45'|
                                              spadeR_all_with$Treatment =='T2_0.8'|
                                              spadeR_all_with$Treatment =='T3_1.2'),]

model_SpadeR_P1 = lm(transformTukey(Estimate) ~ Treatment + Pond + Pond:Treatment,
                     data=spadeR_all_with_P1)

Anova(model_SpadeR_P1, type="II") 




spadeR_all_with_P2 <- spadeR_all_with[which(spadeR_all_with$Treatment =='T1_0.45'|
                                              spadeR_all_with$Treatment =='T4_Sterivex'),]


model_SpadeR_P2 = lm(transformTukey(Estimate) ~ Treatment + Pond + Pond:Treatment,
                     data=spadeR_all_with_P2)

Anova(model_SpadeR_P2, type="II") 




spadeR_all_with_P3 <- spadeR_all_with[which(spadeR_all_with$Treatment =='T1_0.45'|
                                              spadeR_all_with$Treatment =='T5_PF_0.45'),]

model_SpadeR_P3 = lm(transformTukey(Estimate)~ Treatment + Pond + Pond:Treatment,
                     data=spadeR_all_with_P3)

Anova(model_SpadeR_P3, type="II") 





#############correlation#######


Filtration_summary_rep <-merge(Filtration_rs_sample_outlier, Pond_summary_Aug, by=c("Pond","Species"))


for (c in 1:nrow(Filtration_summary_rep)) {
  if(Filtration_summary_rep$Reads[c] == 0){Filtration_summary_rep$Reads[c] <- 1} 
}



Filtration_summary_rep$lgReads <- log10(Filtration_summary_rep$Reads)
Filtration_summary_rep$lgNum <- log10(Filtration_summary_rep$Number)
Filtration_summary_rep$lgBio <- log10(Filtration_summary_rep$Biomass)


correlation_rep_num <- ddply(Filtration_summary_rep, .(Pond,Treatment,SampleID), summarise, s=corfun(lgNum,lgReads)$statistic,
                        cor.est=corfun(lgNum,lgReads)$estimate,
                        pval=corfun(lgNum,lgReads)$p.value,
                        alt=corfun(lgNum,lgReads)$alternative,
                        df=corfun(lgNum,lgReads)$parameter)


correlation_rep_num[is.na(correlation_rep_num)] <- 0



qqnorm(transformTukey(correlation_rep_num$cor.est,plotit=FALSE))


ddply(correlation_rep_num, .(Pond,Treatment), summarize,
                       cor_mean = round(mean(cor.est), 2))


model_num = lm(transformTukey(cor.est) ~ Treatment + Pond + Pond:Treatment,
                  data=correlation_rep_num)

Anova(model_num, type="II") 


correlation_rep_num_P0 <- correlation_rep_num[which(correlation_rep_num$Treatment !='TP5_PF'),]


model_num_0 = lm(transformTukey(cor.est) ~ Treatment + Pond + Pond:Treatment,
                    data=correlation_rep_num_P0)

Anova(model_num_0, type="II") 



correlation_rep_num_P1 <- correlation_rep_num[which(correlation_rep_num$Treatment =='T1_0.45'|
                                              correlation_rep_num$Treatment =='T2_0.8'|
                                              correlation_rep_num$Treatment =='T3_1.2'),]

model_num_P1 = lm(transformTukey(cor.est) ~ Treatment + Pond + Pond:Treatment,
                     data=correlation_rep_num_P1)

Anova(model_num_P1, type="II") 



correlation_rep_num_P4 <- correlation_rep_num[which(correlation_rep_num$Treatment =='T1_0.45'|
                                                      correlation_rep_num$Treatment =='T3_1.2'),]

model_num_P4 = lm(transformTukey(cor.est) ~ Treatment + Pond + Pond:Treatment,
                  data=correlation_rep_num_P4)

Anova(model_num_P4, type="II") 





correlation_rep_num_P2 <- correlation_rep_num[which(correlation_rep_num$Treatment =='T1_0.45'|
                                              correlation_rep_num$Treatment =='T4_Sterivex'),]


model_num_P2 = lm(transformTukey(cor.est) ~ Treatment + Pond + Pond:Treatment,
                     data=correlation_rep_num_P2)

Anova(model_num_P2, type="II") 




correlation_rep_num_P3 <- correlation_rep_num[which(correlation_rep_num$Treatment =='T1_0.45'|
                                              correlation_rep_num$Treatment =='T5_PF_0.45'),]

model_num_P3 = lm(transformTukey(cor.est)~ Treatment + Pond + Pond:Treatment,
                     data=correlation_rep_num_P3)

Anova(model_num_P3, type="II") 






correlation_rep_bio <- ddply(Filtration_summary_rep, .(Pond,Treatment,SampleID), summarise, s=corfun(lgBio,lgReads)$statistic,
                                  cor.est=corfun(lgBio,lgReads)$estimate,
                                  pval=corfun(lgBio,lgReads)$p.value,
                                  alt=corfun(lgBio,lgReads)$alternative,
                                  df=corfun(lgBio,lgReads)$parameter)

correlation_rep_bio[is.na(correlation_rep_bio)] <- 0

qqnorm(transformTukey(correlation_rep_bio$cor.est,plotit=FALSE))


ddply(correlation_rep_bio, .(Pond,Treatment), summarize,
      cor_mean = round(mean(cor.est), 2))


model_bio = lm(transformTukey(cor.est) ~ Treatment + Pond + Pond:Treatment,
               data=correlation_rep_bio)

Anova(model_bio, type="II") 


correlation_rep_bio_P0 <- correlation_rep_bio[which(correlation_rep_bio$Treatment !='TP5_PF'),]


model_bio_0 = lm(transformTukey(cor.est) ~ Treatment + Pond + Pond:Treatment,
                 data=correlation_rep_bio_P0)

Anova(model_bio_0, type="II") 



correlation_rep_bio_P1 <- correlation_rep_bio[which(correlation_rep_bio$Treatment =='T1_0.45'|
                                                      correlation_rep_bio$Treatment =='T2_0.8'|
                                                      correlation_rep_bio$Treatment =='T3_1.2'),]

model_bio_P1 = lm(transformTukey(cor.est) ~ Treatment + Pond + Pond:Treatment,
                  data=correlation_rep_bio_P1)

Anova(model_bio_P1, type="II") 



correlation_rep_bio_P4 <- correlation_rep_bio[which(correlation_rep_bio$Treatment =='T2_0.8'|
                                                      correlation_rep_bio$Treatment =='T3_1.2'),]

model_bio_P4 = lm(transformTukey(cor.est) ~ Treatment + Pond + Pond:Treatment,
                  data=correlation_rep_bio_P4)

Anova(model_bio_P4, type="II") 





correlation_rep_bio_P2 <- correlation_rep_bio[which(correlation_rep_bio$Treatment =='T1_0.45'|
                                                      correlation_rep_bio$Treatment =='T4_Sterivex'),]


model_bio_P2 = lm(transformTukey(cor.est) ~ Treatment + Pond + Pond:Treatment,
                  data=correlation_rep_bio_P2)

Anova(model_bio_P2, type="II") 




correlation_rep_bio_P3 <- correlation_rep_bio[which(correlation_rep_bio$Treatment =='T1_0.45'|
                                                      correlation_rep_bio$Treatment =='T5_PF_0.45'),]

model_bio_P3 = lm(transformTukey(cor.est)~ Treatment + Pond + Pond:Treatment,
                  data=correlation_rep_bio_P3)

Anova(model_bio_P3, type="II") 


