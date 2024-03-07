# This code is to replicate analysis and figures for the paper "Preservation of agri-food byproducts by acidification and fermentation in black soldier fly larvae bioconversion", Alciatore et. al 2024
# Code developed by Giacomo Alciatore





library(pracma)

library(corrr)
library(ggcorrplot)
library(ggthemes)
library(FactoMineR)
library(factoextra)


library(car)
library(sjPlot)
library(gridExtra) 
library(gtable)
library(grid)  


library(RColorBrewer)
library(emmeans)
library(vegan)



library(phyloseq)
library(ggrepel)


# Necessary packages
library(scales)
library(rje)
library(xlsx)
library(ggplot2)
library(tidyr)
library(cowplot)
library(dplyr)
library(pheatmap)
library(lubridate)
library(ampvis2)
library(adegenet)
library(ggforce)
# Define colors, shapes and parameters for the whole study 
# Colors are picked to be clearly different also on greyscale

show_col(cubeHelix(9,start = 0, r = -1, hue = 3, gamma = 1))
cubeHelix(9,start = 0, r = -1, hue = 3, gamma = 1)
ACD = "#00AC00"
CTR = "#FF7822"
FRM = "#FFBCFF"
PRE = "#00325A" 

D1 = "#008140"
D7 = "#8E9B00"
D14 = "#FF7EBF"

sAC = 21 # Shapes 
sWF = 22
sIF = 23
sPR = 24

NcP = 4.67 # Nitrogen Protein conversion factor

# WDir<-choose.dir() #choose folder of the project (alternatively write path)
WDir <- "H:/.shortcut-targets-by-id/1KyBOLIX6Yss5hEEQ9-lJ3Q-Jg0LoRIbU/Fermentation study - Giacomo"

## Import and prepare data ##

# Rearing data
Data <- read.xlsx(file.path(WDir,"2023_05_30_MainExp/DataExp.xlsx"),sheetName="Experiment",stringsAsFactors=TRUE)
Eliminate <- Data$sample[Data$day==8&Data$survival<0.51] # eliminate replicates with too few larva at the beginning of the Exp
Data <- subset(Data, sample!=Eliminate)
dataDay0 <- subset(Data, !duplicated(sample)) # create new rows for the beginning of the experiment, to use in plots
dataDay0$day <- 0
dataDay0$avg_mass_larva_mg <- dataDay0$larvae_start 
Data <- rbind(Data,dataDay0) # merge data together
DataH <- subset(Data, day==8) # only data from Harvest day

# HPLC data
HPLC <- read.xlsx(file.path(WDir,"2023_05_30_MainExp/LabAnalysis/HPLC_Pharmabiome/HPLC.xlsx"),sheetIndex=1,stringsAsFactors=TRUE)
HPLC$int <- interaction(HPLC$treat,HPLC$preserv)
HPLC$ID <- row.names(HPLC)

# Microbial data
FileList <- read.xlsx(file.path(WDir,"2023_05_30_MainExp/LabAnalysis/Microbial_Community/Old+New/NamingMicrobialCommunity.xlsx"),sheetIndex=1,stringsAsFactors=TRUE)

Files_feed <- subset(FileList,type=="feed") # Only files with feedstock data
Files_feed <- Files_feed[with(Files_feed, order(day, preserv)),]

for (file in 1:nrow(Files_feed)){ 
  name.file = file.path(WDir,"2023_05_30_MainExp/LabAnalysis/Microbial_Community/Old+New/Raw_Data",as.character(Files_feed$file)[file])
  
  d <- read.delim(name.file,sep="\t")
  t <- d[,c("tax_id","lineage")]
  f <- d[,1:2]
  ec <- d[,c(1,4)]
  colnames(f)[2] <- as.character(Files_feed$name[file])
  colnames(ec)[2] <- as.character(Files_feed$name[file])
  
  
  t$lineage = substr(t$lineage,1,nchar(t$lineage)-1)
  t2<-t |>
    separate_wider_delim(lineage, delim = ";", names = c("Kingdom", "Phylum","Class","Order","Family","Genus","Species"),too_few="align_start",too_many="merge")
  
  if (file==1){df = f
  tree = t2
  ECC = ec} else {
    df = merge(df,f,by="tax_id",no.dups=TRUE,all=TRUE)
    tree = rbind(tree,subset(t2, !(tax_id %in% tree$tax_id)))
    ECC = merge(ECC,ec,by="tax_id",no.dups=TRUE,all=TRUE)
  }
}

colnames(df)[1] <- "OTU"
colnames(tree)[1] <- "OTU"
colnames(ECC)[1] <- "OTU"

ECC <- ECC %>% 
  mutate_if(is.numeric, round)

SumCounts <- rowSums(ECC[,2:9],na.rm=T)

ECC <- ECC[SumCounts>10,] # For next steps, only select taxas with at least 10 counts
df <- df[SumCounts>10,]

Microb <- merge(df,tree,by=c("OTU")) 

# Temperature data

T <- read.xlsx(file.path(WDir,"2023_05_30_MainExp/Environmental_Data/TempSubstrate.xlsx"), sheetName = "Sheet1", stringsAsFactors = TRUE)
Tcontainer <- read.xlsx(file.path(WDir, "2023_05_30_MainExp/Environmental_Data/TempContainer.xlsx"), sheetName = "Sheet1", stringsAsFactors = TRUE)
Start <- read.xlsx(file.path(WDir, "2023_05_30_MainExp/StartExp.xlsx"), sheetName="Sheet1", stringsAsFactors =TRUE)
Tcontainer$dash = "Dashed"

Tlong <- T %>% 
  pivot_longer(
    cols = c(2:10), 
    names_to = "ID",
    values_to = "Temp"
  )

Tlong = drop_na(Tlong)

Tlong$Treat = as.factor(substring(Tlong$ID, first = 1, last = 3))
Tlong$Exp = as.factor(substring(Tlong$ID, first = 5, last = 6))

Temp1D<- Tlong[Tlong$Exp=='1D' & Tlong$Time>=Start$start[Start$Exp=="1D"] & Tlong$Time<=Start$end[Start$Exp=="1D"],]
Temp1W<- Tlong[Tlong$Exp=='1W' & Tlong$Time>=Start$start[Start$Exp=="1W"] & Tlong$Time<=Start$end[Start$Exp=="1W"],]
Temp2W<- Tlong[Tlong$Exp=='2W' & Tlong$Time>=Start$start[Start$Exp=="2W"] & Tlong$Time<=Start$end[Start$Exp=="2W"],]

## Figures ##
 
# Figure 1 - Larval growth & Temperature

# Fig.1ABC
WWmax=max(Data$avg_mass_larva_mg)

p1 <- ggplot(data=Data[Data$treat == "1D",]) + 
  geom_jitter(aes(x = day,y = avg_mass_larva_mg, fill = preserv, shape = preserv), width = 0.15, size = 1.5, alpha = 0.7, stroke = 1) +
  ylim(0,WWmax) + stat_summary(aes(day, avg_mass_larva_mg,col = preserv), fun = mean,geom = "line", linewidth = .8, linetype = "dashed", show.legend=FALSE) +
  stat_summary(aes(day,avg_mass_larva_mg,fill=preserv,shape=preserv), fun=mean,geom="point",size=2,stroke=.8,show.legend=FALSE)+
  labs(y = "Avg larval weight (mg)", x = "Rearing duration (days)") + ggtitle("A   1-day storage") +
  scale_color_manual(values = c("ACD" = ACD,"CTR" = CTR,"FRM" = FRM), labels = c("AC","WF","IF")) +
  scale_shape_manual(name = "Treatment", values = c("ACD" = sAC, "CTR" = sWF, "FRM" = sIF))+
  scale_fill_manual(name = "Treatment", values = c("ACD" = ACD, "CTR" = CTR, "FRM" = FRM))+
  scale_x_continuous(breaks = seq(0,8,2)) + theme_minimal() + theme(legend.position = 'none') 

p2 <- ggplot(data=Data[Data$treat == "1W",]) + 
  geom_jitter(aes(x = day,y = avg_mass_larva_mg, fill = preserv, shape = preserv), width = 0.15, size = 1.5, alpha = 0.7, stroke = 1) +
  ylim(0,WWmax) + stat_summary(aes(day, avg_mass_larva_mg,col = preserv), fun = mean,geom = "line", linewidth = .8, linetype = "dashed", show.legend=FALSE) +
  stat_summary(aes(day,avg_mass_larva_mg,fill=preserv,shape=preserv), fun=mean,geom="point",size=2,stroke=.8,show.legend=FALSE)+
  labs(y = "", x = "Rearing duration (days)") + ggtitle("B   7-day storage") +
  scale_color_manual(values = c("ACD" = ACD,"CTR" = CTR,"FRM" = FRM), labels = c("AC","WF","IF")) +
  scale_shape_manual(name = "Treatment", values = c("ACD" = sAC, "CTR" = sWF, "FRM" = sIF))+
  scale_fill_manual(name = "Treatment", values = c("ACD" = ACD, "CTR" = CTR, "FRM" = FRM))+
  scale_x_continuous(breaks = seq(0,8,2)) + theme_minimal() + theme(legend.position = 'none') 

p3 <- ggplot(data=Data[Data$treat == "2W",]) + 
  geom_jitter(aes(x = day,y = avg_mass_larva_mg, fill = preserv, shape = preserv), width = 0.15, size = 1.5, alpha = 0.7, stroke = 1) +
  ylim(0,WWmax) + stat_summary(aes(day, avg_mass_larva_mg,col = preserv), fun = mean,geom = "line", linewidth = .8, linetype = "dashed", show.legend=FALSE) +
  stat_summary(aes(day,avg_mass_larva_mg,fill=preserv,shape=preserv), fun=mean,geom="point",size=2,stroke=.8)+
  labs(y = "", x = "Rearing duration (days)") + ggtitle("C   14-day storage") +
  scale_color_manual(values = c("ACD" = ACD,"CTR" = CTR,"FRM" = FRM), labels = c("AC","WF","IF")) +
  scale_shape_manual(name = "Treatment", values = c("ACD" = sAC, "CTR" = sWF, "FRM" = sIF), labels = c("AC", "WF", "IF"))+
  scale_fill_manual(name = "Treatment", values = c("ACD" = ACD, "CTR" = CTR, "FRM" = FRM), labels = c("AC", "WF", "IF"))+
  scale_x_continuous(breaks = seq(0,8,2)) + theme_minimal()

legend <- get_legend(p3 + theme_minimal(base_size = 10)) 

theme_set(theme_cowplot(font_size = 12) + 
            theme(text = element_text(size = 12)))

P <- cowplot::plot_grid(p1, NULL, p2, NULL, p3 + theme(legend.position = 'none'), legend, ncol = 6, rel_widths = c(1,-.1,1,-.1,1, .3))
file_name = paste(WDir,"/Manuscript/Plots/Fig1_ABC.tiff", sep = "")
tiff(file_name, units = "mm", width = 182, height = 60, res = 300)
print(P)
dev.off()


# Fig.1DEF
T1 <- ggplot() + geom_line(data = Temp1D, aes(Time, Temp, group = Treat, colour = Treat), linewidth = 1.2, show.legend = FALSE) +
  geom_line(data = Tcontainer[Tcontainer$Date>=Start$start[Start$Exp == "1D"] & Tcontainer$Date<=Start$end[Start$Exp == "1D"],], 
            aes(Date, Temperature, col = dash), linewidth = 1.2, show.legend = F) + theme_minimal() +
  scale_color_manual(values = c("ACD" = ACD,"CTR" = CTR,"FRM" = FRM,"Dashed" = "black")) +
  labs(y = "Temperature (°C)",x = "Rearing duration (days)") + ggtitle("D   1-day storage") + theme(legend.position="none") + ylim(20,50) + 
  scale_x_datetime(breaks = as.POSIXct(c("2023-06-06","2023-06-08","2023-06-10","2023-06-12", "2023-06-14")), labels = c("0","2","4","6","8 "))

T2 <- ggplot() + geom_line(data = Temp1W, aes(Time, Temp, group = Treat, colour = Treat), linewidth = 1.2, show.legend = FALSE) +
  geom_line(data = Tcontainer[Tcontainer$Date>=Start$start[Start$Exp == "1W"] & Tcontainer$Date<=Start$end[Start$Exp == "1W"],], 
            aes(Date, Temperature, col = dash), linewidth = 1.2, show.legend = F) + theme_minimal() +
  scale_color_manual(values = c("ACD" = ACD,"CTR" = CTR,"FRM" = FRM,"Dashed" = "black")) +
  labs(y = "",x = "Rearing duration (days)") + ggtitle("E   7-day storage") + theme(legend.position="none") + ylim(20,50) + 
  scale_x_datetime(breaks = as.POSIXct(c("2023-06-19","2023-06-21","2023-06-23","2023-06-25", "2023-06-27")),labels = c("0","2","4","6","8 "))

T3 <- ggplot() + geom_line(data = Temp2W, aes(Time, Temp, group = Treat, colour = Treat), linewidth = 1.2) +
  geom_line(data = Tcontainer[Tcontainer$Date>=Start$start[Start$Exp == "2W"] & Tcontainer$Date<=Start$end[Start$Exp == "2W"],], 
            aes(Date, Temperature, col = dash), linewidth = 1.2, show.legend = F) + theme_minimal() +
  scale_color_manual(values = c("ACD" = ACD,"CTR" = CTR,"FRM" = FRM,"Dashed" = "black"),labels = c("ACD" = "AC", "CTR" = "WF", "FRM" = "IF", "Dashed" = "Ambient\ntemperature")) +
  labs(y = "",x = "Rearing duration (days)") + ggtitle("F   14-day storage") + ylim(20,50) + 
  scale_x_datetime(breaks = as.POSIXct(c("2023-06-13","2023-06-15","2023-06-17","2023-06-19", "2023-06-21")),labels = c("0","2","4","6","8 "))

legend <- get_legend(T3 + theme_minimal(base_size = 10)) 

theme_set(theme_cowplot(font_size = 12) + 
            theme(text = element_text(size = 12)))

T<-cowplot::plot_grid(T1, NULL, T2, NULL, T3 + theme(legend.position = 'none'), legend, ncol = 6, rel_widths = c(1, -.15, 1, -.15, 1, .35))
file_name = paste(WDir,"/Manuscript/Plots/Fig1_DEF.tiff", sep = "")
tiff(file_name, units = "mm", width = 182, height = 60, res = 300)
print(T)
dev.off()

# Figure 2 - Heatmap of HPLC

HPLC$Sample <- paste0(HPLC$name, " ", HPLC$treat)
HPLC$name <- factor(HPLC$name, levels = c("Pre", "AC", "IF", "WF"))
HPLC$treat <- factor(HPLC$treat, levels=c("1D", "7D", "14D"))

HPLC_heat <- HPLC %>%
  group_by(Sample, name, treat) %>%
  summarise(across(7:21, list(mean = mean), .names = "{.col}"))

HPLC_sd <- HPLC %>%
  group_by(Sample, name, treat) %>%
  summarise(across(7:21, list(sd = sd), .names = "{.col}"))


col_tre = as.data.frame(HPLC_heat[, 2:3])
rownames(col_tre) <- HPLC_heat$Sample
colnames(col_tre) <- c("Treatment", "Duration")

HPLC_heat <- HPLC_heat[with(HPLC_heat, order(treat, name)),]
pHPLC <- as.data.frame(t(HPLC_heat[, c(4:11)]))

colnames(pHPLC) <- HPLC_heat$Sample

pHPLC <- pHPLC[order(rowSums(pHPLC), decreasing=TRUE), ]
pHPLC <- pHPLC[1:6, ]

colfunc_grey <- colorRampPalette(c("#FEFFBF", "#FC8D59"))

Fig2 <- pheatmap(pHPLC, treeheight_row = 0, treeheight_col = 0, angle_col = 315, labels_col = HPLC_heat$Sample, display_numbers = round(pHPLC),
               annotation_col = col_tre, color = colfunc_grey(100), fontsize = 12,
               annotation_colors = list(Treatment = c(Pre = PRE,AC = ACD,IF = FRM,WF = CTR),Duration=c("1D" = D1,"7D" = D7,"14D" = D14)),
               cluster_cols = FALSE,cluster_rows = FALSE)


file_name = paste(WDir,"/Manuscript/Plots/Fig_2.tiff", sep= "")
tiff(file_name, units = "mm", width = 182, height = 91, res = 300)
print(Fig2)
dev.off()

# Figure 3 - Heatmaps of microbial communities

Amp <- amp_load(otutable = df,taxonomy = tree, metadata = Files_feed)
Amp_counts<- amp_load(otutable = ECC, taxonomy = tree, metadata = Files_feed)


Amp_Filt = amp_filter_taxa(Amp, tax_vector = c("Cyanobacteria","Mitochondria"), normalise = T, remove = TRUE)
Amp_counts_F=amp_filter_taxa(Amp_counts,tax_vector = c("Cyanobacteria", "Mitochondria"), normalise = F, remove = TRUE)


Amp_Filt$metadata$name <- factor(Amp_Filt$metadata$name, levels=c("AC1D", "IF1D", "WF1D", "AC7D", "IF7D", "WF7D", "AC14D", "IF14D", "WF14D"))

G_MC<-amp_heatmap(Amp_Filt,
                  tax_aggregate = "Genus",
                  tax_show = 15,
                  plot_na = TRUE,
                  normalise = TRUE,
                  group_by = "name",
                  color = colfunc_grey(100))+
  ggtitle("B") +
  theme(text = element_text(size = 12), axis.text.x = element_text(angle = 45, vjust = 1),
        axis.text.y = element_text(face = "italic"),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0, vjust = 0))

file_name = paste(WDir, "/Manuscript/Plots/Fig3_B.tiff", sep = "")
tiff(file_name, units = "mm", width = 182, height = 91, res = 300)
print(G_MC)
dev.off()

P_MC<-amp_heatmap(Amp_Filt,
                  tax_aggregate = "Phylum",
                  tax_show = 5,
                  plot_na = TRUE,
                  normalise = TRUE,
                  group_by = "name",
                  color = colfunc_grey(100))+
  ggtitle("A") + labs(fill = "Relative Abundance") +
    theme(text = element_text(size = 12), axis.text.x = element_text(angle = 45, vjust = 1),
        axis.text.y = element_text(face = "italic"),
        legend.position = "bottom", legend.justification = "center", legend.title = element_text(size = 12, margin = margin(r = 20)))

file_name = paste(WDir, "/Manuscript/Plots/Fig3_A.tiff", sep = "")
tiff(file_name, units = "mm", width = 182, height = 91, res = 300)
print(P_MC)
dev.off()

# Microbial diversity

AAD<-amp_alpha_diversity(Amp_counts_F, measure = c("observed","shannon","simpson"), richness = TRUE, rarefy = NULL)
write.xlsx2(AAD,file.path(WDir,"2023_05_30_MainExp/MicrobialDiversity.xlsx") , sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)

# Figure 4 - DPCA Metabolites

metab <- row.names(pHPLC)#Use the same 6 metabolites that are in the heatmap
m <- HPLC[,metab] 
norm_m <- scale(m)

dapc1 <- dapc(norm_m,HPLC$int)  #PC=2 DF=2
tab = as.data.frame(dapc1$ind.coord)
tab$ID <- row.names(tab)
dpca_HPLC <- merge(x = tab, y=HPLC, by=c("ID"))

maxX = max(dpca_HPLC$LD1)*1.2
minX = min(dpca_HPLC$LD1)*1.2
maxY = max(dpca_HPLC$LD2)*1.2
minY = min(dpca_HPLC$LD2)*1.2

DPCA <- ggplot(dpca_HPLC, aes(x = LD1, y = LD2, fill = preserv, shape = preserv))+
  geom_point(size=4,stroke=2)+geom_mark_ellipse(aes(x=LD1,y=LD2,fill=preserv,group=int),con.cap=0,show.legend=FALSE)+
  labs(x="DPC1 (57.5%)", y="DPC2 (29.6%)")+theme_minimal(base_size = 20)+ylim(minY,maxY)+xlim(minX,maxX)+
  theme(axis.text=element_text(face="bold"),legend.title=element_blank())+
  scale_fill_manual(values=c("ACD"=ACD,"CTR"=CTR,"FRM"=FRM,"PRE"=PRE),labels=c("ACD"="AC","CTR"="WF","FRM"="IF","PRE"="Pre"))+
  scale_shape_manual(values=c("ACD"=sAC,"CTR"=sWF,"FRM"=sIF,"PRE"=sPR),labels=c("ACD"="AC","CTR"="WF","FRM"="IF","PRE"="Pre"))

file_name = paste(WDir,"/Manuscript/Plots/HPCL.pdf", sep="") #Unfortunately, best way I found to make labels legible was to add them with a design software 
pdf(file_name, width=8, height=4.5)
print(DPCA)
dev.off()

## Statistical results as reported in paper ##

shapiro.test(DataH$avg_dry_larva)
hist(DataH$avg_dry_larva)
leveneTest(DataH$avg_dry_larva, group = DataH$preserv)

LDm_1D<-aov(larvae_dry ~ preserv, data=DataH[DataH$treat=="1D",])
summary(LDm_1D)
plot(LDm_1D)
TukeyHSD(LDm_1D)

LDm_1W<-aov(larvae_dry ~ preserv, data=DataH[DataH$treat=="1W",])
summary(LDm_1W)
TukeyHSD(LDm_1W)

LDm_2W<-aov(larvae_dry ~ preserv, data=DataH[DataH$treat=="2W",])
summary(LDm_2W)
TukeyHSD(LDm_2W)



emmeans(m1, specs = pairwise ~ Treatment|Age)
emmeans(m1, specs = pairwise ~ Age)

LDm<-aov(avg_dry_larva ~ preserv*treat, data=DataH)
summary(LDm)
emmeans(LDm,specs=pairwise ~ preserv|treat)

CRm<-aov(dry_ConversionRate ~ preserv*treat, data=DataH)
summary(CRm)
emmeans(CRm,specs=pairwise ~ preserv|treat)

WRm<-aov(dry_WasteReduction ~ preserv*treat, data=DataH)
summary(WRm)
emmeans(WRm,specs=pairwise ~ preserv|treat)

DataH$survival[DataH$survival > 1] <- 1
Sm<-aov(survival ~ preserv*treat, data=DataH)
summary(Sm)
emmeans(Sm,specs=pairwise ~ preserv|treat)

PPm<-aov(PP ~ preserv*treat, data=DataH)
summary(PPm)
emmeans(PPm,specs=pairwise ~ preserv|treat)

Ashm<-aov(ADry ~ preserv*treat, data=NCFH)
summary(Ashm)
emmeans(Ashm,specs=pairwise ~ preserv|treat)

Fatm<-aov(Fdrydry ~ preserv*treat, data=NCFH)
summary(Fatm)
emmeans(Fatm,specs=pairwise ~ preserv|treat)

Protm<-aov(PDry ~ preserv*treat, data=NCFH)
summary(Protm)
emmeans(Protm,specs=pairwise ~ preserv|treat)


#Wet Weight treat 1D day 4
WWm1Dday4<-aov(avg_mass_larva_mg ~ preserv, data=Data[Data$treat=="1D"&Data$day==4,])
summary(WWm1Dday4)
emmeans(WWm1Dday4,specs=pairwise ~ preserv)

#Wet Weight treat 1D day 6
WWm1Dday6<-aov(avg_mass_larva_mg ~ preserv, data=Data[Data$treat=="1D"&Data$day==6,])
summary(WWm1Dday6)
emmeans(WWm1Dday6,specs=pairwise ~ preserv)

#Wet Weight treat 1D day 8
WWm1Dday8<-aov(avg_mass_larva_mg ~ preserv, data=Data[Data$treat=="1D"&Data$day==8,])
summary(WWm1Dday8)
emmeans(WWm1Dday8,specs=pairwise ~ preserv)

#Wet Weight treat 1W day 4
WWm1Wday4<-aov(avg_mass_larva_mg ~ preserv, data=Data[Data$treat=="1W"&Data$day==4,])
summary(WWm1Wday4)
emmeans(WWm1Wday4,specs=pairwise ~ preserv)

#Wet Weight treat 1W day 6
WWm1Wday6<-aov(avg_mass_larva_mg ~ preserv, data=Data[Data$treat=="1W"&Data$day==6,])
summary(WWm1Wday6)
emmeans(WWm1Wday6,specs=pairwise ~ preserv)

#Wet Weight treat 1W day 8
WWm1Wday8<-aov(avg_mass_larva_mg ~ preserv, data=Data[Data$treat=="1W"&Data$day==8,])
summary(WWm1Wday8)
emmeans(WWm1Wday8,specs=pairwise ~ preserv)


#Wet Weight treat 2W day 4
WWm2Wday4<-aov(avg_mass_larva_mg ~ preserv, data=Data[Data$treat=="2W"&Data$day==4,])
summary(WWm2Wday4)
emmeans(WWm2Wday4,specs=pairwise ~ preserv)

#Wet Weight treat 2W day 6
WWm2Wday6<-aov(avg_mass_larva_mg ~ preserv, data=Data[Data$treat=="2W"&Data$day==6,])
summary(WWm2Wday6)
emmeans(WWm2Wday6,specs=pairwise ~ preserv)

#Wet Weight treat 2W day 8
WWm2Wday8<-aov(avg_mass_larva_mg ~ preserv, data=Data[Data$treat=="2W"&Data$day==8,])
summary(WWm2Wday8)
emmeans(WWm2Wday8,specs=pairwise ~ preserv)