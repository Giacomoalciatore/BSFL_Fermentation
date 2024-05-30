# This code is to replicate analysis and figures for the paper "Preservation of agri-food byproducts by acidification and fermentation in black soldier fly larvae bioconversion", Alciatore et. al 2024
# Code developed by Giacomo Alciatore

#Packages
library(emmeans)
library(car)
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
library(ggthemes)

# Define colors, shapes and parameters for the publication 

AC = "#cc1122" 
WF = "#000ca7"
IF = "#f1a226"
PRE = "#9d02d7" 

D1 = "#298c8c"
D7 = "#b8b8b8"
D14 = "#5e4c5f"

micro = c("#9E5C00", "#7D3200", "#74c476" , "#238b45", "#fdae6b", "#ff7f00", "#4292c6", "#9ecae1", "#6a51a3", "#9e9ac8", "#969696", "#525252")

sAC = 21 # Shapes 
sWF = 22
sIF = 23
sPR = 24

NcP = 4.67 # Nitrogen Protein conversion factor

# WDir<-choose.dir() #choose folder of the project (alternatively write path)
WDir <- "Downloads/Fermentation_BSFL"

## Import and prepare data ##

# Rearing data
Data <- read.xlsx(file.path(WDir, "DataExp.xlsx"), sheetName = "Experiment", stringsAsFactors = TRUE)
Eliminate <- Data$sample[Data$day==8&Data$survival<0.51] # eliminate replicates with too few larva at the beginning of the Exp
Data <- subset(Data, sample!=Eliminate)
dataDay0 <- subset(Data, !duplicated(sample)) # create new rows for the beginning of the experiment, to use in plots
dataDay0$day <- 0
dataDay0$avg_mass_larva_mg <- dataDay0$larvae_start 
Data <- rbind(Data,dataDay0) # merge data together
DataH <- subset(Data, day==8) # only data from Harvest day

# HPLC data
HPLC <- read.xlsx(file.path(WDir, "HPLC.xlsx"), sheetIndex=1, stringsAsFactors=TRUE)
HPLC$int <- interaction(HPLC$treat,HPLC$preserv)
HPLC$ID <- row.names(HPLC)

# Microbial data
FileList <- read.xlsx(file.path(WDir, "NamingMicrobialCommunity.xlsx"), sheetIndex=1, stringsAsFactors=TRUE)

Files_feed <- subset(FileList, type=="feed") # Pick only files with feedstock data 
Files_feed <- Files_feed[with(Files_feed, order(day, preserv)),]

for (file in 1:nrow(Files_feed)){ #Build a tree and merge OTUs tables
  name.file = file.path(WDir, "Raw_Data", as.character(Files_feed$file)[file])
  
  d <- read.delim(name.file,sep="\t")
  t <- d[, c("tax_id","lineage")]
  f <- d[, 1:2]
  ec <- d[, c(1,4)]
  colnames(f)[2] <- as.character(Files_feed$name[file])
  colnames(ec)[2] <- as.character(Files_feed$name[file])
  
  
  t$lineage = substr(t$lineage, 1, nchar(t$lineage)-1)
  t2<-t |>
    separate_wider_delim(lineage, delim = ";", names = c("Kingdom", "Phylum","Class","Order","Family","Genus","Species"),too_few="align_start",too_many="merge")
  
  if (file==1){df = f
  tree = t2
  ECC = ec} else {
    df = merge(df, f, by="tax_id", no.dups=TRUE, all=TRUE)
    tree = rbind(tree, subset(t2, !(tax_id %in% tree$tax_id)))
    ECC = merge(ECC, ec, by="tax_id", no.dups=TRUE, all=TRUE)
  }
}

colnames(df)[1] <- "OTU"
colnames(tree)[1] <- "OTU"
colnames(ECC)[1] <- "OTU"

ECC <- ECC %>% 
  mutate_if(is.numeric, round)

SumCounts <- rowSums(ECC[,2:9], na.rm = TRUE)

ECC <- ECC[SumCounts>10,] # For next steps, only select taxas with at least 10 counts
df <- df[SumCounts>10,]

Microb <- merge(df,tree,by=c("OTU")) 

# Temperature data

Temp <- read.xlsx(file.path(WDir, "TempSubstrate.xlsx"), sheetName = "Sheet1", stringsAsFactors = TRUE)
Tcontainer <- read.xlsx(file.path(WDir, "TempContainer.xlsx"), sheetName = "Sheet1", stringsAsFactors = TRUE)
Start <- read.xlsx(file.path(WDir, "StartExp.xlsx"), sheetName="Sheet1", stringsAsFactors =TRUE)
Tcontainer$dash = "Dashed"

Tlong <- Temp %>% 
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
 
## Fig.1 Larval growth and temperatures

# Fig.1ACE

WWmax = max(Data$avg_mass_larva_mg) #max larval weight

F1A <- ggplot(data=Data[Data$treat == "1D",]) + 
  geom_jitter(aes(x = day,y = avg_mass_larva_mg, fill = preserv, shape = preserv), width = 0.15, size = 1.5, alpha = 0.7, stroke = 1) +
  ylim(0,WWmax) + stat_summary(aes(day, avg_mass_larva_mg,col = preserv, linetype = preserv), fun = mean,geom = "line", linewidth = 1.5, show.legend=FALSE, alpha=0.6) +
  stat_summary(aes(day,avg_mass_larva_mg, fill=preserv, shape=preserv), fun=mean, geom="point",size=2,stroke=.8,show.legend=FALSE) +
  labs(y = "Avg larval weight (mg)", x = "Rearing duration (days)") + ggtitle("A   1-day storage") +
  scale_color_manual(values = c("ACD" = AC,"CTR" = WF,"FRM" = IF)) +
  scale_shape_manual(name = "Treatment", values = c("ACD" = sAC, "CTR" = sWF, "FRM" = sIF)) +
  scale_fill_manual(name = "Treatment", values = c("ACD" = AC, "CTR" = WF, "FRM" = IF)) +
  scale_linetype_manual(name = "Treatment", values=c("ACD" = "dashed","CTR" = "twodash", "FRM" = "dotdash")) +
  scale_x_continuous(breaks = seq(0,8,2)) + geom_rangeframe() + theme_classic() + theme(legend.position = 'none') 

F1C <- ggplot(data=Data[Data$treat == "1W",]) + 
  geom_jitter(aes(x = day,y = avg_mass_larva_mg, fill = preserv, shape = preserv), width = 0.15, size = 1.5, alpha = 0.7, stroke = 1) +
  ylim(0,WWmax) + stat_summary(aes(day, avg_mass_larva_mg,col = preserv, linetype = preserv), fun = mean, geom = "line", linewidth = 1.5, show.legend=FALSE, alpha=0.6) +
  stat_summary(aes(day,avg_mass_larva_mg,fill=preserv,shape=preserv), fun = mean, geom = "point", size = 2, stroke = .8, show.legend=FALSE)+
  labs(y = "Avg larval weight (mg)", x = "Rearing duration (days)") + ggtitle("C   7-day storage") +
  scale_color_manual(values = c("ACD" = AC,"CTR" = WF,"FRM" = IF)) +
  scale_shape_manual(name = "Treatment", values = c("ACD" = sAC, "CTR" = sWF, "FRM" = sIF)) +
  scale_fill_manual(name = "Treatment", values = c("ACD" = AC, "CTR" = WF, "FRM" = IF)) +
  scale_linetype_manual(name = "Treatment", values=c("ACD" = "dashed","CTR" = "twodash", "FRM" = "dotdash")) +
  scale_x_continuous(breaks = seq(0,8,2)) + theme_classic() + theme(legend.position = 'none') + geom_rangeframe() 

F1E <- ggplot(data=Data[Data$treat == "2W",]) + 
  geom_jitter(aes(x = day,y = avg_mass_larva_mg, fill = preserv, shape = preserv), width = 0.15, size = 1.5, alpha = 0.7, stroke = 1) +
  ylim(0,WWmax) + stat_summary(aes(day, avg_mass_larva_mg,col = preserv, linetype = preserv), fun = mean, geom = "line", linewidth = 1.5, show.legend=FALSE, alpha=0.6) +
  stat_summary(aes(day,avg_mass_larva_mg,fill=preserv,shape = preserv), fun = mean, geom="point",size=2,stroke=.8)+
  labs(y = "Avg larval weight (mg)", x = "Rearing duration (days)") + ggtitle("E   14-day storage") +
  scale_color_manual(values = c("ACD" = AC,"CTR" = WF,"FRM" = IF), labels = c("AC","WF","IF")) +
  scale_shape_manual(name = "Treatment", values = c("ACD" = sAC, "CTR" = sWF, "FRM" = sIF), labels = c("AC", "WF", "IF")) +
  scale_fill_manual(name = "Treatment", values = c("ACD" = AC, "CTR" = WF, "FRM" = IF), labels = c("AC", "WF", "IF")) +
  scale_linetype_manual(name = "Treatment", values=c("ACD" = "dashed","CTR" = "twodash", "FRM" = "dotdash")) +
  scale_x_continuous(breaks = seq(0,8,2)) + theme_classic() + geom_rangeframe()

legend <- get_legend(F1E + theme_minimal(base_size = 10) + theme(legend.position="bottom")) 

theme_set(theme_cowplot(font_size = 12) + 
            theme(text = element_text(size = 12)))

F1ACE <- cowplot::plot_grid(F1A, F1C, F1E + theme(legend.position = 'none'), legend, ncol = 1, rel_heights = c(1,1,1, .1))
file_name = paste(WDir,"/Plots/Fig1_ACE.tiff", sep = "")
tiff(file_name, units = "mm", width = 91, height = 200, res = 400)
print(F1ACE)
dev.off()

# Fig.1BDF
F1B <- ggplot() + geom_line(data = Temp1D, aes(Time, Temp, group = Treat, colour = Treat), linewidth = 1.2, show.legend = FALSE, alpha= 0.7) +
  geom_line(data = Tcontainer[Tcontainer$Date>=Start$start[Start$Exp == "1D"] & Tcontainer$Date<=Start$end[Start$Exp == "1D"],], 
            aes(Date, Temperature, col = dash), linewidth = 1.2, show.legend = F, alpha= 0.7) + theme_classic() +
  scale_color_manual(values = c("ACD" = AC,"CTR" = WF,"FRM" = IF,"Dashed" = "black")) +
  labs(y = "Temperature (°C)",x = "Rearing duration (days)") + ggtitle("B   1-day storage") + theme(legend.position="none") + ylim(20,50) + 
  scale_x_datetime(breaks = as.POSIXct(c("2023-06-06","2023-06-08","2023-06-10","2023-06-12", "2023-06-14")), labels = c("0","2","4","6","8 "))

F1D <- ggplot() + geom_line(data = Temp1W, aes(Time, Temp, group = Treat, colour = Treat), linewidth = 1.2, show.legend = FALSE, alpha= 0.7) +
  geom_line(data = Tcontainer[Tcontainer$Date>=Start$start[Start$Exp == "1W"] & Tcontainer$Date<=Start$end[Start$Exp == "1W"],], 
            aes(Date, Temperature, col = dash), linewidth = 1.2, show.legend = F, alpha= 0.7) + theme_classic() +
  scale_color_manual(values = c("ACD" = AC,"CTR" = WF,"FRM" = IF,"Dashed" = "black")) +
  labs(y = "Temperature (°C)",x = "Rearing duration (days)") + ggtitle("D   7-day storage") + theme(legend.position="none") + ylim(20,50) + 
  scale_x_datetime(breaks = as.POSIXct(c("2023-06-19","2023-06-21","2023-06-23","2023-06-25", "2023-06-27")),labels = c("0","2","4","6","8 "))

F1F <- ggplot() + geom_line(data = Temp2W, aes(Time, Temp, group = Treat, colour = Treat), linewidth = 1.2, alpha= 0.7) +
  geom_line(data = Tcontainer[Tcontainer$Date>=Start$start[Start$Exp == "2W"] & Tcontainer$Date<=Start$end[Start$Exp == "2W"],], 
            aes(Date, Temperature, col = dash), linewidth = 1.2, show.legend = F, alpha= 0.7) + theme_classic() +
  scale_color_manual(values = c("ACD" = AC,"CTR" = WF,"FRM" = IF,"Dashed" = "black"), labels = c("ACD" = "AC", "CTR" = "WF", "FRM" = "IF", "Dashed" = "Ambient\ntemperature")) +
  labs(y = "Temperature (°C)",x = "Rearing duration (days)") + ggtitle("F   14-day storage") + ylim(20,50) + 
  scale_x_datetime(breaks = as.POSIXct(c("2023-06-13","2023-06-15","2023-06-17","2023-06-19", "2023-06-21")),labels = c("0","2","4","6","8 "))

legend <- get_legend(F1F + theme_minimal(base_size = 10) + theme(legend.position="bottom"))

theme_set(theme_cowplot(font_size = 12) + 
            theme(text = element_text(size = 12)))

F1BDF <- cowplot::plot_grid(F1B, F1D, F1F + theme(legend.position = 'none'), legend, ncol = 1, rel_heights = c(1,1,1, .1))
file_name = paste(WDir,"/Plots/Fig1_BDF.tiff", sep = "")
tiff(file_name, units = "mm", width = 91, height = 200, res = 400)
print(F1BDF)
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
                 annotation_colors = list(Treatment = c(Pre = PRE,AC = AC,IF = IF,WF = WF),Duration=c("1D" = D1,"7D" = D7,"14D" = D14)),
                 cluster_cols = FALSE,cluster_rows = FALSE)


file_name = paste(WDir,"/Plots/Fig_2.tiff", sep= "")
tiff(file_name, units = "mm", width = 182, height = 91, res = 300)
print(Fig2)
dev.off()

# Figure 3 - Barplots of microbial communities

MC <- amp_load(otutable = df,taxonomy = tree, metadata = Files_feed)
MC_counts<- amp_load(otutable = ECC, taxonomy = tree, metadata = Files_feed)

MC_Filt = amp_filter_taxa(MC, tax_vector = c("Cyanobacteria","Mitochondria"), normalise = T, remove = TRUE)
MC_counts_F=amp_filter_taxa(MC_counts,tax_vector = c("Cyanobacteria", "Mitochondria"), normalise = F, remove = TRUE)


MC_Filt$metadata$name <- factor(MC_Filt$metadata$name, levels=c("AC1D", "IF1D", "WF1D", "AC7D", "IF7D", "WF7D", "AC14D", "IF14D", "WF14D"))

MC_long <- amp_export_long(MC_Filt)

MC_long[MC_long$Phylum == ""] <- NA
phylums <- MC_long %>% drop_na(Phylum) %>%  # get the 5 most abundant Phylum
  group_by(Phylum) %>% 
  summarise(count = sum(count)) %>% 
  slice_max(count, n = 5) %>% pull(Phylum)
MC_phylum <- MC_long[MC_long$Phylum %in% phylums,]

MC_phylum$name <- factor(MC_phylum$name, levels = c( "AC1D", "IF1D", "WF1D", "AC7D",  "IF7D",  "WF7D", "AC14D", "IF14D",  "WF14D"))


F3A <- ggplot(MC_phylum, aes(name, count, fill = factor(Phylum))) + 
  geom_bar(stat="identity", position="stack") + ggtitle("A") + 
  labs(y = "Relative Abundance (%)", x = "", fill = "Phylum") + theme_classic() +
  theme(legend.key = element_blank(),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 45, vjust = 1, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") + scale_fill_manual(values = micro[1:5])

file_name = paste(WDir, "/Plots/Fig3_A.tiff", sep = "")
tiff(file_name, units = "mm", width = 182, height = 100, res = 400)
print(F3A)
dev.off()

MC_long[MC_long$Genus == ""] <- NA
genuses <- MC_long %>% drop_na(Genus) %>%  # get the 12 most abundant Genus
  group_by(Genus) %>% 
  summarise(count = sum(count)) %>% 
  slice_max(count, n = 12) %>% pull(Genus)


MC_genus <- MC_long [ MC_long$Genus %in% genuses, ]

MC_genus$name <- factor(MC_genus$name, levels = c( "AC1D", "IF1D", "WF1D", "AC7D",  "IF7D",  "WF7D", "AC14D", "IF14D",  "WF14D"))


F3B <- ggplot(MC_genus, aes(name, count, fill = factor(Genus))) + 
  geom_bar(stat="identity", position="stack") + ggtitle("B") + 
  labs(y = "Relative Abundance (%)", x = "", fill = "Genus") + theme_classic() +
  theme(legend.key = element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 45, vjust = 1, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") + scale_fill_manual(values = micro)

file_name = paste(WDir, "/Plots/Fig3_B.tiff", sep = "")
tiff(file_name, units = "mm", width = 182, height = 100, res = 400)
print(F3B)
dev.off()

# Microbial diversity

Microb_div <- amp_alpha_diversity(MC_counts_F, measure = c("observed","shannon","simpson"), richness = TRUE, rarefy = NULL)
write.xlsx2(Microb_div, file.path(WDir, "MicrobialDiversity.xlsx") , sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE)

# Figure 4 - DPCA Metabolites

metab <- row.names(pHPLC) #Use the same 6 metabolites that are in the heatmap
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
  geom_point(size = 4, stroke = 2) + geom_mark_ellipse(aes(x = LD1, y = LD2, fill = preserv, group = int), con.cap = 0, show.legend = FALSE)+
  labs(x = "DPC1 (57.5%)", y = "DPC2 (29.6%)") + theme_classic(base_size = 20) + ylim(minY, maxY) + xlim(minX, maxX) +
  theme(axis.text = element_text(face = "bold"), legend.title = element_blank()) +
  scale_fill_manual(values = c("ACD" = AC, "CTR" = WF, "FRM" = IF, "PRE" = PRE), labels = c("ACD" = "AC", "CTR" = "WF", "FRM" = "IF", "PRE" = "Pre"))+
  scale_shape_manual(values = c("ACD" = sAC, "CTR" = sWF, "FRM" = sIF, "PRE" = sPR), labels=c("ACD" = "AC", "CTR" = "WF", "FRM" = "IF", "PRE" = "Pre"))

file_name = paste(WDir,"/Plots/HPCL.pdf", sep="") #Unfortunately, best way I found to make labels legible was to add them with a design software 
pdf(file_name, width=8, height=4.5)
print(DPCA)
dev.off()

## Statistical results as reported in paper ##

## Larvae NC
#N&C&F
NC <- read.xlsx(paste(WDir,"NCLarvae.xlsx", sep=""),sheetIndex = 1) #Load data on larvae nitrogen and carbon content
NC$int <- interaction(NC$treat,NC$preserv)
NC$treat <- as.factor(NC$treat)

NCavg<-NC %>% drop_na(number) %>%
  group_by(treat,preserv,sample) %>%
  summarise(avg_N = mean(N, na.rm = TRUE), sd_N = sd(N, na.rm = TRUE), avg_C = mean(C, na.rm = TRUE), sd_C = sd(C, na.rm = TRUE),
            avg_DM = mean(DM, na.rm = TRUE), sd_DM = sd(DM, na.rm = TRUE), avg_A = mean(AshDry, na.rm = TRUE),sd_A = sd(AshDry, na.rm = TRUE))

Fat<- read.xlsx(paste(WDir,"LarvaeFat.xlsx", sep=""),sheetIndex = 2) #Load data on larvae fat content
Fat$int<- interaction(Fat$treat,Fat$preserv)
Fat$treat<-as.factor(Fat$treat)       


NCF<-merge(NCavg,Fat,by=c("sample","preserv","treat"))   
NCFH<-merge(NCF,DataH,by=c("sample","preserv","treat"))

shapiro.test(DataH$avg_dry_larva) #Compare larval dry weight
hist(DataH$avg_dry_larva)
leveneTest(DataH$avg_dry_larva, group = DataH$preserv)

LDm_1D<-aov(larvae_dry ~ preserv, data=DataH[DataH$treat=="1D",]) #Compare larval dry weight for 1 day storage
summary(LDm_1D)
plot(LDm_1D)
TukeyHSD(LDm_1D)

LDm_1W<-aov(larvae_dry ~ preserv, data=DataH[DataH$treat=="1W",]) #Compare larval dry weight for 1 week storage
summary(LDm_1W)
TukeyHSD(LDm_1W)

LDm_2W<-aov(larvae_dry ~ preserv, data=DataH[DataH$treat=="2W",]) #Compare larval dry weight for 2 weeks storage
summary(LDm_2W)
TukeyHSD(LDm_2W)

LDm<-aov(avg_dry_larva ~ preserv*treat, data=DataH) #Dry weight
summary(LDm)
emmeans(LDm,specs=pairwise ~ preserv|treat)

CRm<-aov(dry_ConversionRate ~ preserv*treat, data=DataH) #Conversion rate
summary(CRm)
emmeans(CRm,specs=pairwise ~ preserv|treat)

WRm<-aov(dry_WasteReduction ~ preserv*treat, data=DataH) #Waste reduction
summary(WRm)
emmeans(WRm,specs=pairwise ~ preserv|treat)

DataH$survival[DataH$survival > 1] <- 1  #Survival
Sm<-aov(survival ~ preserv*treat, data=DataH)
summary(Sm)
emmeans(Sm,specs=pairwise ~ preserv|treat)

PPm<-aov(PP ~ preserv*treat, data=DataH) #prePupae
summary(PPm)
emmeans(PPm,specs=pairwise ~ preserv|treat)

Ashm<-aov(avg_A ~ preserv*treat, data=NCFH) #Ash content
summary(Ashm)
emmeans(Ashm,specs=pairwise ~ preserv|treat)

Fatm<-aov(fat ~ preserv*treat, data=NCFH) #Fat content
summary(Fatm)
emmeans(Fatm,specs=pairwise ~ preserv|treat)

Protm<-aov(avg_N * NcP ~ preserv*treat, data=NCFH) #Protein content
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
