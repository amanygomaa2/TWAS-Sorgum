#!/usr/bin/env Rscript

setwd("/work/schnablelab/vladimir/HarshitaTWAS/FarmCPU")
rm(list = ls())

chooseCRANmirror(ind=72)
source("http://zzlab.net/GAPIT/gapit_functions.txt") #install this directly in R and not here, it will give error
#source("/work/schnablelab/vladimir/GAPIT_functions/gapit_functions.txt")

install.packages("data.table")
install.packages("tidyverse")
library("data.table")
library("tidyverse")

analysis <- "TWAS_farmCPU_RMIP"
out <- paste0("", analysis)

# Ensure 'TWASres' folder exists
if (!dir.exists(out)) {
  dir.create(out)
}

dateDIR <- as.character(Sys.Date())
pathout <- file.path(out, dateDIR)

if (!dir.exists(pathout)) {
  dir.create(pathout)
}




#load the new phe data
phe <- fread("phenotypesmatched.csv",data.table = F) #pheno_693.txt 
phe <- phe[,-1]
colnames(phe)
trait=which(colnames(phe) == "blues") #args[1] #blues
colnames(phe[trait])

colnames(phe)[1] <- "taxa"

###calculate means by taxa from phe
# phe.NEmedians <- phe %>%
#   group_by(taxa) %>%
#   dplyr::select(taxa, colnames(phe[trait])) %>%
#   summarise(across(everything(), ~mean(., na.rm = TRUE)))
# 
# phe.NEmedians <- data.frame(phe.NEmedians)


#load counts data
counts <- fread("tpmfinal.csv", data.table = F)
counts <- counts[,-1]
colnames(counts)[1:5]
colnames(counts)[1] <- "taxa"

row.names(counts) <- counts$taxa

# Exclude the first column and convert to numeric
counts_numeric <- counts[,-1]
counts_numeric <- lapply(counts_numeric, as.numeric)
counts_numeric <- as.data.frame(counts_numeric)

keep <- colSums(counts_numeric>.1) >= nrow(counts_numeric)/2
counts_combined.filtered <- counts_numeric[,keep]
print(paste0("Genes keeped: ", ncol(counts_combined.filtered)))


# Combine with the first column
counts_combined <- cbind(counts[,1], counts_combined.filtered)
colnames(counts_combined)[1] <- "taxa"
row.names(counts_combined) <- counts_combined$taxa


prcompResult<-prcomp(counts_combined[,-1],center=TRUE,scale.=TRUE) #This should take less than a minute.
PCs<-prcompResult$x
PCsTop <- data.frame(PCs[,1:100])
PCs.2020 <- cbind(row.names(counts_combined),PCsTop)
colnames(PCs.2020)[1] <- "taxa"

plot(PCs.2020$PC1, PCs.2020$PC2)


NROW(merge(counts_combined[,1], phe, by = 1))


#use quantile method to handle outliers
Quantile<- apply(counts_combined[,-1],2,  # 2 indicates it is for column and 1 indicates it is for row
                 function(A){min_x=as.numeric(quantile(A,0.05));
                 max_x=as.numeric(quantile(A,0.95));
                 out<-(2*(A-min_x)/(max_x-min_x));
                 out[out>2]<-2;out[out< 0]<- 0;return(out)})

#Quantile.t <-  as.data.frame(t(Quantile))
Quantile.t <- as.data.frame(Quantile)
Quantile.t$taxa <- row.names(Quantile.t)
myGD <-  Quantile.t[,c(ncol(Quantile.t),1: (ncol(Quantile.t)-1))]


myY <- phe
colnames(myY)[1] <- "taxa"

myGM <- fread("geneposition.csv", data.table = F)
unique(myGM$chr) #only cromosomes 
myGM <- myGM[,-1]
myGM <- myGM[myGM$GeneID %in% colnames(counts_combined[,-1]),]

# #covariates
# covariates <- read.csv("sampling_693.order.csv", head = TRUE)
# #colnames(covariates) #DaysToPollen
# myCV <- covariates
# colnames(myCV)
# myCV2 <- myCV[,-c(2)]
# colnames(myCV2)
# 
# ###calculate means by taxa from anthesis as covariate
# cov.ft <- phe %>%
#   group_by(taxa) %>%
#   dplyr::select(taxa, 'Anthesis.sp') %>%
#   summarise(across(everything(), ~mean(., na.rm = TRUE)))
# 
# cov.ft <- data.frame(cov.ft)
# 
# cvF <- merge(myCV2, cov.ft, by ='taxa')
# cvF <- merge(cvF, PCs.2020[,c(1:4)], by ='taxa')



#make a conditional statement to remove 'anthesis.sp' as covariate in the trait is a flowering time trait
# if (colnames(phe[trait]) == 'Anthesis.sp') {
#   cvF <- cvF[,-3]
#   print("flowering time trait")
# } else if (colnames(phe[trait]) == 'Silking.sp') {
#   cvF <- cvF[,-3]
#   print("flowering time trait")
# } else {
#   cvF <- cvF
#   print("including flowering time as covariate")
# }
# 
# colnames(cvF)

#create empty data to run resampling
gen <- data.frame()
sampled_rows_list <- list()
install.packages("GAPIT")
library(GAPIT)


for (i in 1:100) {
  # Sample 90% of the rows
  sampled_rows <- sort(sample(x = length(phe$taxa), size = round(0.9 * length(phe$taxa))))
  sampled_rows_list[[paste0("sampled_rows", i)]] <- sampled_rows
  
  # Extract the sampled rows
  sampled_data <- phe$taxa[sampled_rows_list[[paste0("sampled_rows", i)]]]
  
  
  myGAPIT <- GAPIT(Y=myY[myY$taxa %in% sampled_data,],
                   GD=myGD,
                   GM=myGM,
                   #CV=cvF,
                   PCA.total=3,
                   model= "FarmCPU",
                   SNP.MAF=0,
                   file.output=F
  )
  #warnings()
  
  #getting the important genes and Manhattan plots
  values <- data.frame(myGAPIT$GWAS)
  values$log <- -log10(values$P.value)
  threshold <- -log10(0.05/round(0.9 * length(phe$taxa)))
  
  values2 <- values[which(values$log > threshold),]
  
  gen <- rbind(gen, values2[c(1,4,8)])
}


results <- data.frame(table(gen$SNP))
results2 <- merge(myGM, results, by = 1 ,all.x = T)

fwrite(results2, paste0(pathout,"/TWAS.farmCPU_",colnames(phe[trait]),".csv"))


###line to stop 
stop("Stopping the script here")

####

pathout <- "/work/schnablelab/vladimir/HarshitaTWAS/FarmCPU/TWAS_farmCPU_RMIP/2024-05-27"
gwas.datTWAS<- fread(paste0(pathout,"/TWAS.farmCPU_",colnames(phe[trait]),".csv"), data.table = F)

theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
             plot.title = element_text(hjust = 0.5), plot.subtitle=element_text(hjust=0.5))

colnames(gwas.datTWAS)
colnames(gwas.datTWAS)[c(2:3)] <- c("chr", "start")


nCHR <- length(unique(gwas.datTWAS$chr))
gwas.datTWAS$BPcum <- NA
s <- 0
nbp <- c()

for (i in sort(unique(gwas.datTWAS$chr))){
  nbp[i] <- max(gwas.datTWAS[gwas.datTWAS$chr == i,]$start)
  gwas.datTWAS[gwas.datTWAS$chr == i,"BPcum"] <- gwas.datTWAS[gwas.datTWAS$chr == i,"start"] + s
  s <- s + nbp[i]
}

axis.set <- gwas.datTWAS %>% 
  group_by(chr) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

#IDs for genes with FDR less than 0.05
treshold <- .1 #.01 .1

hits <- gwas.datTWAS[which(gwas.datTWAS$Freq/100 > treshold),]

###merge with maize data
maize <- fread("B73_geneModels_v5_v2.csv", data.table = F)
colnames(maize)


ortholog <- fread("1_ortholog_list.csv", data.table = F)
colnames(ortholog)

ortholog2 <- merge(ortholog, maize , by.x = 2, by.y = 1 )
colnames(ortholog2)
ortholog2 <- ortholog2[,c(1,2,9,10,11,12,17:19)]


#to name hits in the plot 
genes.desc <- fread("Sbicolor_730_v5.1.annotation_info.txt", data.table = F)
colnames(genes.desc)

genes.desc <- genes.desc[!duplicated(genes.desc$locusName),]

genes.desc2 <- merge(genes.desc[,c(2,11,12)], ortholog2 , by.x = 1, by.y = 2, all.x = T)

setdiff(hits$GeneID, genes.desc2$locusName)
setdiff(hits$GeneID, genes.desc$locusName)


results3 <- merge(hits, genes.desc2, by = 1 )

length(unique(results3$GeneID))
results3$GeneID[duplicated(results3$GeneID)]

results3 <- results3[!duplicated(results3$GeneID),]

TWAS.mlm.hits.txt <- fread("TWAS.mlm.hits.txt", data.table = F)
colnames(TWAS.mlm.hits.txt)

results3 <- merge(results3, TWAS.mlm.hits.txt, by = 1, all.x = T)

fwrite(results3,"GenesFarmCPU.csv")

gwas.datTWAS <- merge(gwas.datTWAS, results3[,-c(2:5)], by = 1, all.x = T)


colnames(gwas.datTWAS)
nameonplot <- paste0("TWAS farmCPU ", colnames(phe[trait])," Sorghum 2021")
plot1 <- ggplot(data=gwas.datTWAS, aes(BPcum, Freq/100, label=name, colour=factor(chr, levels = c(1:10)))) + 
  ggrastr::rasterise(geom_point(size = 2), dpi = 900, dev = "ragg_png") + 
  geom_text( vjust="inward",hjust="inward", fontface="italic", size = 3) +
  geom_hline(yintercept = 0.1 , linetype=2) + #-log10(0.05/nrow(gwas.datTWAS)) fdr.treshold
  scale_color_manual(values = c('#03193F','#28708C','#BF930F','#0f3bbf','#295E52','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')) + 
  annotate("text", label= nameonplot, y=1, x=100, size=5,hjust="inward") +
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center) + 
  theme(legend.position = "none") + 
  ylab("RMIP") + 
  xlab("Chromosome") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, .1))
#+ labs(title="", subtitle=paste0(trait))

plot1


folder <- "ManhattanPlots"
out2 <- paste0("", pathout,"/",folder)

# Ensure 'TWASres' folder exists
if (!dir.exists(out2)) {
  dir.create(out2)
}

#ggsave(paste0(pathout,"TWAS.res_farmCPU_",colnames(phe[trait]),".png"), height = 5, width = 10)
ggsave(paste0(out2,"/TWAS_farmCPU_",colnames(phe[trait]),".pdf"), height = 5, width = 10)

# Correlation Matrix TWAS mlm hits 
TWAS.mlm <- counts_combined[,c("taxa",intersect(TWAS.mlm.hits.txt$GeneID, colnames(counts_combined)))]

# Find the intersection of GeneIDs with column names of counts_combined
matching_genes <- intersect(TWAS.mlm.hits.txt$GeneID, colnames(counts_combined))

# Subset TWAS.mlm.hits.txt to include only rows with matching GeneIDs
matching_rows <- TWAS.mlm.hits.txt[TWAS.mlm.hits.txt$GeneID %in% matching_genes, ]

colnames(TWAS.mlm) <- c("taxa", matching_rows$name)
TWAS.mlm <- TWAS.mlm[,-1]

library(corrplot)
M<-cor(TWAS.mlm)

# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(TWAS.mlm)
head(p.mat[, 1:5])


col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         #p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)


genes <- counts_combined[,c("taxa","Sobic.001G121400", "Sobic.003G017200", "Sobic.004G003434", "Sobic.002G312300", "Sobic.001G227900", "Sobic.010G085400")]
colnames(genes) <- c("taxa", "Sobic.001G121400", "pebp14","mads76","sbp29", "glk45", "mads74")

genes2 <- merge(genes, phe, by = "taxa")
colnames(genes2)

theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
             plot.title = element_text(hjust = 0.5), plot.subtitle=element_text(hjust=0.5))

Sobic.001G121400 <- ggplot() +
  rasterise(geom_point(data= genes2, aes(x=blues, y=Sobic.001G121400, fill="red4"), size = 3, shape = 21),dpi=600, dev = "ragg_png") +
  labs(x = "flowering time",
       y = expression(italic("Sobic.001G121400")),
       fill = "Group") +
  theme(legend.position = "none", strip.background = element_blank()) +
  #facet_grid(cols = vars(Group))+  
  #scale_x_continuous(breaks = seq(1000, 1900, 200), labels = seq(1000, 1900, 200),limits = c(1000, 1900)) +
  stat_poly_line(data= genes2, aes(x=blues, y=Sobic.001G121400), method=lm,  linetype="dashed", se=FALSE, linewidth = 1,
                 color="black", fill="blue") +
  stat_poly_eq(data= genes2, aes(x=blues, y=Sobic.001G121400),size = 5)
Sobic.001G121400



pebp14 <- ggplot() +
  rasterise(geom_point(data= genes2, aes(x=blues, y=pebp14, fill="red4"), size = 3, shape = 21),dpi=600, dev = "ragg_png") +
  labs(x = "flowering time",
       y = expression(italic("pebp14")),
       fill = "Group") +
  theme(legend.position = "none", strip.background = element_blank()) +
  #facet_grid(cols = vars(Group))+  
  #scale_x_continuous(breaks = seq(1000, 1900, 200), labels = seq(1000, 1900, 200),limits = c(1000, 1900)) +
  stat_poly_line(data= genes2, aes(x=blues, y=pebp14), method=lm,  linetype="dashed", se=FALSE, linewidth = 1,
                 color="black", fill="blue") +
  stat_poly_eq(data= genes2, aes(x=blues, y=pebp14),size = 5)
pebp14


sbp29 <- ggplot() +
  rasterise(geom_point(data= genes2, aes(x=blues, y=sbp29, fill="red4"), size = 3, shape = 21),dpi=600, dev = "ragg_png") +
  labs(x = "flowering time",
       y = expression(italic("sbp29")),
       fill = "Group") +
  theme(legend.position = "none", strip.background = element_blank()) +
  #facet_grid(cols = vars(Group))+  
  #scale_x_continuous(breaks = seq(1000, 1900, 200), labels = seq(1000, 1900, 200),limits = c(1000, 1900)) +
  stat_poly_line(data= genes2, aes(x=blues, y=sbp29), method=lm,  linetype="dashed", se=FALSE, linewidth = 1,
                 color="black", fill="blue") +
  stat_poly_eq(data= genes2, aes(x=blues, y=sbp29),size = 5)
sbp29


top_row <- plot_grid(Sobic.001G121400, pebp14, sbp29, nrow=1)
#top_row
ggsave(paste0(out2,"/TWAS_farmCPU_genes",colnames(phe[trait]),".pdf"), height = 3.8, width = 10)


#### scatter plot R2 v Freq

# Correlation Matrix TWAS mlm hits 
hits.exp <- counts_combined[,c("taxa",hits$GeneID)]
hits.exp2 <- merge(phe, hits.exp, by = 1)

### histograms of frequency
colnames(hits.exp)
correlation.TPMcu <- data.frame()

for (i in 2:ncol(hits.exp)) {
  #with raw TPM data
  x1 <- merge(phe, hits.exp[,c(1,i)], by = 1)
  r1 <- cor.test(x1[,2], x1[,3])$estimate^2
  
  ##ssave data
  dta <- data.frame(colnames(hits.exp)[i], r1)
  
  correlation.TPMcu <- rbind(correlation.TPMcu,dta)
}

colnames(correlation.TPMcu) <- c("taxa","squared.R")

his <- ggplot(correlation.TPMcu, aes(x = squared.R)) +
  geom_histogram( fill = "skyblue", color = "black") +
  labs( x = expression("R"^2),
        y = "Frequency")
his 


sc <- merge(correlation.TPMcu, hits, by =1)
sc <- merge(sc, results3[,c(1,16)], by =1, all.x = T)

colnames(sc)
sc.plot <- ggplot() +
  rasterise(geom_point(data= sc, aes(x=Freq/100, y=squared.R), fill="grey", size = 4, shape = 21),dpi=600, dev = "ragg_png") +
  geom_text(data = sc, aes(x = Freq / 100, y = squared.R, label = name), size = 3, vjust = -1) + # Add labels to points
  labs(x = "RMIP",
       y = expression("R"^2),
       fill = "Group") +
  theme(legend.position = "none", strip.background = element_blank()) +
  #facet_grid(cols = vars(Group))+  
  #scale_x_continuous(breaks = seq(1000, 1900, 200), labels = seq(1000, 1900, 200),limits = c(1000, 1900)) +
  stat_poly_line(data= sc, aes(x=Freq/100, y=squared.R), method=lm,  linetype="dashed", se=FALSE, linewidth = 1,
                 color="black", fill="blue") +
  stat_poly_eq(data= sc, aes(x=Freq/100, y=squared.R),size = 5)
sc.plot

ggsave(paste0(out2,"/TWAS_farmCPU_genes_scatter.pdf"), height = 8, width = 8)

######
#correlation plot all genes
TWAS.farmcpu <- counts_combined[,hits$GeneID]


i=1
for (i in 1:nrow(results3)) {
  if(is.na(results3[i,16])){
    results3[i,16] <- results3[i,11]
  }
}

## some mistakes in the labels, they are not NA
results3[which(results3[,16] == ""),16] <- NA
results3[,16]


##fill thos NAs with gene ID
for (i in 1:nrow(results3)) {
  if(is.na(results3[i,16])){
    results3[i,16] <- results3[i,1]
  }
}


##name column in hits

for (i in 2:ncol(TWAS.farmcpu)) {
  colnames(TWAS.farmcpu)[i][which(colnames(TWAS.farmcpu)[i] == results3[i,1])] <- results3[i,16]
}



library(corrplot)
Mf<-cor(TWAS.farmcpu)

# # mat : is a matrix of data
# # ... : further arguments to pass to the native R cor.test function
# cor.mtest <- function(mat, ...) {
#   mat <- as.matrix(mat)
#   n <- ncol(mat)
#   p.mat<- matrix(NA, n, n)
#   diag(p.mat) <- 0
#   for (i in 1:(n - 1)) {
#     for (j in (i + 1):n) {
#       tmp <- cor.test(mat[, i], mat[, j], ...)
#       p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
#     }
#   }
#   colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
#   p.mat
# }
# matrix of the p-value of the correlation
p.mat.cpu <- cor.mtest(TWAS.farmcpu)
head(p.mat.cpu[, 1:5])


col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(Mf, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45,number.cex=0.5, #Text label color and rotation
         # Combine with significance
         #p.mat = p.mat.cpu, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)


fwrite(results3,paste0(pathout,"/GenesFarmCPU.csv"))


##### model james asked

#ft = gene expression + PC1 + PC2 + PC3
head(hits.exp2)

mldta <- merge(hits.exp2[,c("taxa","blues","Sobic.001G121400")], PCs.2020[,c("taxa","PC1","PC2","PC3","PC4","PC5")], by = "taxa")
lm = lm(blues~Sobic.001G121400+PC1+PC2+PC3+PC4+PC5, data = mldta) #Create the linear regression
summary(lm) #Review the results




# #####
# 
# counts2 <- counts[colnames(counts) %in% c("taxa",results3$B73_V5)]
# dta <- merge(myY ,counts2, by =1)
# 
# 
# grupos <- fread("/work/schnablelab/vladimir/NE2020_HISAT/2020_RNA/out.TWAS/metadata.csv", data.table = F)
# grupos$taxa[which(grupos$taxa == "INB_101LFY_LFY_(A632_X_M16_S5)")] <- "INB_101LFY/LFY_(A632_X_M16_S5)"
# 
# 
# dta.snps <- merge(grupos[,c(2,3)], dta, by = "taxa")
# dta.snps$Group[which(dta.snps$taxa == "INB_101LFY/LFY_(A632_X_M16_S5)")] <- "SS"
# 
# unique(dta.snps$Group)
# 
# ### get correlation persentage   
# #Function to determine if values in a vector are all positive, all negative, or contrasting
# check_values <- function(x) {
#   if (all(x > 0)) {
#     return("positive")
#   } else if (all(x < 0)) {
#     return("negative")
#   } else {
#     return("contrasting")
#   }
# }
# 
# correlation.farmCPU <- data.frame()
# i=4
# for (i in 4:ncol(dta.snps)) {
#   #with raw TPM data
#   x1 <- dta.snps[,c(2,3,i)]
#   r1 <- cor.test(x1[,2], x1[,3])$estimate
#   
#   x2 <- x1[which(x1$Group == "SS"),]
#   r2 <- cor.test(x2[,2], x2[,3])$estimate
#   
#   x3 <- x1[which(x1$Group == "NSS"),]
#   r3 <- cor.test(x3[,2], x3[,3])$estimate
#   
#   x4 <- x1[which(x1$Group == "IDT"),]
#   r4 <- cor.test(x4[,2], x4[,3])$estimate
#   
#   x5 <- x1[which(x1$Group == "Sweet corn"),]
#   r5 <- cor.test(x5[,2], x5[,3])$estimate
#   
#   x6 <- x1[which(x1$Group == "Popcorn"),]
#   r6 <- cor.test(x6[,2], x6[,3])$estimate
#   
#   x7 <- x1[which(x1$Group == "Tropical"),]
#   r7 <- cor.test(x7[,2], x7[,3])$estimate
#   
#   ##ssave data
#   dtaxx <- data.frame(colnames(dta.snps)[i], r1, r2 , r3, r4, r5, r6, r7)
#   
#   # Selecting a row from dtaxx for checking
#   row_to_check <- dtaxx[1, c(2:7)]
#   
#   # Applying the function to the selected row and printing the result
#   dtaxx$summary <- check_values(row_to_check)
#   
#   correlation.farmCPU <- rbind(correlation.farmCPU,dtaxx)
# }
# 
# colnames(correlation.farmCPU) <- c("gene","overall", "SS","NSS", "IDT","SweetCorn", "Popcorn", "Tropical", "summary")
# 
# #correlation.farmCPU[,c(2:6)] <- correlation.farmCPU[,c(2:6)]^2
# 
# correlation.farmCPU2 <- merge(hits, correlation.farmCPU, by = 1)
# correlation.farmCPU2.extend <- merge(results3, correlation.farmCPU, by = 1)
# 
# 
# fwrite(correlation.farmCPU2.extend,"res_farmCPU_RMIP_hits_silking.extended.csv")
# 
# 
# 
# ########
# groups <- c("SS","IDT", "NSS", "Sweet corn", "Popcorn", "Tropical", "Others")
# groups.col <- c("#D7DCB8", "#53413B", "#1E4261", "#C98D49", "#DBD68F", "darkblue", "white")
# 
# # Define the desired order of groups
# desired_order <- c("NSS", "Popcorn", "IDT", "SS", "Sweet corn", "Tropical", "Others")
# # Reorder the levels of the "Group" factor
# dta.snps$Group <- factor(dta.snps$Group, levels = desired_order)
# 
# dta.snps <- dta.snps[which(dta.snps$Group != "Others"),]
# 
# 
# order.freq <- results3$B73_V5[order(results3$Freq, decreasing = T)]
# 
# 
# # Load necessary libraries
# #library(ggplot2)
# library(gridExtra)
# 
# # Initialize a PDF file
# pdf("plots_res_farmCPU.pdf", width = 15, height = 10)  # Adjust width and height as needed
# 
# # Create an empty list to store plots
# plot_list <- list()
# 
# # Plotting loop
# for (i in order.freq) {
#   p <- ggplot() +
#     rasterise(geom_point(data = dta.snps, aes(x = colnames(phe[trait]), y = .data[[i]], fill = Group), size = 3, shape = 21), dpi = 600, dev = "ragg_png") +
#     labs(x = "days to silking",
#          y = paste0(i, " (", results3$Freq[which(results3$B73_V5 == i)], ")"),
#          fill = "Group",
#          subtitle = paste0(round(correlation.farmCPU2$overall[which(correlation.farmCPU2$B73_V5 == i)], 3), " ", results3$symbol[which(results3$B73_V5 == i)]," ", results3$`Gene Symbol`[which(results3$B73_V5 == i)], " ",
#                            results3$interproscan[which(results3$B73_V5 == i)])) +
#     theme(legend.position = "none", strip.background = element_blank(), plot.subtitle = element_text(size = 13)) +
#     scale_fill_manual(labels = groups, 
#                       values = groups.col, 
#                       breaks = groups,
#                       name = "") +
#     facet_grid(cols = vars(Group)) +  
#     #scale_x_continuous(breaks = c(65, 85), labels = c(65, 85), limits = c(65, 85)) +
#     stat_poly_line(data = dta.snps, aes(x = colnames(phe[trait]), y = .data[[i]]), method = lm,  linetype = "dashed", se = FALSE, linewidth = 1,
#                    color = "black", fill = "blue") +
#     stat_poly_eq(data = dta.snps, aes(x = colnames(phe[trait]), y = .data[[i]]), size = 5)
#   
#   plot_list[[i]] <- p
# }
# 
# # Function to arrange plots into three per page
# arrange_plots <- function(plot_list) {
#   for (i in seq(1, length(plot_list), 3)) {
#     plots_to_arrange <- plot_list[i:min(i + 2, length(plot_list))]
#     grid.arrange(grobs = plots_to_arrange, nrow = 3)
#   }
# }
# 
# # Arrange plots into three per page
# arrange_plots(plot_list)
# 
# # Close the PDF device
# dev.off()
