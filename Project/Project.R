setwd("/Users/sanch/Documents/WSU/Spring 2017/CSC 7300/Project/")

library(Biobase)
source("https://bioconductor.org/biocLite.R")

library(GEOquery)
library(limma)
library(ggplot2)
library(RColorBrewer) # For heatMap
library(randomGLM)
library(gtools)
#biocLite("e1071")
#library(xgboost)  
#library("archdata") 
library(caret)    
library(dplyr) 
library(e1071)
#biocLite("rpart")
library(rpart)
library(annotate)
library(hgu133a.db)


# load series and platform data from GEO

gset <- getGEO("GSE6613", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
var2 <- gset

SampleTable <- pData(gset)
group1 <- subset(SampleTable, SampleTable$characteristics_ch1 == "healthy control")
group2 <- subset(SampleTable, SampleTable$characteristics_ch1 == "Parkinson\'s disease")
group3 <- subset(SampleTable, SampleTable$characteristics_ch1 == "neurological disease control")

colmns <- colnames(gset)
gsms <- paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXX")

# healthy control versus Parkinson's disease
for(i in 1:nchar(gsms))
{
  if(colmns[i] %in% group1$geo_accession){
    gsms <- paste(substring(gsms, 1, i-1), substring(gsms, i+1), sep ="1")
  }
}
for(i in (1:105))
{
  if(colmns[i] %in% group2$geo_accession){
    gsms <- paste(substring(gsms, 1, i-1), substring(gsms, i+1), sep ="0")
  }
}

sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

var1 <- exprs(gset) #var1 would contain Healthy and PD samples


#Boxplot
ex2 <- var1[ , order(sml)]
sml2 <- sml[order(sml)]
fl2 <- as.factor(sml2)
labels <- c("PD","Healthy")

palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(ex2)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE6613", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex2, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl2)
legend("topleft", labels, fill=palette(), bty="n")


# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=22283)
write.table(tT, file=stdout(), row.names=F, sep="\t")

library(affy)
#MAS5 normalization
affy.data = ReadAffy()
save(affy.data, file = 'AffyData.Rdata')
eset.mas5 = mas5(affy.data)

exprSet.nologs = exprs(eset.mas5)
colnames(exprSet.nologs)
exprSet = log(exprSet.nologs, 2)
newVar1 <- exprSet
save(newVar1, file = 'newVar1.Rdata')

#HeatMap
library(ComplexHeatmap)

ImpG <- c("SERPINB9", "ST13", "SRP46", "PRPF4B", "PRPF4B", "RAP2B", "KIAA1128", "BUB3", "PURA", "SEC13L", "PRKBCB", "FLJ22222", "NAP1L1", "LRPPRC", "BCL2", "BCL11B", "SEPT6", "SH3BGRL", "SEPT6", "TCEA1", "UBE1B1", "TBPL1", "PEA15" )

Probes <- as.character(row.names(newVar1))
ProToGene <- select(hgu133a.db, Probes, c("SYMBOL", "GENENAME"))

ImpP <- ProToGene[which(ProToGene$SYMBOL %in% ImpG), 1]
ImpG2 <- ProToGene[which(ProToGene$SYMBOL %in% ImpG), 2]

smpls <- as.character(ImpG2)

dataHeatmap = newVar1[ImpP,  ]
rownames(dataHeatmap) = smpls
colnames(dataHeatmap) = sapply(strsplit(colnames(dataHeatmap), split='.', fixed=TRUE), function(x) (x[1]))
dataHeatmap = dataHeatmap[,order(sml)] 

Heatmap(scale(dataHeatmap), col = rev(brewer.pal(11, "RdBu")), column_title = "Healthy and PD samples", 
        column_title_side = "bottom", name="Expression")

#heatmap((newVar1[ImpP, ]), col = rev(brewer.pal(11, "RdBu")), name="Correlation", labRow = smpls)


#Volcano Plot, using ggplot 

tT$color <- "black"

tT[which( -log10(tT$P.Value)>2 & tT$logFC>1 ),"color"] <- "red"
tT[which( -log10(tT$P.Value)>2 & tT$logFC< -1 ),"color"] <- "green"
tT[which( -log10(tT$P.Value)>2 & tT$logFC> -1 & tT$logFC< 1),"color"] <- "grey"

L1 <- length(which(-log10(tT$P.Value) >2 & tT$logFC < -1))
L2 <- length(which(-log10(tT$P.Value) >2 & tT$logFC > 1))
L3 <- length(which(-log10(tT$P.Value) >2 & tT$logFC > -1 & tT$logFC < 1))
L4 <- length(which(-log10(tT$P.Value) <2 & tT$logFC > 1))
L5 <- length(which(-log10(tT$P.Value) <2 & tT$logFC < -1))
L6 <- length(which(-log10(tT$P.Value) <2 & tT$logFC > -1 & tT$logFC < 1))


volc = ggplot(tT, aes(logFC, -log10(P.Value))) +
  geom_point(aes(col=color)) +
  geom_vline(xintercept = c(-1,1), colour = "black") + #Add horizontal lines
  geom_hline(yintercept = 2, colour = "black") +  # Add vertical lines
  annotate(geom = "text", label = paste("N = ", L1), x = -2, y = 5, size = 5, colour = "black") + 
  annotate(geom = "text", label = paste("N = ", L2), x = 2, y = 5, size = 5, colour = "black") + 
  annotate(geom = "text", label = paste("N = ", L3), x = 0, y = 5, size = 5, colour = "black") +
  annotate(geom = "text", label = paste("N = ", L4), x = 2, y = 0, size = 5, colour = "black") + 
  annotate(geom = "text", label = paste("N = ", L5), x = -2, y = 0, size = 5, colour = "black") +
  annotate(geom = "text", label = paste("N = ", L6), x = -0.5, y = 0, size = 5, colour = "black") +
  scale_color_manual(values = c("red" = "red", 
                                "black" = "black", 
                                "green" = "green",
                                "grey" = "grey"))


#plot volcano
volc  


#Classification
randomPD = sample(1:50, 31)
randomGroup2 <- group2[randomPD,] 

randomHealthy = sample(1:22, 17)
randomGroup1 <- group1[randomHealthy,] 

randomND = sample(1:33, 18)
randomGroup3 <- group3[randomND,] 


mysamp1 <- c(as.character(randomGroup1$geo_accession),
            as.character(randomGroup2$geo_accession),
            as.character(randomGroup3$geo_accession)
           )
 

mysamp2 <- which(as.character(colnames(var2)) %in% mysamp1)
train = var2[,(mysamp2)]
test = var2[,-c(mysamp2)]

x = 1*(train$characteristics_ch1 == "Parkinson's disease")
z = 2*(train$characteristics_ch1 == "neurological disease control")

y = x + z;

length(which(y == 1) )
length(which(y == 2) )
length(which(y == 0) )


tT2 <- topTable(fit2, adjust="fdr", number=10000)

expression_mnz <- as.data.frame(t(exprs(var2)))
expression_mnz <- expression_mnz[,rownames(tT2)]
# add a column with the label of each sample
expression_mnz$label <- var2$characteristics_ch1

train = expression_mnz[sort(mysamp2),]
test = expression_mnz[-sort(mysamp2),]

#SVM with reduced features

svm.model <- svm(label ~ ., data = train,  scale = F) 
svm.pred  <- predict(svm.model, test[,-ncol(test)]) #Excluding the colunm label

## compute svm confusion matrix
table(pred = svm.pred, true = test$label)

## rpart
rpart.model <- rpart(label ~ ., data = train)
rpart.pred  <- predict(rpart.model, test[,-ncol(test)], type = "class")


## compute rpart confusion matrix
table(pred = rpart.pred, true = test$label)
