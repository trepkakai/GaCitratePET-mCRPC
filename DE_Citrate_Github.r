# Imports
library(EnhancedVolcano) # for volcano plot
library(magrittr)
library(edgeR) # for GEX analysis

## Set working directory to filepath with relevant sequencing TPM file  
## and Design Matrix that will be used in EnhancedVolcano, created in Python
# setwd(filepath)

## STANDARDS
fccut = 1
alpha = 0.1

## BONE ONLY

# Import the gene expression and design matrices
tpm <- read.csv("tpm_citrate_bone.csv", header = TRUE, row.names="gene")
des <- read.csv("des_citrate_bone.csv", header = TRUE, row.names="samples")

# Analysis - low vs high GaCitrate uptake. High GaCitrate uptake appears on the right. 
y <- DGEList(counts=tpm, group=des$ga_high)
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)

# Change names for compatibility with Enhanced Volcano, including adding FDR
et <- data.frame(et)
names(et)[names(et) == "logFC"] <- "log2FoldChange"
et$PValue = p.adjust(et$PValue, method="fdr")
names(et)[names(et) == "PValue"] <- "FDR"

# Plot
pdf("GaCitrate_bone.pdf", width = 10, height = 10) # Initialize the export. 
EnhancedVolcano(data.frame(et),
                lab = rownames(et),
                x = 'log2FoldChange',
                y = 'FDR',
                ylim = c(0,2.5), # Set ylimit
                xlab = 'log2FC',
                ylab = '-log10(FDR)',
                title = 'Low vs high 68Ga-citrate uptake in bone metastases',
                legendLabels=c('Not significant','abs(log2FC) > 1','FDR < 0.1',
                               'FDR < 0.1 and abs(log2FC) > 1'),
                pCutoff = alpha,
                cutoffLineCol = 'grey',
                FCcutoff = fccut,
                pointSize = 3.0,
                labSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                #colAlpha = 1,
                lengthConnectors = unit(0, "npc"),
                arrowheads = FALSE,
                colConnectors = 'black')
dev.off() # Close the export 

# Export the significant results
filt = et[et$FDR < alpha,]
colnames(filt) <- c('log2FoldChange', 'logTPM', 'FDR')
filt = filt[filt$log2FoldChange < -fccut | filt$log2FoldChange > fccut,]
write.csv(filt,"jupyter/GaCitrate_bone.csv", row.names = TRUE)




