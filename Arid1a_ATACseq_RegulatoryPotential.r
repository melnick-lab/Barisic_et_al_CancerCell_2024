## TF regulatory potential from differential accessibility
## For "ARID1A Orchestrates SWI/SNF-mediated Sequential Binding of Transcription Factors and Its Deficiency Drives Pre-memory B-cell Fate and Lymphomagenesis"
## Code by Cem Meydan

library(motifmatchr)
library(JASPAR2018)
library(chromVAR)
library(chromVARmotifs)
library(gsubfn)
library(SummarizedExperiment)
library(broom)
library(tidyr)
library(lmtest)
library(sandwich)
library(modelr)
library(tidyverse)



## Function (borrowed from chromVAR/TFBStools) to access JASPAR motifs (here the CORE vertebrate collection):
jaspar = function (collection = "CORE", ...) 
{
  opts = list()
  opts["tax_group"] = "vertebrates"
  opts["collection"] = collection
  opts = c(opts, list(...))
  out = TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out)))) 
    names(out) = paste(names(out), TFBSTools::name(out), 
                        sep = "_")
  return(out)
}



filterPeaks = FALSE
GenomeRef = "mm10"

## Genome reference
if(GenomeRef == "mm10")
{
    library(BSgenome.Mmusculus.UCSC.mm10)
    Genome = BSgenome.Mmusculus.UCSC.mm10
}else if(GenomeRef == "hg38")
{
    library(BSgenome.Hsapiens.UCSC.hg38)
    Genome = BSgenome.Hsapiens.UCSC.hg38
}


## Get differential outputs for each comparison. Each file should have genomic coordinates and log2FoldChange column, and optionally counts for filtering
outPath = "."

diffFiles = list.files(outPath, ".*results_deseq2.txt", recursive=T, full.names=T)
diffFiles = data.frame(File=diffFiles, Filename=basename(diffFiles))
diffFiles$ComparisonName = gsub(".results_deseq2.txt", "", diffFiles$Filename)


## Loop over every comparison
outFilename = ""
allRegsMV = data.frame()
for(curComparison in unique(diffFiles$ComparisonName))
{
    
    curDiffFile = diffFiles$File[ diffFiles$ComparisonName == curComparison]
    diffPeaks = read.tsv2(curDiffFile)
    
    counts = setRemoveRownames(diffPeaks[, c(1, grep(".raw", colnames(diffPeaks)))], 1)
	## Optional: Filter by low abundance, or "singleton" peaks. If this is not needed, then counts can be a dummy matrix of zeroes
	if(filterPeaks)
	{
		keepPeak = (rowSums(counts > 20) > 2)
		counts = counts[ keepPeak, ]
		diffPeaks = diffPeaks[ keepPeak, ]
	}
    rownames(diffPeaks) = diffPeaks$Geneid
	design = data.frame(colnames(counts), Type="A")
	peaksGR = makeGRangesFromDataFrame(diffPeaks, seqnames.field = "Chr", start.field = "Start", end.field = "End", keep.extra.columns = F, ignore.strand = T)	
    
    ## Match TF motifs inside each peak and make a binary matrix, and also get the GC content of every peak for correction purposes
	fragment_counts = SummarizedExperiment(assays = list(counts = as.matrix(counts)), rowRanges = peaksGR, colData = design)
	fragment_counts2 = addGCBias(fragment_counts, genome = Genome)
	gcBias = data.frame(PeakID=diffPeaks$Geneid, GC_Bias=fragment_counts2@rowRanges$bias)
	filtered_counts = filterPeaks(fragment_counts2)
	motifs_JASPAR2018 = jaspar()
	motifDat = matchMotifs(motifs_JASPAR2018, filtered_counts@rowRanges, genome = Genome)
	motifDatAssay = data.frame(as.matrix(assay(motifDat)))
	
	## Optional: Filter the TFs to ones that are expressed in the tissue based on RNAseq. This may require manual oversight to match the TF names to gene names based on the organism.
	## ....
	
	
	olap = motifDatAssay*1		
	gcBiasCorrection = c(TRUE, FALSE)
	for(CorrectGCBias in gcBiasCorrection)
	{
		print(paste0("Running multivariate TF analysis for ", curComparison))
		curDat = diffPeaks[ diffPeaks$ComparisonName %in% curComparison, ]
		peakData = merge(curDat[, c("Geneid", "log2FoldChange")], olap, by.x="Geneid", by.y="row.names")
		
		if(CorrectGCBias) ## If GC correction is on, add another covariate
		{
			peakData = merge(peakData, gcBias, by.x="Geneid", by.y="PeakID")
		}
		
		## Multivariate testing with all TFs
		if(F)
		{
			tfMv = lm(log2FoldChange ~ ., data = peakData[, -1])		
			tfMvSs = coeftest(tfMv, vcov. = vcovHC)
			tfMvSs2 = as.data.frame(tfMvSs)
			tfMvSs2 = addRownameColumn(tfMvSs2, "TF")
			colnames(tfMvSs2) = c("TF", "estimate", "std.error", "t", "p.value")
			Multiplier = qnorm(1 - 0.05 / 2)
			tfMvSs2$p.adj = p.adjust(tfMvSs2$p.value, "BH")
			tfMvSs2$ymin = tfMvSs2$estimate - Multiplier * tfMvSs2$std.error
			tfMvSs2$ymax = tfMvSs2$estimate + Multiplier * tfMvSs2$std.error
		}else ## Version with lm_robust
		{
			tfMv = estimatr::lm_robust(log2FoldChange ~ ., data = peakData[, -1])	
			tfMvSs2 = data.frame(TF=names(tfMv$coefficients), estimate=tfMv$coefficients, std.error=tfMv$std.error, t=tfMv$statistic, p.value=tfMv$p.value, p.adj = p.adjust(tfMv$p.value, "BH"), ymin=tfMv$conf.low, ymax=tfMv$conf.high)
		}
		
		tfMvSs2 = tfMvSs2[ ! tfMvSs2$TF %in% c("(Intercept)", "GC_Bias"), ]
		allRegsMV = rbind(allRegsMV, data.frame(ComparisonName = curComparison, GC_BiasCorrection=CorrectGCBias, tfMvSs2))
	}
}
write.tsv(allRegsMV, paste0(outFilename, "_JASPAR_motifMultivariateLinearModel.txt"), row.names = F)


## Plot results as bars
allRegsMV2 = allRegsMV[ order(allRegsMV$p.adj), ]
allRegsMV2 = allRegsMV2[ allRegsMV2$p.adj < 0.05, ]

topN = 30
for(CorrectGCBias in c(TRUE, FALSE))
{
    for(curComparison in unique(allRegsMV2$ComparisonName))
    {
        curReg = allRegsMV2[ allRegsMV2$ComparisonName %in% curComparison & allRegsMV2$GC_BiasCorrection %in% CorrectGCBias, ]
        if(nrow(curReg) < 1) next
        
        curReg = curReg[ 1:min(topN, nrow(curReg)), ]
        curReg$TF = splitGet(curReg$TF, "_", 2)
        curReg$TF = factor(curReg$TF, levels = curReg$TF[order(curReg$estimate)])
        curReg$p.adj[curReg$p.adj == 0] = 1e-300
        ggplot(curReg, aes(x=TF, y=estimate, fill=-log10(p.adj))) + geom_hline(yintercept=0, color="black") + geom_bar(stat="identity") + geom_errorbar(aes(ymin=ymin, ymax=ymax), color="black", width=0.2) + theme_cem + coord_flip() + scale_fill_gradientn(colors=rev(viridis::inferno(10)[2:8])) + ggtitle(paste0("Multivariate TF regulatory potential\n", curComparison))
        ggsave(paste0(outFilename, "_JASPAR_motifMultivariateLinearModel_", make.names(curComparison), ifelse(CorrectGCBias, "_GCbiasCorrection", "_NoCorrection"), ".png"), width=12, height=12)
        ggsave(paste0(outFilename, "_JASPAR_motifMultivariateLinearModel_", make.names(curComparison), ifelse(CorrectGCBias, "_GCbiasCorrection", "_NoCorrection"), ".pdf"), width=12, height=12)
    }
}

## Plot results as dots for selected TFs at q < 0.01
{
	allRegsMV  = read.tsv2(paste0(outFilename, "_JASPAR_motifMultivariateLinearModel.txt"))
	allRegsMV2 = allRegsMV[ order(allRegsMV$p.adj), ]
	allRegsMV2 = allRegsMV2[ allRegsMV2$p.adj < 0.01, ]
	selTF = c("MA0655.1_JDP2", "MA0476.1_FOS", "MA0490.1_JUNB", "MA0808.1_TEAD3", "MA0089.1_MAFG..NFE2L1", "MA0833.1_ATF4", "MA0604.1_Atf1", "MA0018.3_CREB1", "MA1099.1_Hes1", "MA1148.1_PPARA..RXRA", "MA0599.1_KLF5", "MA0038.1_Gfi1", "MA0512.2_Rxra", "MA0050.2_IRF1", "MA0114.3_Hnf4a", "MA1109.1_NEUROD1", "MA0139.1_CTCF", "MA0119.1_NFIC..TLX1")

	allRegsMV3 = allRegsMV2[ allRegsMV2$TF %in% selTF, ]
	CorrectGCBias = T

	curReg = allRegsMV3[  allRegsMV3$GC_BiasCorrection %in% CorrectGCBias, ]
	curReg = curReg[ grep("DESeq", curReg$ComparisonName), ]
	curReg = curReg[ grep("AA", curReg$ComparisonName), ]

	curReg$TF = splitGet(curReg$TF, "_", 2)
	curReg$TF = factor(curReg$TF, levels = unique(curReg$TF[order(sign(curReg$estimate)* -log10(curReg$p.adj))]))
	curReg$p.adj[curReg$p.adj == 0] = 1e-300
	curReg$ComparisonName = splitGet(curReg$ComparisonName, ":", 2)
		
	ggplot(curReg, aes(x=TF, y=ComparisonName, fill=estimate, size=-log10(p.adj))) + geom_point(shape=21) + theme_cem + coord_flip() + scale_fill_gradientn(colors=rev(brewer.pal(9, "RdYlBu")),  limits = c(-0.1,0.1), oob=scales::squish) + ggtitle(paste0("Multivariate TF regulatory potential")) + scale_size(range=c(2.5,7), breaks=c(2,5,10,20,30,40,50)) + facet_wrap(~ComparisonName)
	ggsave(paste0(outFilename, "_JASPAR_motifMultivariateLinearModel_SelectedDot_AA_", ifelse(CorrectGCBias, "_GCbiasCorrection", "_NoCorrection"), ".png"), width=4, height=5)
	ggsave(paste0(outFilename, "_JASPAR_motifMultivariateLinearModel_SelectedDot_AA_", ifelse(CorrectGCBias, "_GCbiasCorrection", "_NoCorrection"), ".pdf"), width=4, height=5)
}
