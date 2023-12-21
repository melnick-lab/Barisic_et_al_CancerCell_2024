library(GenomicRanges)
library(UpSetR)
library(grid)
library(Rtsne)
library(tidyr)
library(pheatmap)
library(purrr)
library(cowplot)
library(tidyverse)
library(magrittr)

# Helper functions
{
           
    scaledRank = function(x) { return(scales::rescale(rank(x))); } 


	getListIndex = function(x, n)
	{
		if(n > length(x)) return("") else return(x[[n]])
	}
	splitGet=function (strings, sep, n)
	{
		splitVals = strsplit(as.character(strings), sep)
		return ( sapply(splitVals, getListIndex, n) )
	}

	splitGetFromEnd=function (strings, sep, n)
	{
		splitVals = strsplit(strings, sep)
		strLength = sapply(splitVals,length) - n + 1
		return(mapply("[[", splitVals, strLength))
	}

	
	read.tsv2 = function(file, header=T, row.names=F, ...)
	{
		if(row.names==F)
		{
			library(data.table)
			df = fread(file, stringsAsFactors=F, na.strings = c("NA","#N/A","N/A", "NULL"), data.table=F, ...)
			colnames(df) = make.names(colnames(df))
			return(df)
		}else
		{
			return(read.table(file, sep="\t", header=header, stringsAsFactors=F,  na.strings = c("NA","#N/A","N/A", "NULL"), ...))
		}
	}

	write.tsv = function(object, file, row.names=T, quote=F, col.names=NA)
	{
		if(row.names==F & is.na(col.names)) col.names=T
		write.table(object, file, sep="\t", quote=quote, row.names=row.names, col.names=col.names)
	}

    getChipDesign = function(rootPath, peakType="narrowPeak")
    {
        design = list.files(paste0(rootPath, "/bwa/mergedLibrary/bigwig/"), "*(bigWig|bigwig|bw)$", full.names = T)
        design = data.frame(FileName = basename(design), BigwigFile = design)
        design$Name = gsub(".mLb.clN", "", design$FileName)
        design$Name = gsub("(.bw|.bigWig|.bigwig)$", "", design$FileName)
        design$Sample = design$Name
        design = design[, c("Sample", "Name", "BigwigFile")]

        bedFiles = list.files(paste0(rootPath, "/bwa/mergedLibrary/macs/"), paste0("*.", peakType, "$"), full.names = T, recursive=T)
        bedFiles = data.frame(BedFile = bedFiles)
        bedFiles$Name = basename(bedFiles$BedFile)
        bedFiles$Name = gsub("(.mLb.clN_peaks|_peaks)", "", bedFiles$Name)
        bedFiles$Name = gsub(paste0(".", peakType, "$"), "", bedFiles$Name)
        bedFiles$Factor = splitGet(bedFiles$Name, "[_-]", 2)
        bedFiles$PeakCaller = paste0("macs_", peakType)
        

        bamFiles = list.files(paste0(rootPath, "/bwa/mergedLibrary/"), "*.bam$", full.names = T)
        bamFiles = data.frame(BamFile = bamFiles)
        bamFiles$Name = basename(bamFiles$BamFile)
        bamFiles$Name = gsub(".mLb.clN.sorted.bam", "", bamFiles$Name)

        design = merge(design, bamFiles, by.x="Name", all.x=T)
        design = merge(design, bedFiles, by.x="Name", all.x=T)
        #design$Sample = paste(splitGet(design$Name, "_", 1), splitGet(design$Name, "_", 2), splitGet(design$Name, "_", 3), sep="_")

        design = design[ !is.na(design$BedFile) & !is.na(design$BigwigFile) & !is.na(design$BamFile), ]
        
        return(design)
    }

    getCnrDesign = function(rootPath)
    {
        design = list.files(paste0(rootPath, "/cnr_output/S4_A_aln_bigWig.all/"), "*bigWig$", full.names = T)
        design = data.frame(FileName = basename(design), BigwigFile = design)
        design$Name = gsub("(.bigWig)", "", design$FileName)
        design$Sample = design$Name
        design = design[, c("Sample", "Name", "BigwigFile")]

        bedFiles = list.files(paste0(rootPath, "/cnr_output/S5_B_peaks_seacr.all/"), "*.all.peaks.seacr.stringent.bed$", full.names = T)
        bedFiles = data.frame(BedFile = bedFiles)
        bedFiles$Name = basename(bedFiles$BedFile)
        bedFiles$Name = gsub(".all.peaks.seacr.stringent.bed", "", bedFiles$Name)
        bedFiles$Factor = splitGet(bedFiles$Name, "_", 1)
        bedFiles$PeakCaller = "seacr"
        

        bamFiles = list.files(paste0(rootPath, "/cnr_output/S2_C_aln_bdg.all/"), "*.bam$", full.names = T)
        bamFiles = data.frame(BamFile = bamFiles)
        bamFiles$Name = basename(bamFiles$BamFile)
        bamFiles$Name = gsub("_sort_byname.bam", "", bamFiles$Name)

        design = merge(design, bamFiles, by.x="Name", all.x=T)
        design = merge(design, bedFiles, by.x="Name", all.x=T)
        design$Sample = paste(splitGet(design$Name, "_", 1), splitGet(design$Name, "_", 2), splitGet(design$Name, "_", 3), sep="_")

        design = design[ !is.na(design$BedFile) & !is.na(design$BigwigFile) & !is.na(design$BamFile), ]
        
        return(design)
    }

    getAtacDesign = function(rootPath, assay = "ATACseq", peakType="narrowPeak")
    {
        design = list.files(paste0(rootPath, "/bwa/mergedLibrary/bigwig/"), "*(bigWig|bigwig|bw)$", full.names = T)
        design = data.frame(FileName = basename(design), BigwigFile = design)
        design$Name = gsub(".mLb.clN", "", design$FileName)
        design$Name = gsub("(.bw|.bigWig|.bigwig)$", "", design$Name)
        design$Sample = design$Name
        design = design[, c("Sample", "Name", "BigwigFile")]

        bedFiles = list.files(paste0(rootPath, "/bwa/mergedLibrary/macs/"), paste0("*.", peakType, "$"), full.names = T, recursive=T)
        bedFiles = data.frame(BedFile = bedFiles)
        bedFiles$Name = basename(bedFiles$BedFile)
        bedFiles$Name = gsub("(.mLb.clN_peaks|_peaks)", "", bedFiles$Name)
        bedFiles$Name = gsub(paste0(".", peakType, "$"), "", bedFiles$Name)
        bedFiles$Factor = assay


        bedFiles$PeakCaller = "macs"

        bamFiles = list.files(paste0(rootPath, "/bwa/mergedLibrary/"), "*.bam$", full.names = T)
        bamFiles = data.frame(BamFile = bamFiles)
        bamFiles$Name = basename(bamFiles$BamFile)
        bamFiles$Name = gsub(".mLb.clN.sorted.bam", "", bamFiles$Name)

        design = merge(design, bamFiles, by.x="Name", all.x=T)
        design = merge(design, bedFiles, by.x="Name", all.x=T) 

        design = design[ !is.na(design$BedFile), ]

        return(design)
    }
    getAtacPooledDesign = function(rootPath, assay = "ATACseq", peakType="narrowPeak")
    {
        design = list.files(paste0(rootPath, "/bwa/mergedReplicate/bigwig/"), "*(bigWig|bigwig|bw)$", full.names = T)
        design = data.frame(FileName = basename(design), BigwigFile = design)
        design$Name = gsub(".mRp.clN", "", design$FileName)
        design$Name = gsub("(.bw|.bigWig|.bigwig)$", "", design$Name)
        design$Sample = design$Name
        design = design[, c("Sample", "Name", "BigwigFile")]

        bedFiles = list.files(paste0(rootPath, "/bwa/mergedReplicate/macs/"), paste0("*.", peakType, "$"), full.names = T, recursive=T)
        bedFiles = data.frame(BedFile = bedFiles)
        bedFiles$Name = basename(bedFiles$BedFile)
        bedFiles$Name = gsub("(.mRp.clN_peaks|_peaks)", "", bedFiles$Name)
        bedFiles$Name = gsub(paste0(".", peakType, "$"), "", bedFiles$Name)
        bedFiles$Factor = assay


        bedFiles$PeakCaller = "macs"

        bamFiles = list.files(paste0(rootPath, "/bwa/mergedReplicate/"), "*.bam$", full.names = T)
        bamFiles = data.frame(BamFile = bamFiles)
        bamFiles$Name = basename(bamFiles$BamFile)
        bamFiles$Name = gsub(".mRp.clN.sorted.bam", "", bamFiles$Name)

        design = merge(design, bamFiles, by.x="Name", all.x=T)
        design = merge(design, bedFiles, by.x="Name", all.x=T) 

        design = design[ !is.na(design$BedFile), ]

        return(design)
    }
     
    getMarkAnnotations = function(filePath, hasType=T, annotationPrefix="",  chrCol=1, startCol=2, endCol=3, typeCol=4)
    {
        curAnn = read.tsv2(filePath)
        if(hasType)
        {
            curAnn2 = curAnn[, c(chrCol, startCol, endCol, typeCol)]
            colnames(curAnn2) = c("chr", "start", "end", "Type")
            annPrefix = ""
            if(annotationPrefix != "" & !grepl("_$", annotationPrefix)) annPrefix = paste0(annotationPrefix, "_")
            curAnn2$Type = paste0(annPrefix, curAnn2$Type)
        }else
        {
            curAnn2 = curAnn[, c(chrCol, startCol, endCol)]
            colnames(curAnn2) = c("chr", "start", "end")
            annPrefix = annotationPrefix
            if(annotationPrefix == "") annPrefix="TRUE"
            curAnn2$Type = annPrefix
        }
        return(curAnn2)
    }

    getMarkAnnotationsMulti = function(fileList, chrCol=1, startCol=2, endCol=3)
    {
        curAnn = list()
        for(i in 1:length(fileList))
        {
            curPath = fileList[[i]]
            curName = names(fileList)[i]
            if(is.null(curName)) curName = file_path_sans_ext(basename(curPath))
            curAnn[[curName]] = getMarkAnnotations(curPath, chrCol=chrCol, startCol=startCol, endCol=endCol, hasType=F, annotationPrefix=curName)
        }
        return(curAnn)
    }

    smoothScorePredict = function(trainDat, newGrid, knn=30)
    {
        scorePredicted = data.frame(score = as.numeric(NA), 
                                    x = newGrid[, 1], 
                                    y = newGrid[, 2])
        # run KKNN
        scoreKknn = kknn::kknn(score ~ ., 
                            train = trainDat, 
                            test = scorePredicted, 
                            kernel = "gaussian", 
                            k = knn)

        scorePredicted %<>% mutate(score = fitted(scoreKknn))
        return(scorePredicted)
    }

    smoothScore2d = function(score, x, y, type=NULL, numGrid = 100, knn = 100, m = 2, expand=0.05, xrng=NULL, yrng=NULL)
    {
        library(ash)
        curDat = data.frame(score=score, x=x, y=y)
        if(is.null(xrng)) xrng = range(curDat$x)
        if(is.null(yrng)) yrng = range(curDat$y)
        xdiff = xrng[2] - xrng[1]
        ydiff = yrng[2] - yrng[1]
        xrng[1] = xrng[1] - xdiff*expand
        xrng[2] = xrng[2] + xdiff*expand
        yrng[1] = yrng[1] - ydiff*expand
        yrng[2] = yrng[2] + ydiff*expand
        bins = bin2(cbind(curDat$x, curDat$y), ab = rbind(xrng, yrng), nbin = c(numGrid, numGrid))
        binCounts = ash2(bins, m = c(m, m))
        gridDat = data.frame( expand.grid( x = binCounts$x, y = binCounts$y), density = melt(binCounts$z)[,3] )
        gridDat2 = gridDat[ gridDat$density > 0, ]

        if(is.null(type))
        {
            return( smoothScorePredict(curDat, gridDat2, knn = knn) )
        }else
        {
            curDat$type = type
            allPredicted = ddply(curDat, "type", function(x) { smoothScorePredict( x[, colnames(x) != "type"], gridDat2, knn = knn) } )
            return(allPredicted)
        }
    }
	
	## Augment Plot functions and other related functions are borrowed/adapted from the Seurat package used under MIT License. 
	## Used to rasterize the plot with high number of points so they can be visualized better in PDF format while keeping the axes in vector
	AugmentPlot2 = function(plot1, imgFile)
    {
      range.values = c(
        ggplot_build(plot = plot1)$layout$panel_ranges[[1]]$x.range,
        ggplot_build(plot = plot1)$layout$panel_ranges[[1]]$y.range
      )
      img = png::readPNG(source = imgFile)
      p1mod = plot1 + annotation_raster(
        img,
        xmin = range.values[1],
        xmax = range.values[2],
        ymin = range.values[3],
        ymax = range.values[4]
      )
      return(p1mod)
    }

    GetXYAesthetics = function(plot, geom = 'GeomPoint', plot.first = TRUE) 
    {
      geoms = sapply(
        X = plot$layers,
        FUN = function(layer) {
          return(class(x = layer$geom)[1])
        }
      )
      # handle case where raster is set to True
      if (geom == "GeomPoint" && "GeomScattermore" %in% geoms){
        geom = "GeomScattermore"
      }
      geoms = which(x = geoms == geom)
      if (length(x = geoms) == 0) {
        stop("Cannot find a geom of class ", geom)
      }
      geoms = min(geoms)
      if (plot.first) {
        # x = as.character(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)[2]
        x = as_label(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)
        # y = as.character(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)[2]
        y = as_label(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)
      } else {
        x = as_label(x = plot$layers[[geoms]]$mapping$x %||% plot$mapping$x)
        y = as_label(x = plot$layers[[geoms]]$mapping$y %||% plot$mapping$y)
      }
      return(list('x' = x, 'y' = y))
    }

    NoLegend = function(...) 
    {
      no.legend.theme = theme(
        # Remove the legend
        legend.position = 'none',
        # Validate the theme
        validate = TRUE,
        ...
      )
      return(no.legend.theme)
    }

    NoAxes = function(..., keep.text = FALSE, keep.ticks = FALSE) 
    {
      blank = element_blank()
      no.axes.theme = theme(
        # Remove the axis elements
        axis.line.x = blank,
        axis.line.y = blank,
        # Validate the theme
        validate = TRUE,
        ...
      )
      if (!keep.text) {
        no.axes.theme = no.axes.theme + theme(
          axis.text.x = blank,
          axis.text.y = blank,
          axis.title.x = blank,
          axis.title.y = blank,
          validate = TRUE,
          ...
        )
      }
      if (!keep.ticks){
        no.axes.theme = no.axes.theme + theme(
          axis.ticks.x = blank,
          axis.ticks.y = blank,
          validate = TRUE,
          ...
        )
      }
      return(no.axes.theme)
    }

    AugmentPlot = function(plot, width = 10, height = 10, dpi = 300) 
    {
      pbuild.params = ggplot_build(plot = plot)$layout$panel_params[[1]]
      range.values = c(
        pbuild.params$x.range,
        pbuild.params$y.range
      )
      xyparams = GetXYAesthetics(
        plot = plot,
        geom = class(x = plot$layers[[1]]$geom)[1]
      )
      title = plot$labels$title
      tmpfile = tempfile(fileext = '.png')
      ggsave(
        filename = tmpfile,
        plot = plot + NoLegend() + NoAxes() + theme(plot.title = element_blank()),
        width = width,
        height = height,
        dpi = dpi
      )
      img = png::readPNG(source = tmpfile)
      file.remove(tmpfile)
      blank = ggplot(
        data = plot$data,
        mapping = aes_string(x = xyparams$x, y = xyparams$y)
      ) + geom_blank()
      blank = blank + plot$theme + ggtitle(label = title)
      blank = blank + annotation_raster(
        raster = img,
        xmin = range.values[1],
        xmax = range.values[2],
        ymin = range.values[3],
        ymax = range.values[4]
      )
      return(blank)
    }

    AugmentedGG = function(ggPlot, width=12, height=9, dpi=150)
    {
        previous_call = blank_call = png_call =  ggPlot #match.call()
        blank_call$pt.size = -1
        blank_call$do.return = TRUE
        blank_call$vector.friendly = FALSE
        png_call$no.axes = TRUE
        png_call$no.legend = TRUE
        png_call$do.return = TRUE
        png_call$vector.friendly = FALSE
        png_call$plot.title = NULL
        blank_plot = eval(blank_call, sys.frame(sys.parent()))
        png_plot = eval(png_call, sys.frame(sys.parent()))
        png.file = paste0(tempfile(), ".png")
        ggsave(
          filename = png.file,
          plot = png_plot,
          width = width,
          height = height,
          dpi = dpi
        )
        to_return = AugmentPlot(plot1 = blank_plot, imgFile = png.file)
    }

    AugmentPlotList = function(ggPlotList, width=12, height=9, dpi=150)
    {
        newList = list()
        for(curPlot in ggPlotList)
        {
            if(is.null(curPlot))
            {
                newList = append(newList, list(x=NULL))
            }else
            {
                curAug = AugmentPlot(curPlot, width, height, dpi)
                newList[[ length(newList) + 1]] = curAug
            }
        }
        return(newList)
    }
}
 
setwd("~/Projects/Arid1a/MultiMarkModel_RIVA")

chrList = paste0("chr", c(1:22, "X"))
genomeRef = "hg38"
GenomeRef = genomeRef

if(GenomeRef == "hg38")
{
    GenomeDat = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    greatDomainsFile = "~/great/GREAT_hg38_gencodeSelected_2020_RegDom.txt"
}else if(GenomeRef == "mm10")
{
    GenomeDat = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    greatDomainsFile = "~/great/GREAT_mm10_gencodeSelected_2020_RegDom.txt"
}

## Peak merging parameters
minPeakSize = 250 # Any peak smaller than this will be extended on both sides. Note that the merged/split peaks may result in regions be smaller than this.
maxPeakSize = 100000 # Any peak larger than this will not be allowed to be merged in the second step, but they won't actively be filtered out if they have enough replicate overlaps
minPeakMergeOverlapPercent = 0.5 # 50% overlap to create merged peaks
minPeakOverlapBp = 50 # 50bp overlap while counting which merged peaks overlap with which bed peaks
minPeakReps = 3 # min number of replicates to overlap with a merged peak to keep. Depends on the number of input samples and factors.


outFilename = "MultiMarkModel_RIVA_PeakMerge"
outTitle = "RIVA MMM"


##### Input Data #####
{
    ## RNAseq as readout
	## This uses DESeq2 output to correlate differential gene expression with the chromatin clusters
    rna = read.tsv2("../RNAseq_RIVA/RIVA_DEG.txt")
    rna = rna[ rna$Coefficient %in% "Het vs WT, RIVA postselection", ]
    rna$Comparison = "Arid1aHetVsWT"
    rna = rna[, c("Gene", "Comparison", "baseMean", "log2FoldChange", "padj")]

    ## Create a design matrix that has paths for peak beds, bam files, and other info (Mark/Factor, Genotype...) for all samples.
	## For each Factor (TF or Histone Mark name) you should have bed files for peaks and bam files for quantifications, can have multiple replicates per mark 
    
    # Peaks/bams from ATACseq
    designAtac = getAtacDesign("../ATACseq_RIVA/results/")
    designAtac$CellType = "RIVA"

    # Peaks/bams from 20+X Ab C&R 
    designCnr1 = getCnrDesign("../CutRun_RIVA/cnr")
    designCnr1$CellType = "RIVA"
    designCnr1$Genotype = splitGet(designCnr1$Name, "_", 2)
	
	# Peaks/bams from ChIPseq
	designCh1 = getChipDesign("../ChIPseq_RIVA/results")
	designCh1$CellType = "RIVA"
	designCh1$Factor = splitGet(designCh1$Name, "_", 1)
	designCh1$Genotype  = splitGet(designCh1$Name, "_", 2)

    # Combined list for MMM
	design = rbind.fill(designAtac, designCnr1, designCh1)
	design$Type = design$Genotype
	design$Type2 = paste0(design$Factor, "_", design$Type)
	design$Sample = design$Name
	design = design[ order(design$Factor), ]
    
}

## use all factors for clustering...
# factorsForClustering = unique(design$Factor) 

## OR
## Subset factors for clustering as needed. The rest will be quantified
## but will not participate in peak/region definition, or cluster definition. 
## They will be plotted as read-outs (along with RNAseq and annotations).
factorsForClustering = c("ATACseq", "H3K27ac", "H3K27me3")

## "WT" Genotype will be quantified and added by default.
wtColumn = "WT"

## Which genotypes to compare against WT while clustering the signal. 
## WT base signal plus the log2FC of these genotypes against WT will be added as a signal for clustering
## Log2FC values will also be plotted in the outputs.
genotypeList = c("Arid1aHet", "Arid1aKD")

## !!! Note that grep is used for finding columns, so if the name of the "WT" column and other genotypes are similar this will cause issues (e.g. "Bcl6" as baseline, "Bcl6KO" as mutant will not work)



##### External Annotation Data and Colors #####
{
    # Custom annotations for overlapping with clusters
    # If the file has multiple types of ranges define by hasType and Type columns.
    # Note: If a peak overlaps with multiple ranges of different types the behaviour is undefined, it will count against only one. Use combinatorial function to see multiple types.
    tobiasPath1 = "../ATACseq_RIVA/TOBIAS/"

    customAnnotations = list()
    customAnnotations[["hg38"]] = list()
    customAnnotations[["hg38"]][["GCB_Enhancer"]] = getMarkAnnotations("~/ReferenceFiles/Peaks/hg38/GCB_H3K27ac.scores.bed", hasType=T, chrCol=1, startCol=2, endCol=3, typeCol=4)


    # Combinatorial custom annotations for overlapping with clusters. Use only on single-type (hasType=F) ranges.
    customCombAnnotations = list()
    customCombAnnotations[["hg38"]] = list()
    customCombAnnotations[["hg38"]][["TOBIAS_TFs_combined"]] = getMarkAnnotationsMulti(list(
                                                        NFKB="TOBIAS_RIVA_NFKB12RELA_combined.bed",
                                                        SPI1_SPIB="/TOBIAS_RIVA_SPI1SPIB_combined.bed"
                                                        ))



    annotationColors = list()
    annotationColors[["GCB_Enhancer"]] = c("No annotation"="#bbbbbb", "promoter"="#a6cee3", "enhancer"="#ff7f00", "SE"="#e31a1c"))
}



#####  Colors  #####
{
    colorsDef = list()

    colorsDef[["Genotype"]] = c(
        "WT" = "#377eb8",
        "Arid1aHet" = "#e41a1c",
        "Arid1aKD" = "#9c1795"
    )

    colorsDef[["Factor"]] = c(brewer.pal(10, "Paired"), brewer.pal(10, "Set3"), iwanthue(10))
    names(colorsDef[["Factor"]]) = unique(design$Factor)


    colorsDiff = c(
        "Arid1aHetvsWT" = "#FF2600"
    )


    ## Color Scales
    colScale = rev(brewer.pal(11, "Spectral"))
    colScaleDiff = c("#2d004b", "#542788", "#8073AC", "#B2ABD2", "#ffffff", "#FDB863", "#E08214", "#B35806", "#7f3b08") 
    colScaleDiff2 = rev(c("#543005", "#8c510a", "#bf812d", "#dfc27d", "#ffffff", "#80cdc1", "#35978f", "#01665e", "#003c30")) 
    colScaleDiff3 = rev(c("#7a0177", "#c51b8a", "#f768a1", "#ffffff", "#d9f0a3", "#78c679", "#006837"))
    colScaleDiffRNA = c("royalblue4", "royalblue4", "royalblue4", "royalblue2", "cornflowerblue", "lightskyblue",  "#ffffff", "darkorange", "firebrick2", "firebrick3", "firebrick4", "firebrick4", "firebrick4") 

    ## Add colors to the design matrix
    colorsDefList = list()
    for(curColumn in names(colorsDef))
    {
        curColors = colorsDef[[curColumn]]
        if(curColumn %in% colnames(design))
        {
            design[, curColumn] = factor(design[, curColumn], levels = names(curColors))
            design[, paste0("color", curColumn)] = curColors[ as.numeric(design[, curColumn]) ]

            #colorsDef[[curColumn]] = as.list(curColors)
            ## in case colors above are not named manually name them
            curColorsDf = unique(design[, c(curColumn, paste0("color", curColumn))])
            colorsDefList[[curColumn]] = list()
            colorsDefList[[curColumn]][ as.character(curColorsDf[,1]) ] = curColorsDf[,2]
        }
    }

    colorsDefContinuous = list()
}

########################################
##### Peak calling, quantification #####
########################################

regionsFile = paste0(outFilename, "_obj_allFactorRegions.rds")
regionsFileDF = paste0(outFilename, "_obj_allFactorRegionsDF.rds")



## For next step either manually define common peaks to quantify the signal in
## Or use the code below to automatically define regions based on some criteria
## Merge all factor peaks
{
    curPeakCaller = design$PeakCaller[1]
    if( !file.exists(regionsFileDF) )
    {
        allFactorRegions = list()
        allFactorRegionsDF = list()
        
         
        ReadExtendBed = function(curFile, minPeakSize)
        {
            curFileDat = read.tsv2(curFile, header = F)
            curFileDat = curFileDat[, c(1:3)]
            colnames(curFileDat) = c("chr", "start", "end")
            curFileDat = curFileDat[ curFileDat$chr %in% chrList, ]
            curFileDat$width = curFileDat$end - curFileDat$start
            curFileDat$extend = minPeakSize - curFileDat$width
            curFileDat$extend[ curFileDat$extend < 0] = 0
            curFileDat$start = curFileDat$start - curFileDat$extend/2
            curFileDat$end = curFileDat$end + curFileDat$extend/2
            curFileDat = curFileDat[, c("chr", "start", "end")]
            return(curFileDat)
       }
            
    
        curDesign = design[ design$Factor %in% factorsForClustering, ]
        curRegions = GRanges()
        
        # Combine all peaks from all files
        for(curFile in curDesign$BedFile)
        {
            curFileDat = ReadExtendBed(curFile, minPeakSize)
            curFileGR = makeGRangesFromDataFrame(curFileDat)
            curRegions = c(curRegions, curFileGR)
        }
        curRegions = unique(curRegions)
        
        # Find overlaps between all pairwise set of peaks across replicates and marks
        peakOverlap = findOverlaps(curRegions, curRegions)
        peakOverlapDf = data.frame(queryHits(peakOverlap), subjectHits(peakOverlap), curRegions[queryHits(peakOverlap)], curRegions[subjectHits(peakOverlap)], width(pintersect(curRegions[queryHits(peakOverlap)], curRegions[subjectHits(peakOverlap)])))
        colnames(peakOverlapDf) = c("i", "j", "chr1", "start1", "end1", "width1", "strand1", "chr2", "start2", "end2", "width2", "strand2", "overlapWidth")
        peakOverlapDf = peakOverlapDf[ peakOverlapDf$i != peakOverlapDf$j, ]
        
        # Find which peaks have more than 50% overlap (of both peaks), these will have edges in the graph
        peakOverlapDf2 = peakOverlapDf[ peakOverlapDf$overlapWidth >= minPeakMergeOverlapPercent*peakOverlapDf$width1 & peakOverlapDf$overlapWidth >= minPeakMergeOverlapPercent*peakOverlapDf$width2, ]
        
        library(igraph)
        
        # Create a graph for all the peaks, significantly overlapping peaks will have edges between them
        vertices = data.frame(row.names=1:length(curRegions), i=1:length(curRegions), curRegions)
        edges = data.frame(from=peakOverlapDf2$i, to=peakOverlapDf2$j)
        G = graph_from_data_frame(edges, directed = FALSE, vertices = vertices)
        
        
        # Find connected components in the graph. 
        connectedComponents = components(G)$membership
        vertices$Component = connectedComponents
        
        # merge connected components into merged peaks, non-connected peaks will also be included
        mergedPeaks = vertices %>% group_by(Component) %>% summarise(chr=unique(seqnames), start=min(start), end=max(end))
        mergedPeaks$width = mergedPeaks$end - mergedPeaks$start
        
        
        # mergedPeaks will include peak clusters, but may still have overlaps (e.g. peaks A=10-20, B=8-22, C=9-18, D=10-50, E=11-49. A,B,C are a component. D,E are a component. There is no edge between A-B-C and D-E clusters even though they overlap, because the overlap is not 50%.
        # We will use the boundaries to make a set of disjoint peaks from these clusters
        chromLengths = seqlengths(GenomeDat)
        chromLengths = chromLengths[ names(chromLengths) %in% chrList ]
        

        rangesList = GRanges()

        # Loop through each chromosome and create the ranges from the boundaries
        for (curChr in names(chromLengths)) 
        {
          chromLen = chromLengths[curChr]
          
          curBoundaries = sort(unique(c(mergedPeaks$start[mergedPeaks$chr %in% curChr], mergedPeaks$end[mergedPeaks$chr %in% curChr])))
          coords = c(0, curBoundaries, chromLen)
          
          ranges = GRanges(seqnames=curChr, IRanges(start = coords[-length(coords)], end = coords[-1]))
          rangesList = c(rangesList, ranges)
        }

        # Combine all the ranges into a single object. 
        # Note: This includes non-peak (inbetween) regions as well! Next step will remove those regions
        curRegions = rangesList
             
        # Find overlap of input bed peaks with the merged regions
        overlapDF = as.data.frame( curRegions )
        for(i in 1:nrow(curDesign))
        {
            curFile = curDesign$BedFile[i]
            curName = curDesign$Sample[i]
            
            curFileDat = ReadExtendBed(curFile, minPeakSize)

            curFileGR = makeGRangesFromDataFrame(curFileDat)
            
            overlapDF[, curName ] = 0
            
            overlap = findOverlaps(curFileGR, curRegions, minoverlap=minPeakOverlapBp)
            overlapDF[ subjectHits(overlap), curName ] = 1
        }
        
        rownames(overlapDF) = paste0(overlapDF$seqnames, ".", overlapDF$start, ".", overlapDF$end)
        
        
        
        
        # We may get consecutive peaks that have the same peak overlap pattern. Merge any consecutive peaks if they have the exact or approximately the same overlap profile
        overlapDF2 = overlapDF
        overlapDF2tmp = overlapDF2
        overlapDF2$OverlapCount = rowSums(overlapDF2tmp[, -c(1:5)] == 1)
        overlapDF2$OverlapStatus = do.call(paste0, c(overlapDF2tmp[, -c(1:5)])) 
        overlapDF2$OverlapDiff = 0
        # Bitwise difference/XOR of peak presence
        for(i in 6:(ncol(overlapDF2tmp)))
        {
            curDiff = abs(diff(overlapDF2tmp[, i]))
            overlapDF2$OverlapDiff = overlapDF2$OverlapDiff + c(0, curDiff[-length(curDiff)]) # shift one
        }
        
        # Merge peaks that are mostly similar based on their peak presence patterns
        overlapDF2$OverlapDiffBreakpoint = 0 
        overlapDF2$OverlapDiffBreakpoint[  setdiff(which(overlapDF2$OverlapCount < minPeakReps) + 1, nrow(overlapDF2)+1) ] = 1 # Do not merge any region that has less than minPeakReps (shifted one position for lag, remove the last position if out of bounds)
        overlapDF2$OverlapDiffBreakpoint[ overlapDF2$OverlapDiff >= minPeakReps] = 1 # Do not merge any consecutive peaks that have more than minPeakReps difference in their peak composition either
        overlapDF2$OverlapDiffBreakpoint[ setdiff(which( overlapDF2$width > maxPeakSize) + 1, nrow(overlapDF2)+1) ] = 1 # Do not merge any peaks that are above the max peak size limit (shifted one position for lag, remove the last position if out of bounds)
        
        
        overlapDF2$OverlapDiffGroup = cumsum(overlapDF2$OverlapDiffBreakpoint)
        overlapDF3 = overlapDF2 %>% group_by(seqnames, grp = OverlapDiffGroup) %>% summarise(start = min(start), end = max(end)) %>%  select(-grp)
        overlapDF3$width = overlapDF3$end - overlapDF3$start
        
        
        
        # Repeat the overlap of base bed peaks with the new merged/split peaks
        newPeaks = overlapDF3
        newPeaksGR = makeGRangesFromDataFrame(newPeaks)
        overlapDF = as.data.frame( newPeaks )
        for(i in 1:nrow(curDesign))
        {
            curFile = curDesign$BedFile[i]
            curName = curDesign$Sample[i]
            curFileDat = ReadExtendBed(curFile, minPeakSize)
            curFileGR = makeGRangesFromDataFrame(curFileDat)
            
            overlapDF[, curName ] = 0
            
            overlap = findOverlaps(curFileGR, newPeaksGR, minoverlap=minPeakOverlapBp)
            overlapDF[ subjectHits(overlap), curName ] = 1
        }
        
        
        rownames(overlapDF) = paste0(overlapDF$seqnames, overlapDF$start, overlapDF$end, sep=":")
        # Keep only the regions with more than X replicate overlap
        commonPeaks = overlapDF
        commonPeaks = commonPeaks[, -c(1:4)]
        commonPeaks = commonPeaks[ rowSums(commonPeaks == 1) >= minPeakReps, ]

        overlapFiltered = overlapDF[ rownames(overlapDF) %in% rownames(commonPeaks), ]
        
        allFactorRegionsDF[[ curPeakCaller ]] = overlapFiltered

        saveRDS(allFactorRegionsDF, regionsFileDF)
    }else
    {
        allFactorRegionsDF = readRDS(regionsFileDF)
    }
}

## Quantify merged peaks
{
	overlapDF = allFactorRegionsDF[[ curPeakCaller ]]
	curCols = setdiff(colnames(overlapDF), c("seqnames", "start", "end", "width", "strand"))
	
    curDesign = design # This time quantify all factors
	
	ctFile = paste0(outFilename, "_PeakCounts_AllFactorPeaks.txt")
	if(file.exists(ctFile))
	{
		countsDat = read.tsv2(ctFile)
		countsDat = setRemoveRownames(countsDat, 1)
	}else
	{
		curAnnot = overlapDF[, 1:4]
		colnames(curAnnot) = c("Chr", "Start", "End", "Width")
        curAnnot$Strand = "."
		curAnnot$GeneID = paste0(curAnnot$Chr, ".", curAnnot$Start, ".", curAnnot$End)
		
		counts = Rsubread::featureCounts(curDesign$BamFile, annot.ext=curAnnot, minOverlap=10, isPairedEnd=T, requireBothEndsMapped=T, nthreads=20) # change parameters as necessary
		countsDat = counts$counts
		colnames(countsDat) = curDesign$Name
		write.tsv(countsDat, ctFile, row.names=T)
	}
		

    if(!file.exists(paste0(outFilename, "_AllFactorPeaks_log2cpmData.txt")))
    {
        counts = read.tsv2(paste0(outFilename, "_PeakCounts_AllFactorPeaks.txt"))
        counts = setRemoveRownames(counts, 1)
             
        curDesign = design
        curDesign$SampleID = make.names(curDesign$Sample)
        count3 = counts[, curDesign$SampleID]
        
        coldata = data.frame(row.names = colnames(count3), temp = rep("A", ncol(count3)))
        dds = DESeqDataSetFromMatrix(countData = round(count3), colData = coldata, design = ~1)
        rldDat = log2(fpm(dds) + 1)
        
        rldDat2 = rldDat[, curDesign$SampleID]
        write.tsv(rldDat2, paste0(outFilename, "_AllFactorPeaks_log2cpmData.txt"))
        rldDatsAll = rldDat2
    }
    
}

## Average replicates
{
    rldDatsAll = read.tsv2(paste0(outFilename, "_AllFactorPeaks_log2cpmData.txt"))
    rldDatsAll = setRemoveRownames(rldDatsAll, 1)
    
    rldDatsAll2 = rldDatsAll   
    rldM = melt(addRownameColumn(rldDatsAll2, "PeakID"), variable.name="SampleID", value.name="Score")
    rldM$SampleType = mapvalues(rldM$SampleID, make.names(design$Name), make.names(paste0(design$Factor, "_", design$Genotype)))
    
    rldMS = rldM %>% group_by(PeakID, SampleType) %>% summarise(Score=median(Score, na.rm=T))
    rldW = dcast(rldMS, PeakID~SampleType, value.var="Score")
    rldW = setRemoveRownames(rldW, 1)
}

## Filter to rows that have a minimum signal level (log2 FPM) in at least one mark. 
rldW = rldW[ rowSums( rldW[, grep(wtColumn, colnames(rldW), value=T)] >= 2) >= 1, ]

## Add WT as baseline, and each log2FC difference as another signal, rather than including baselines for all genotypes
## Calculate approximate difference in log scale
## (since peaks are merged/collapsed across marks we cannot take the pre-calculated shrunken log2FC any more)
{
    dat = rldW[, grep(wtColumn, colnames(rldW), value=T)]
    for(curFactor in make.names(unique(design$Factor)))
    {
        for(curGenotype in genotypeList)
        {
            if(paste0(curFactor, "_", curGenotype) %in% colnames(rldW) & paste0(curFactor, "_", wtColumn) %in% colnames(rldW))
            {
                dat[, paste0(curFactor, "_", curGenotype, "vs", wtColumn)] = rldW[, paste0(curFactor, "_", curGenotype)] - rldW[, paste0(curFactor, "_", wtColumn)] 
            }
        }
    }
}


datOrg = dat
datM = melt(addRownameColumn(dat, "PeakID"), "PeakID", variable.name="SampleType", value.name="Score")
selFactorCols = grep(paste0("^(", paste0(factorsForClustering, collapse="|"), ")_"), colnames(dat), value=T)
datScale = scale(dat[, selFactorCols])


if(!file.exists(paste0(outFilename, "_scaledAbs_coords.txt")))
{  
    datDimRed_scaled = Rtsne(datScale, verbose = TRUE, dims = 2, initial_dims = min(100, ncol(datScale)), perplexity = 50, theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000, is_distance= FALSE)
    datDimRed2_scaled = data.frame(dat, datDimRed_scaled$Y)
    write.tsv(datDimRed2_scaled, paste0(outFilename, "_scaledAbs_coords.txt"), row.names=T)
}else
{
    datDimRed2_scaled = read.tsv2(paste0(outFilename, "_scaledAbs_coords.txt"))
    datDimRed2_scaled = setRemoveRownames(datDimRed2_scaled, 1)
}

datM_scaled = merge(datM, datDimRed2_scaled[, c("X1", "X2")], by.x="PeakID", by.y="row.names")


######################
##### Clustering #####
######################


k = 16
curOutName = paste0(outFilename, "_2023.05.16_k", k)
Qlim = function(x) return(quantile(x, c(0.02, 0.98))) # scale by top-bottom Nth percentile




if(!file.exists(paste0(curOutName, "_ordered.txt")))
{
    km = kmeans(datScale, k, iter.max = 30)
    #pm = cluster::pam(datScale, k, metric = "euclidean", stand = FALSE)
    #pmk = cluster::pamk(datScale, metric = "euclidean", stand = FALSE)

    centerOrder = getTSPOrder(km$centers)
    names(centerOrder) = paste0("C_", sprintf("%02d", 1:length(centerOrder)))
    dat$Cluster = mapvalues(km$cluster, centerOrder, names(centerOrder))

    ## Order the peaks
    datL5 = data.frame()
    clustBreaks = c()
    for(i in sort(unique(dat$Cluster)))
    {
        curDat = dat[ dat$Cluster == i, !colnames(dat) %in% c("X1", "X2", "Cluster")]
        curDatAll = dat[ dat$Cluster == i, ]
        if(nrow(curDat)==0) next
        clustBreaks = c(clustBreaks, nrow(curDat))

        if(F)
        {
            curDat2 = Rtsne(curDat, verbose = TRUE, dims = 1, initial_dims = 100, perplexity = 30, check_duplicates = FALSE, pca = TRUE, max_iter = 1000, is_distance= FALSE, theta = 0) #
            curDatAll = curDatAll[ order(curDat2$Y), ]
        }else
        {
            dists = dist(curDat)
            hc = hclust(dists, method="ward.D")
            curDatAll = curDatAll[hc$order, ]
        }
        
        datL5 = rbind(datL5, curDatAll)
    }

    write.tsv(datL5, paste0(curOutName, "_ordered.txt"))
}else
{
    datL5 = setRemoveRownames(read.tsv2(paste0(curOutName, "_ordered.txt")), 1)
    dat$Cluster = mapvalues(rownames(dat), rownames(datL5), datL5$Cluster)
}

dat2 = data.frame(dat, datDimRed2_scaled[, c("X1", "X2")])

bedOut = paste0(curOutName, "_clusters.bed")
bedHeader = paste0("track name=\"", as.character(outTitle), "\" visibility=2 itemRgb=\"On\"")
writeLines(paste(bedHeader, sep = "\n"), bedOut)

colorsList = c(brewer.pal(9, "Set1"), brewer.pal(9, "Pastel1"), brewer.pal(12, "Set3"))
bedDat = data.frame(chr=splitGet(rownames(dat), "\\.", 1), start=as.numeric(splitGet(rownames(dat), "\\.", 2)), end=as.numeric(splitGet(rownames(dat), "\\.", 3)), Cluster=dat$Cluster)
bedDat$Cluster = factor(bedDat$Cluster, levels = sort(unique(bedDat$Cluster)))
bedDat$Color = colorsList[ as.integer(bedDat$Cluster) ]
bedDat$Color2 = t(col2rgb(bedDat$Color))
bedDat2 = data.frame(bedDat$chr, bedDat$start, bedDat$end, bedDat$Cluster, 0, ".", bedDat$start, bedDat$end, paste(bedDat$Color2[,1], bedDat$Color2[,2], bedDat$Color2[,3], sep=","))
write.table(bedDat2, bedOut, sep = "\t", append = TRUE, col.names = FALSE, quote=FALSE, row.names=FALSE)

clustBreakSpec = cumsum(rle(datL5$Cluster)$lengths)


## Match with diff gene expression
{
    great = read.tsv2(greatDomainsFile)
    colnames(great) = c("chr", "start", "end", "Gene", "TSS", "strand")
    great$strand = NULL
    greatGR = makeGRangesFromDataFrame(great, ignore.strand = T, keep.extra.columns = T)

    tss = great[, c("Gene", "chr", "TSS")]
    
    #!! Find closest gene present in RNAseq. Can be removed depending on the goals.
    tss = tss[ tss$Gene %in% rna$Gene, ]
    tss$start = tss$TSS - 5000
    tss$end = tss$TSS + 5000
    tssGR = makeGRangesFromDataFrame(tss, ignore.strand = T, keep.extra.columns = T)

    tssExact = great[, c("Gene", "chr", "TSS")]
    #!! Find closest gene present in RNAseq. Can be removed depending on the goals.
    tssExact = tssExact[ tssExact$Gene %in% rna$Gene, ]
    tssExact$start = tssExact$TSS - 1
    tssExact$end = tssExact$TSS + 1
    tssGRExact = makeGRangesFromDataFrame(tssExact, ignore.strand = T, keep.extra.columns = T)


    datCoord = addRownameColumn(dat2, "PeakID")
    datCoord$chr = splitGet(datCoord$PeakID, "\\.", 1)
    datCoord$start = as.numeric(splitGet(datCoord$PeakID, "\\.", 2))
    datCoord$end = as.numeric(splitGet(datCoord$PeakID, "\\.", 3))
    datCoordGR = makeGRangesFromDataFrame(datCoord, ignore.strand = T, keep.extra.columns = T)



    ### GREAT
    {
        ov = findOverlaps(datCoordGR, greatGR, ignore.strand = T)
        ovDat = data.frame(datCoord[queryHits(ov), ], great[subjectHits(ov), c("Gene", "TSS")])
        datExp = merge(ovDat, rna, by.x="Gene", by.y="Gene")
        write.tsv(datExp, paste0(curOutName, "_MarkScoreCluster_vs_Expression_GREAT.txt"))



        ggplot(datExp, aes(x=Cluster, y=log2FoldChange, fill=Comparison)) + geom_hline(yintercept=0) + geom_violin() + stat_summary(fun.data = "median.quartile", geom = "pointrange", color = "white", alpha=0.8, size=1, fatten=2) + theme_cem + facet_grid(Comparison~.) + scale_fill_manual(values=colorsDiff) #+ coord_trans(ylim=c(-5,5))
        ggsave(paste0(curOutName, "_ExpByCluster_GREAT_1.png"), width=18, height=12)
        ggsave(paste0(curOutName, "_ExpByCluster_GREAT_1.pdf"), width=18, height=12)

        ggplot(datExp, aes(x=Comparison, y=log2FoldChange, fill=Comparison)) + geom_hline(yintercept=0) + geom_violin() + stat_summary(fun.data = "median.quartile", geom = "pointrange", color = "white", alpha=0.8, size=0.1) + theme_cem + facet_grid(~Cluster) + scale_fill_manual(values=colorsDiff) #+ coord_trans(ylim=c(-5,5))
        ggsave(paste0(curOutName, "_ExpByCluster_GREAT_2.png"), width=18, height=8)
        ggsave(paste0(curOutName, "_ExpByCluster_GREAT_2.pdf"), width=18, height=8)
    }
    ### Closest TSS
    {
        distToNearest = distanceToNearest(datCoordGR, tssGRExact, ignore.strand = T)
        ovDat = data.frame(datCoord[distToNearest@from, ], tssExact[distToNearest@to, c("Gene", "TSS")])
        ovDat$DistanceToNearestTSS = distToNearest@elementMetadata@listData[["distance"]]
        datExp = merge(ovDat, rna, by.x="Gene", by.y="Gene")             
        write.tsv(datExp, paste0(curOutName, "_MarkScoreCluster_vs_Expression_ClosestTSS.txt"))
        

        
        ggplot(datExp, aes(x=Cluster, y=log2FoldChange, fill=Comparison)) + geom_hline(yintercept=0) + geom_violin() + stat_summary(fun.data = "median.quartile", geom = "pointrange", color = "white", alpha=0.8, size=1, fatten=2) + theme_cem + facet_grid(Comparison~.) + scale_fill_manual(values=colorsDiff) #+ coord_trans(ylim=c(-5,5))
        ggsave(paste0(curOutName, "_ExpByCluster_ClosestTSS_1.png"), width=18, height=12)
        ggsave(paste0(curOutName, "_ExpByCluster_ClosestTSS_1.pdf"), width=18, height=12)
        
        ggplot(datExp, aes(x=Comparison, y=log2FoldChange, fill=Comparison)) + geom_hline(yintercept=0) + geom_violin() + stat_summary(fun.data = "median.quartile", geom = "pointrange", color = "white", alpha=0.8, size=0.1) + theme_cem + facet_grid(~Cluster) + scale_fill_manual(values=colorsDiff) #+ coord_trans(ylim=c(-5,5))
        ggsave(paste0(curOutName, "_ExpByClusterClosestTSS_2.png"), width=18, height=8)
        ggsave(paste0(curOutName, "_ExpByClusterClosestTSS_2.pdf"), width=18, height=8)
        
    }

    ### TSS +- X kb
    {
        ov = findOverlaps(datCoordGR, tssGR, ignore.strand = T)
        ovDat = data.frame(datCoord[queryHits(ov), ], tss[subjectHits(ov), c("Gene", "TSS")])
        datExpTss = merge(ovDat, rna, by.x="Gene", by.y="Gene")
        write.tsv(datExpTss, paste0(curOutName, "_MarkScoreCluster_vs_Expression_TSS5kb.txt"))


        ggplot(datExpTss, aes(x=Cluster, y=log2FoldChange, fill=Comparison)) + geom_hline(yintercept=0) + geom_violin() + stat_summary(fun.data = "median.quartile", geom = "pointrange", color = "white", alpha=0.8, size=1, fatten=2) + theme_cem + facet_grid(Comparison~.) + scale_fill_manual(values=colorsDiff) #+ coord_trans(ylim=c(-5,5))
        ggsave(paste0(curOutName, "_ExpByCluster_TSS5kb_1.png"), width=18, height=12)
        ggsave(paste0(curOutName, "_ExpByCluster_TSS5kb_1.pdf"), width=18, height=12)

        ggplot(datExpTss, aes(x=Comparison, y=log2FoldChange, fill=Comparison)) + geom_hline(yintercept=0) + geom_violin() + stat_summary(fun.data = "median.quartile", geom = "pointrange", color = "white", alpha=0.8, size=0.1) + theme_cem + facet_grid(~Cluster) + scale_fill_manual(values=colorsDiff) #+ coord_trans(ylim=c(-5,5))
        ggsave(paste0(curOutName, "_ExpByCluster_TSS5kb_2.png"), width=18, height=8)
        ggsave(paste0(curOutName, "_ExpByCluster_TSS5kb_2.pdf"), width=18, height=8)
    }
}
#




####################
##### Plotting #####
####################

# Continue with closest to TSS:

datExp = read.tsv2(paste0(curOutName, "_MarkScoreCluster_vs_Expression_ClosestTSS.txt"))

datExp$Comparison2 = paste0("RNA_", datExp$Comparison)
datExp = datExp[order(datExp$Comparison2), ]
datExpTmp1 = datExp[, ! colnames(datExp) %in% c("V1", "chr", "start", "end", "TSS", "Comparison", "baseMean", "padj")]
datExpTmp2 = datExp[, ! colnames(datExp) %in% c("V1", "chr", "start", "end", "TSS", "Comparison", "log2FoldChange", "padj")]

datExpSum = dcast(datExpTmp1, ...~Comparison2, value.var="log2FoldChange", fun.aggregate=median)
datExpSumBaseExp = dcast(datExpTmp2, ...~Comparison2, value.var="baseMean", fun.aggregate=median)
datExpSumBaseExp = datExpSumBaseExp[, c("PeakID", colnames(datExpSumBaseExp)[ncol(datExpSumBaseExp)])] 
colnames(datExpSumBaseExp) = c("PeakID", "RNA_BaseExpression")
datExpSumBaseExp$RNA_BaseExpression = log10(datExpSumBaseExp$RNA_BaseExpression + 1)

datExpSum = merge(datExpSum, datExpSumBaseExp, by="PeakID")
datExpSumM = melt(datExpSum, c("PeakID", "Gene", "Cluster", "X1", "X2", "DistanceToNearestTSS"), variable.name="SampleType", value.name="Score")
datExpSum2 = datExpSum
datExpSum2$chr = splitGet(datExpSum2$PeakID, "\\.", 1)
datExpSum2$start = as.numeric(splitGet(datExpSum2$PeakID, "\\.", 2))
datExpSum2$end = as.numeric(splitGet(datExpSum2$PeakID, "\\.", 3))
datExpSumGR = makeGRangesFromDataFrame(datExpSum2)


datExpAnn = data.frame(PeakID=datExpSum2$PeakID)

## Overlap with custom annotations
for(curAnnName in names(customAnnotations[[GenomeRef]]))
{
    curAnn = customAnnotations[[GenomeRef]][[curAnnName]]
    datExpAnn[, curAnnName] = "No annotation"
    for(curType in unique(curAnn$Type))
    {
        curAnnType = curAnn[curAnn$Type %in% curType, ]
        curGR = makeGRangesFromDataFrame(curAnnType)
        curOv = findOverlaps(datExpSumGR, curGR, ignore.strand=T)
        datExpAnn[queryHits(curOv), curAnnName] = curType 
    }
    annFreq = as.data.frame(sort(table(datExpAnn[,curAnnName]), decreasing = T))
    annFreq = annFreq[ annFreq$Var1 != "No annotation", ]
    annLevels = unique(c(as.character(annFreq$Var1), "No annotation"))
    datExpAnn[,curAnnName] = factor(datExpAnn[,curAnnName], levels = rev(annLevels))
}


## Overlap with combinatorial annotations
maxN = 9
for(curAnnName in names(customCombAnnotations[[GenomeRef]]))
{
    curAnnAll = data.frame(PeakID=datExpSum2$PeakID)

    curAnn = customCombAnnotations[[GenomeRef]][[curAnnName]]

    for(curType in unique(names(curAnn)))
    {
        curAnnType = curAnn[[curType]]
        curGR = makeGRangesFromDataFrame(curAnnType)
        curOv = findOverlaps(datExpSumGR, curGR, ignore.strand=T)
        
        curAnnAll[, curType] = 0
        curAnnAll[queryHits(curOv), curType] = 1
    }
    curAnnAll$Name = ""
    for(curCol in setdiff(colnames(curAnnAll), c("PeakID")))
    {
        curDat = ifelse(curAnnAll[,curCol]==1, paste0(" & ", curCol), "")
        curAnnAll$Name = paste0(curAnnAll$Name , curDat)
    }
    curAnnAll$Name = gsub("^ & ", "", curAnnAll$Name)
    curAnnAll$Name[ curAnnAll$Name == ""] = "No annotation"
    annFreq = as.data.frame(sort(table(curAnnAll$Name), decreasing = T))
    annFreq = annFreq[ annFreq$Var1 != "No annotation", ]
    annFreq$Type2 = as.character(annFreq$Var1)
    
    if(nrow(annFreq) > maxN + 1)
    {
        annFreq$Type2[(maxN+1):nrow(annFreq)] = "Other combinations"
    }
    annLevels = unique(c(annFreq$Type2, "No annotation"))
    curRes = mapvalues(curAnnAll$Name, annFreq$Var1, annFreq$Type2)
    curRes = factor(curRes, levels = rev(annLevels))
    datExpAnn[, curAnnName] = curRes
} 

datExpSum2 = merge(datExpSum2, datExpAnn, by="PeakID")



clusterMedian = datExpSum %>%
                    select(-PeakID,-Gene) %>%
					group_by(Cluster) %>%
					summarise_each(funs(median))

clusterMedian2 = setRemoveRownames(as.data.frame(clusterMedian), 1)
rownames(clusterMedian2) = paste0(rownames(clusterMedian2))
clusterMedian2 = clusterMedian2[, ! colnames(clusterMedian2) %in% c("X1", "X2","Cluster")]
breaks = c(0, seq(0.1, 0.9, length.out=49), 1)

histWtLimit = c(0.5, 6)
histDiffLimit = c(-2, 2)
rnaDiffLimit = c(-2, 2)
selFactorCols = grep(paste0("^(", paste0(factorsForClustering, collapse="|"), ")_"), colnames(clusterMedian2), value=T)

# Plot the cluster median heatmap
plotClusterHeatmap = function(sumDat, rawDat)
{

    plotlist = list()
    plotWidths = c()
    curDat = melt(addRownameColumn(sumDat, "Cluster"), "Cluster", variable.name="Variable", value.name="Score")
    clusters = sort(unique(curDat$Cluster))

    curColNames = colnames(sumDat)
    curColNames = curColNames[ 1:(grep("Distance", curColNames)[1]-1) ]
    
    colTypes = unique(splitGet(selFactorCols, "_", 2))

    for(curColType in colTypes)
    {
        curColumns = grep(paste0("_", curColType), colnames(sumDat), value=T)
        curColumns = curColumns[ curColumns %in% selFactorCols]
        curColor = if(grepl("vs", curColType)) colScaleDiff else colScale
        curLimits = if(grepl("vs", curColType)) histDiffLimit else histWtLimit
        
        curGG = ggplot(curDat[curDat$Variable %in% curColumns, ], aes(x=Variable, y=Cluster, fill=Score)) + geom_tile() + theme_cem + scale_fill_gradientn(colours=curColor, limits = curLimits, oob=scales::squish) + scale_y_discrete(limits=rev(clusters)) + guides(fill = guide_colourbar(label.theme = element_text(angle=-90)))
        if(curColType == colTypes[1]) # first element
        {
            curGG = curGG + theme(legend.position = "bottom", plot.margin = unit(c(0.1,0,0,0.1), "in"))
        }else
        {
            curGG = curGG + theme(legend.position = "bottom", axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0.1,0.1,0,0), "in"), axis.ticks.length.y = unit(0, "pt")) + labs(y = NULL)
        }
        
        plotlist[[ length(plotlist) + 1 ]] = curGG
        plotWidths = c(plotWidths, sqrt(length(curColumns)))
    }

   
    ggRnaHeat = ggplot(curDat[grepl(paste0("RNA_.*Vs", wtColumn), curDat$Variable), ], aes(x=Variable, y=Cluster, fill=Score)) + geom_tile() + theme_cem + scale_fill_gradientn(colours=colScaleDiffRNA, limits = rnaDiffLimit, oob=scales::squish) + theme(legend.position = "bottom", axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0.1,0.1,0,0), "in"), axis.ticks.length.y = unit(0, "pt")) + labs(y = NULL) + scale_y_discrete(limits=rev(clusters)) + guides(fill = guide_colourbar(label.theme = element_text(angle=-90)))
    
    rnaDat = gather(rawDat, "Comparison", "log2FoldChange", grep(paste0("RNA_.*Vs", wtColumn), colnames(sumDat), value=T))
    rnaDat$Comparison = gsub("RNA_", "", rnaDat$Comparison)
    rnaDat$Comparison = factor(rnaDat$Comparison, levels = rev(names(colorsDiff)))

    rnaP = rnaDat %>%
        group_by(Comparison, Cluster) %>%
        rstatix::wilcox_test(log2FoldChange ~ 1) %>%
        rstatix::adjust_pvalue() %>%
        mutate(y.position = 3, p.adj = signif(p.adj, 1))
    rnaP$p.level = "ns"
    rnaP$p.level[ rnaP$p.adj < 0.05] = "*"
    rnaP$p.level[ rnaP$p.adj < 0.01] = "**"
    rnaP$p.level[ rnaP$p.adj < 0.001] = "***"
    rnaP$p.level[ rnaP$p.adj < 0.0001] = "****"
    rnaP$p.level[ rnaP$p.adj < 0.00001] = "*****"

    ggRnaViolin = ggplot(rnaDat, aes(x=Cluster, y=log2FoldChange, fill=Comparison)) + geom_hline(yintercept=0) + geom_violin() + stat_summary(fun.data = "median.quartile", geom = "pointrange", color = "white", alpha=0.8, size=0.4, position=position_dodge(width=0.9)) + theme_cem + scale_fill_manual(values=colorsDiff) + coord_flip(ylim=c(-4,4)) + theme(legend.position = "bottom",  axis.title.y = element_blank(), plot.margin = unit(c(0.1,0.1,0,0), "in"), axis.ticks.length.y = unit(0, "pt")) + scale_x_discrete(limits=rev(clusters)) + guides(fill = guide_legend(ncol=1, title.theme = element_blank())) + ggpubr::stat_pvalue_manual(rnaP, label = "p.level", xmin = "Cluster", xmax = NULL,  position = position_dodge(0.95))

    ggDist = ggplot(rawDat, aes(x=Cluster, y=DistanceToNearestTSS+1)) + geom_violin(fill="lightblue") + stat_summary(fun.data = "median.quartile", geom = "pointrange", color = "black", alpha=0.8, size=1, fatten=2) + scale_y_log10(breaks=c(1,100, 1000,10000,100000,250000), labels=c("TSS", "100bp", "1kb", "10kb", "100kb", "250kb")) + theme_cem + ylab("Distance to TSS") + coord_flip()+ theme(legend.position = "bottom",  axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0.1,0.1,0,0), "in"), axis.ticks.length.y = unit(0, "pt")) + scale_x_discrete(limits=rev(clusters))
    
    rawDat$size = rawDat$end - rawDat$start
    ggSize = ggplot(rawDat, aes(x=Cluster, y=size)) + geom_violin(fill="moccasin") + stat_summary(fun.data = "median.quartile", geom = "pointrange", color = "black", alpha=0.8, size=1, fatten=2) + scale_y_log10(breaks=c(100, 1000,5000,10000,25000,50000,100000), labels=c("100bp", "1kb", "5kb", "10kb", "25kb", "50kb", "100kb")) + theme_cem + ylab("Size of region") + coord_flip()+ theme(legend.position = "bottom",  axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0.1,0.1,0,0), "in"), axis.ticks.length.y = unit(0, "pt")) + scale_x_discrete(limits=rev(clusters))
    
    ggNum = ggplot(rawDat, aes(x=Cluster)) + geom_bar() + theme_cem + ylab("Number of peaks") + coord_flip()+ theme(legend.position = "bottom",  axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0.1,0.1,0,0), "in"), axis.ticks.length.y = unit(0, "pt")) + scale_x_discrete(limits=rev(clusters))
    
    plotlist[[ length(plotlist) + 1 ]] = ggRnaHeat
    plotWidths = c(plotWidths, 1)
    plotlist[[ length(plotlist) + 1 ]] = ggRnaViolin
    plotWidths = c(plotWidths, 2)
    plotlist[[ length(plotlist) + 1 ]] = ggDist
    plotWidths = c(plotWidths, 2)
    plotlist[[ length(plotlist) + 1 ]] = ggSize
    plotWidths = c(plotWidths, 2)
    plotlist[[ length(plotlist) + 1 ]] = ggNum
    plotWidths = c(plotWidths, 2)
    
    for(curColType in colTypes)
    {
        curColumns = grep(paste0("_", curColType), colnames(sumDat), value=T)
        curColumns = curColumns[ ! curColumns %in% selFactorCols]
        if(length(curColumns) == 0) next
        curColor = if(grepl("vs", curColType)) colScaleDiff else colScale
        curLimits = if(grepl("vs", curColType)) histDiffLimit else histWtLimit
        
        curGG = ggplot(curDat[curDat$Variable %in% curColumns, ], aes(x=Variable, y=Cluster, fill=Score)) + geom_tile() + theme_cem + scale_fill_gradientn(colours=curColor, limits = curLimits, oob=scales::squish) + scale_y_discrete(limits=rev(clusters)) + guides(fill = guide_colourbar(label.theme = element_text(angle=-90)))

        curGG = curGG + theme(legend.position = "bottom", axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0.1,0.1,0,0), "in"), axis.ticks.length.y = unit(0, "pt")) + labs(y = NULL)
        
        plotlist[[ length(plotlist) + 1 ]] = curGG
        plotWidths = c(plotWidths, sqrt(length(curColumns)))
    }

    
    
    for(curAnn in setdiff(colnames(datExpAnn), "PeakID"))
    {
        if(curAnn %in% names(annotationColors))
        {
            rawDat[, curAnn] = factor(rawDat[, curAnn], levels = names(annotationColors[[curAnn]]))
        }
        ggAnn = ggplot(rawDat, aes_string(x="Cluster", fill=curAnn)) + geom_bar(position="fill") + theme_cem + scale_y_continuous(name=curAnn, labels=scales::percent) + coord_flip()+ theme(legend.position = "bottom",  axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0.1,0.1,0,0), "in"), axis.ticks.length.y = unit(0, "pt")) + scale_x_discrete(limits=rev(clusters)) + guides(fill = guide_legend(ncol=1, title.theme = element_blank())) 
        
        if(curAnn %in% names(annotationColors))
        {
            ggAnn = ggAnn + scale_fill_manual(values=annotationColors[[curAnn]])
        }else
        {
            ggAnn = ggAnn + scale_fill_manual(values=iwanthue(length(unique(rawDat[, curAnn]) )))
        }
        
        plotlist[[ length(plotlist) + 1 ]] = ggAnn
        plotWidths = c(plotWidths, 2)
    }
    plot_grid(plotlist = plotlist, align = 'hv', axis="tb", ncol = length(plotlist), rel_widths = plotWidths)
}
ggPlot = plotClusterHeatmap(clusterMedian2, datExpSum2)
ggPlot
ggsave(paste0(curOutName, "_cluster_heatmap_All.png"), ggPlot, width=45, height=16, dpi=150)
ggsave(paste0(curOutName, "_cluster_heatmap_All.pdf"), ggPlot, width=45, height=16, dpi=150)


# Plot each cluster as a heatmap
datExpSumOrdered = datExpSum2
for(curClust in sort(unique(datExpSumOrdered$Cluster)))
{
    curWide = datExpSumOrdered[ datExpSumOrdered$Cluster %in% curClust, ]
    
    
    curColNames = colnames(curWide)
    curColNames = curColNames[ 1:(grep("chr", curColNames)[1]-1) ]

    # Order/cluster peaks
    curWideDat = setRemoveRownames(curWide, 1)
    curWideDat = curWideDat[, ! colnames(curWideDat) %in% c("PeakID", "Gene", "X1", "X2", "Cluster", "DistanceToNearestTSS", "chr", "start", "end")]
    curWideDat = curWideDat[, colnames(curWideDat) %in% selFactorCols | grepl(paste0("RNA.*Vs", wtColumn), colnames(curWideDat))]
    dists = dist(curWideDat)
    if(T)
    {
        hc = hclust(dists, method="ward.D")
        datOrd = hc$order
    }else
    {
        datOrd = seriation::seriate(dists)
        datOrd = get_order(datOrd)
    }
    curWideDat = curWideDat[datOrd, ]
    
    curWide$PeakID = factor(curWide$PeakID, levels = rownames(curWideDat))
    curWide = curWide[ order(curWide$PeakID), ]
    
    curDat = melt(curWide[, curColNamesExtra ], c("PeakID", "Gene", "Cluster"), variable.name="Variable", value.name="Score")
    
    colTypes = unique(splitGet(selFactorCols, "_", 2))
    
    plotlist = list()
    plotWidths = c()
    plotlistV = list()
    plotWidthsV = c()
    for(curColType in colTypes)
    {
        curColumns = grep(paste0("_", curColType), colnames(curWide), value=T)
        curColumns = curColumns[ curColumns %in% selFactorCols]
        curColor = if(grepl("vs", curColType)) colScaleDiff else colScale
        curLimits = if(grepl("vs", curColType)) histDiffLimit else histWtLimit
        
        curGG = ggplot(curDat[curDat$Variable %in% curColumns, ], aes(x=Variable, y=PeakID, fill=Score)) + geom_tile() + theme_cem + scale_fill_gradientn(colours=curColor, limits = curLimits, oob=scales::squish) + guides(fill = guide_colourbar(label.theme = element_text(angle=-90))) + theme(legend.position = "bottom", axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0.1,0.1,0,0), "in"), axis.ticks.length.y = unit(0, "pt")) + labs(y = NULL)
        if(curColType == colTypes[1]) # first element
        {
            curGG = curGG + theme(legend.position = "bottom", plot.margin = unit(c(0.1,0,0,0.1), "in")) + ggtitle(curClust)
        }else
        {
            curGG = curGG + theme(legend.position = "bottom", axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0.1,0.1,0,0), "in"), axis.ticks.length.y = unit(0, "pt")) + labs(y = NULL)
        }
        
        plotlist[[ length(plotlist) + 1 ]] = curGG
        plotWidths = c(plotWidths, sqrt(length(curColumns)))
        
        curGG = ggplot(curDat[curDat$Variable %in% curColumns, ], aes(x=Variable, y=Score, fill=splitGet(Variable,"_",1))) + geom_violin() + stat_summary(fun.data = "median.quartile", geom = "pointrange", color = "white", alpha=0.8, size=1, fatten=2) + scale_fill_manual(values=colorsDef[["Factor"]]) + theme_cem + geom_hline(yintercept=0) + expand_limits(y=curLimits) + theme(legend.position="off")
        
        if(curColType == colTypes[1])
        {
            curGG = curGG + ggtitle(curClust, curColType)
        }else
        {
            curGG = curGG + ggtitle("", curColType) + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0.1,0.1,0,0), "in"), axis.ticks.length.y = unit(0, "pt")) + labs(y = NULL)
        }
        
        plotlistV[[ length(plotlistV) + 1 ]] = curGG
        plotWidthsV = c(plotWidthsV, sqrt(length(curColumns)))
    }
    
    ggRnaHeat = ggplot(curDat[grepl(paste0("RNA_.*Vs", wtColumn), curDat$Variable), ], aes(x=Variable, y=PeakID, fill=Score)) + geom_tile() + theme_cem + scale_fill_gradientn(colours=colScaleDiffRNA, limits = rnaDiffLimit, oob=scales::squish) + theme(legend.position = "bottom", axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0.1,0.1,0,0), "in"), axis.ticks.length.y = unit(0, "pt")) + labs(y = NULL) + guides(fill = guide_colourbar(label.theme = element_text(angle=-90)))
    plotlist[[ length(plotlist) + 1 ]] = ggRnaHeat
    plotWidths = c(plotWidths, 1)
    
    ggRnaViolin = ggplot(curDat[grepl(paste0("RNA_.*Vs", wtColumn), curDat$Variable), ], aes(x=Variable, y=Score, fill=splitGet(Variable,"_",2))) + geom_violin() + stat_summary(fun.data = "median.quartile", geom = "pointrange", color = "white", alpha=0.8, size=1, fatten=2) + scale_fill_manual(values=colorsDiff) + ggtitle("", "RNAseq") + theme_cem + geom_hline(yintercept=0) + expand_limits(y=rnaDiffLimit) + theme(legend.position="off") + theme( axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0.1,0.1,0,0), "in"), axis.ticks.length.y = unit(0, "pt")) + labs(y = NULL)
    plotlistV[[ length(plotlistV) + 1 ]] = ggRnaViolin
    plotWidthsV = c(plotWidthsV, 1)


    
    for(curColType in colTypes)
    {
        curColumns = grep(paste0("_", curColType), colnames(curWide), value=T)
        curColumns = curColumns[ ! curColumns %in% selFactorCols]
        if(length(curColumns) == 0) next
        curColor = if(grepl("vs", curColType)) colScaleDiff else colScale
        curLimits = if(grepl("vs", curColType)) histDiffLimit else histWtLimit
        
        curGG = ggplot(curDat[curDat$Variable %in% curColumns, ], aes(x=Variable, y=PeakID, fill=Score)) + geom_tile() + theme_cem + scale_fill_gradientn(colours=curColor, limits = curLimits, oob=scales::squish)

        curGG = curGG + theme(legend.position = "bottom", axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0.1,0.1,0,0), "in"), axis.ticks.length.y = unit(0, "pt")) + labs(y = NULL)
        
        plotlist[[ length(plotlist) + 1 ]] = curGG
        plotWidths = c(plotWidths, sqrt(length(curColumns)))
        
        curGG = ggplot(curDat[curDat$Variable %in% curColumns, ], aes(x=Variable, y=Score, fill=splitGet(Variable,"_",1))) + geom_violin() + stat_summary(fun.data = "median.quartile", geom = "pointrange", color = "white", alpha=0.8, size=1, fatten=2) + scale_fill_manual(values=colorsDef[["Factor"]]) + ggtitle("", curColType) + theme_cem + geom_hline(yintercept=0) + expand_limits(y=curLimits) + theme(legend.position="off") + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0.1,0.1,0,0), "in"), axis.ticks.length.y = unit(0, "pt")) + labs(y = NULL)
        plotlistV[[ length(plotlistV) + 1 ]] = curGG
        plotWidthsV = c(plotWidthsV, sqrt(length(curColumns)))
    }

    ggPlot = plot_grid(plotlist = plotlist, align = 'hv', axis="tb", ncol = length(plotlist), rel_widths = plotWidths)
    png(paste0(curOutName, "_cluster_heatmap_", curClust, ".png"), width=22, height=10, res=150, units="in")
    print(ggPlot)
    dev.off()
    pdf(paste0(curOutName, "_cluster_heatmap_", curClust, ".pdf"), width=22, height=10)
    print(ggPlot)
    dev.off()
    
    ggPlot2 = plot_grid(plotlist = plotlistV, align = 'hv', axis="tb", ncol = length(plotlistV), rel_widths = plotWidthsV)
    png(paste0(curOutName, "_cluster_violin_", curClust, ".png"), width=22, height=10, res=150, units="in")
    print(ggPlot2)
    dev.off()
    pdf(paste0(curOutName, "_cluster_violin_", curClust, ".pdf"), width=22, height=10)
    print(ggPlot2)
    dev.off()
}



## Plot signal on tSNE

ggPlotList = datExpSumM %>% 
  group_split(SampleType) %>% 
  map(
    function(x)
    {
        curType = x$SampleType[1]
        curPlot = NULL
        if(grepl(paste0("H3K4me1.*vs", wtColumn), curType))
        {
            ggplot(x, aes(X1, X2, color = Score)) + geom_point() + theme_cem + scale_color_gradientn(colours=colScaleDiff3, limits=histDiffLimit, oob = scales::squish)   + theme(legend.position = "bottom") + ggtitle(curType)
        }else if(grepl(paste0("H3K4me3.*vs", wtColumn), curType))         
        {                                   
            ggplot(x, aes(X1, X2, color = Score)) + geom_point() + theme_cem + scale_color_gradientn(colours=colScaleDiff2, limits=histDiffLimit, oob = scales::squish)   + theme(legend.position = "bottom") + ggtitle(curType)
        }else 
        if(grepl(paste0(".*vs", wtColumn), curType))      
        {                                          
            ggplot(x, aes(X1, X2, color = Score)) + geom_point() + theme_cem + scale_color_gradientn(colours=colScaleDiff, limits=histDiffLimit, oob = scales::squish)    + theme(legend.position = "bottom") + ggtitle(curType)
        }else if(grepl(paste0("_", wtColumn), curType))                   
        {                                                
           ggplot(x, aes(X1, X2, color = Score)) + geom_point() + theme_cem + scale_color_gradientn(colours=colScale, limits=histWtLimit, oob = scales::squish)          + theme(legend.position = "bottom")  + ggtitle(curType)
        }
		
        else if(grepl(paste0("RNA_.*Vs", wtColumn), curType))             
        {                                              
           ggplot(x, aes(X1, X2, color = Score)) + geom_point() + theme_cem + scale_color_gradientn(colours=colScaleDiffRNA, limits=rnaDiffLimit, oob = scales::squish)  + theme(legend.position = "bottom")  + ggtitle(curType)
        }
		
        else if(grepl("RNA_BaseExpression", curType))       
        {                                                
           ggplot(x, aes(X1, X2, color = Score)) + geom_point() + theme_cem + scale_color_gradientn(colours=colScale, limits=c(0.5,4.5), oob = scales::squish)           + theme(legend.position = "bottom")  + ggtitle(curType)
        }
    }
  ) 

ggPlot = plot_grid(plotlist = ggPlotList, align = 'hv', ncol = 4)
ggPlot
ggsave(paste0(curOutName, "_tSNE_scaled_color_ClosestToTSS.png"), ggPlot, width=22, height=20, dpi=150)
ggsave(paste0(curOutName, "_tSNE_scaled_color_ClosestToTSS.pdf"), ggPlot, width=22, height=20, dpi=150)

ggPlotList2 = AugmentPlotList(ggPlotList, 22/5, 20/5)
ggPlot2 = plot_grid(plotlist = ggPlotList2, align = 'hv', ncol = 4)
ggsave(paste0(curOutName, "_tSNE_scaled_color_ClosestToTSS_aug.png"), ggPlot2, width=22, height=20, dpi=150)
ggsave(paste0(curOutName, "_tSNE_scaled_color_ClosestToTSS_aug.pdf"), ggPlot2, width=22, height=20, dpi=150)


## Plot signal on tSNE, smoothed

ggPlotList = datExpSumM %>% 
  group_split(SampleType) %>% 
  map(
    function(x)
    {
        curType = x$SampleType[1]
        curPlot = NULL
        sm = smoothScore2d(x$Score, x$X1, x$X2, numGrid=100, knn=50, m=2)
        sm$SampleType = curType
        if(grepl(paste0("H3K4me1.*vs", wtColumn), curType))
        {
            ggplot(sm, aes(x, y, fill = score)) + geom_tile() + theme_cem + scale_fill_gradientn(colours=colScaleDiff3, limits=histDiffLimit, oob = scales::squish)   + theme(legend.position = "bottom") + ggtitle(curType)
        }else if(grepl(paste0("H3K4me3.*vs", wtColumn), curType))      
        {                                        
            ggplot(sm, aes(x, y, fill = score)) + geom_tile() + theme_cem + scale_fill_gradientn(colours=colScaleDiff2, limits=histDiffLimit, oob = scales::squish)   + theme(legend.position = "bottom") + ggtitle(curType)
        }else 
        if(grepl(paste0(".*vs", wtColumn), curType))    
        {                                          
            ggplot(sm, aes(x, y, fill = score)) + geom_tile() + theme_cem + scale_fill_gradientn(colours=colScaleDiff, limits=histDiffLimit, oob = scales::squish)    + theme(legend.position = "bottom") + ggtitle(curType)
        }else if(grepl(paste0("_", wtColumn), curType))            
        {                                            
           ggplot(sm, aes(x, y, fill = score)) + geom_tile()  + theme_cem + scale_fill_gradientn(colours=colScale, limits=histWtLimit, oob = scales::squish)          + theme(legend.position = "bottom") + ggtitle(curType)
        }                                            
        else if(grepl(paste0("RNA_.*Vs", wtColumn), curType))   
        {                                          
           ggplot(sm, aes(x, y, fill = score)) + geom_tile()  + theme_cem + scale_fill_gradientn(colours=colScaleDiffRNA, limits=rnaDiffLimit, oob = scales::squish)  + theme(legend.position = "bottom") + ggtitle(curType)
        }                                                                        
        else if(grepl("RNA_BaseExpression", curType))                          
        {                                                                       
           ggplot(sm, aes(x, y, fill = score)) + geom_tile()  + theme_cem + facet_wrap(~SampleType) + scale_fill_gradientn(colours=colScale, limits=c(0.5,4.5), oob = scales::squish)           + theme(legend.position = "bottom") + ggtitle(curType)
        }
    }
  ) 

ggPlot = plot_grid(plotlist = ggPlotList, align = 'hv', ncol = 4)
ggPlot
ggsave(paste0(curOutName, "_tSNEsmooth_scaled_color_ClosestToTSS.png"), ggPlot, width=22, height=20, dpi=150)
ggsave(paste0(curOutName, "_tSNEsmooth_scaled_color_ClosestToTSS.pdf"), ggPlot, width=22, height=20, dpi=150)

ggPlotList2 = AugmentPlotList(ggPlotList, 22/5, 20/5)
ggPlot2 = plot_grid(plotlist = ggPlotList2, align = 'hv', ncol = 4)
ggsave(paste0(curOutName, "_tSNEsmooth_scaled_color_ClosestToTSS_aug.png"), ggPlot2, width=22, height=20, dpi=150)
ggsave(paste0(curOutName, "_tSNEsmooth_scaled_color_ClosestToTSS_aug.pdf"), ggPlot2, width=22, height=20, dpi=150)

