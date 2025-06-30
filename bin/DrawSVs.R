#!/usr/bin/env Rscript

## Heavily modified from ARRIBA's original script draw_fusions.R
## By blancoj@mskcc.org

# Print warnings as they happen instead of collecting them for after a loop ends
options(warn = 1)

# Define valid parameters
parameters <- list(
  SVs = list("SVsFile", "file", "SVs.tsv", TRUE),
  annotation = list("exonsFile", "file", "annotation.gtf", TRUE),
  output = list("outputFile", "string", "output.pdf", TRUE),
  alignments = list("alignmentsFile", "file", "Aligned.sortedByCoord.out.bam"),
  cytobands = list("cytobandsFile", "file", "cytobands.tsv"),
  minConfidenceForCircosPlot = list("minConfidenceForCircosPlot", "string", "Medium"), # nolint
  proteinDomains = list("proteinDomainsFile", "file", "protein_domains.gff3"),
  sampleName = list("sampleName", "string", ""),
  squishIntrons = list("squishIntrons", "bool", TRUE),
  printExonLabels = list("printExonLabels", "bool", TRUE),
  render3dEffect = list("render3dEffect", "bool", TRUE),
  pdfWidth = list("pdfWidth", "numeric", 11.692),
  pdfHeight = list("pdfHeight", "numeric", 8.267),
  color1 = list("color1", "string", "#FA7550"),
  color2 = list("color2", "string", "#9683F5"),
  mergeDomainsOverlappingBy = list("mergeDomainsOverlappingBy", "numeric", 0.9),
  optimizeDomainColors = list("optimizeDomainColors", "bool", FALSE),
  fontSize = list("fontSize", "numeric", 1),
  fontFamily = list("fontFamily", "string", "Helvetica"),
  showIntergenicVicinity = list("showIntergenicVicinity", "string", "0"),
  transcriptSelection = list("transcriptSelection", "string", "coverage"),
  fixedScale = list("fixedScale", "numeric", 0),
  coverageRange = list("coverageRange", "string", "0")
)

# Print help if necessary

args <- commandArgs(trailingOnly = TRUE)
if (any(grepl("^--help", args)) || length(args) == 0) {
  usage <- "Usage: DrawSVs.R"
  for (parameter in names(parameters)) {
    usage <- paste0(usage, " ")
    if (length(parameters[[parameter]]) <= 3 || !parameters[[parameter]][[4]])
      usage <- paste0(usage, "[")
    usage <- paste0(usage, "--", parameter, "=", parameters[[parameter]][[3]])
    if (length(parameters[[parameter]]) <= 3 || !parameters[[parameter]][[4]])
      usage <- paste0(usage, "]")
  }
  message(usage)
  quit("no", ifelse(length(args) == 0, 1, 0))
}

# Make sure mandatory arguments are present

for (parameter in names(parameters)) {
  if (length(parameters[[parameter]]) > 3 && parameters[[parameter]][[4]])
    if (!any(grepl(paste0("^--", parameter, "="), args), perl = TRUE))
      stop(paste0("Missing mandatory argument: --", parameter))
}
# Set default values

for (parameter in names(parameters)) {
  assign(parameters[[parameter]][[1]], ifelse(parameters[[parameter]][[2]] == "file", "", parameters[[parameter]][[3]])) # nolint
}
# Parse command-line parameters

for (arg in args) {
	argName <- sub("=.*", "", sub("^--", "", arg, perl = TRUE), perl = TRUE)
	argValue <- sub("^[^=]*=", "", arg, perl = TRUE)
	if (!(argName %in% names(parameters)) || !grepl("^--", arg, perl = TRUE))
		stop(paste("Unknown parameter:", arg))
	if (parameters[[argName]][[2]] == "bool") {
		if (argValue %in% c("TRUE", "T", "FALSE", "F")) {
			assign(parameters[[argName]][[1]], as.logical(argValue))
		} else {
			stop(paste0("Invalid argument to --", argName))
		}
	} else if (parameters[[argName]][[2]] == "string") {
		assign(parameters[[argName]][[1]], argValue)
	} else if (parameters[[argName]][[2]] == "numeric") {
		if (is.na(suppressWarnings(as.numeric(argValue))))
			stop(paste0("Invalid argument to --", argName))
		assign(parameters[[argName]][[1]], as.numeric(argValue))
	} else if (parameters[[argName]][[2]] == "file") {
		if (file.access(argValue) == -1)
			stop(paste("Cannot read file:", argValue))
		assign(parameters[[argName]][[1]], argValue)
	}
}

# Validate values of parameters

if (cytobandsFile == "")
	warning("Missing parameter '--cytobands'. No ideograms nor circos plots can be drawn.")
if (!(minConfidenceForCircosPlot %in% c("None", "Low", "Medium", "High")))
	stop("Invalid argument to --minConfidenceForCircosPlot")
showIntergenicVicinity <- as.list(unlist(strsplit(showIntergenicVicinity, ",", fixed = TRUE)))
if (!(length(showIntergenicVicinity) %in% c(1,4)))
	stop(paste0("Invalid argument to --showIntergenicVicinity"))
showIntergenicVicinity <- lapply(showIntergenicVicinity, function(x) {
	if (x == "closestGene") {
		return("exon")
	} else if (x == "closestProteinCodingGene") {
		return("CDS")
	} else if (is.na(suppressWarnings(as.numeric(x))) || as.numeric(x) < 0) {
		stop(paste0("Invalid argument to --showIntergenicVicinity"))
	} else {
		return(as.numeric(x))
	}
})
if (length(showIntergenicVicinity) == 1)
	showIntergenicVicinity <- rep(showIntergenicVicinity, 4)
if (squishIntrons)
	if (any(!is.numeric(unlist(showIntergenicVicinity))) || any(showIntergenicVicinity > 0))
		stop("--squishIntrons must be disabled, when --showIntergenicVicinity is > 0")
if (!(transcriptSelection %in% c("coverage", "provided", "canonical")))
	stop("Invalid argument to --transcriptSelection")
if (fixedScale < 0)
	stop("Invalid argument to --fixedScale")
if (!(fontFamily %in% names(pdfFonts())))
	stop(paste0("Unknown font: ", fontFamily, ". Available fonts: ", paste(names(pdfFonts()), collapse=", ")))
coverageRange <- suppressWarnings(as.numeric(unlist(strsplit(coverageRange, ",", fixed = TRUE))))
if (!(length(coverageRange) %in% 1:2) || any(is.na(coverageRange)) || any(coverageRange < 0))
	stop("Invalid argument to --coverageRange")

# Check if required packages are installed

if (!suppressPackageStartupMessages(require(GenomicRanges)))
	warning("Package 'GenomicRanges' is not installed. Circos plots cant be drawn.")
if (!suppressPackageStartupMessages(require(circlize)))
	warning("Package 'circlize' is not installed. Circos plots cant be drawn.")
if (alignmentsFile != "")
	if (!suppressPackageStartupMessages(require(GenomicAlignments)))
		stop("Package 'GenomicAlignments' must be installed when '--alignments' is used")

# Define colors

changeColorBrightness <- function(color, delta) {
	rgb(
		min(255,max(0,col2rgb(color)["red",]+delta)),
		min(255,max(0,col2rgb(color)["green",]+delta)),
		min(255,max(0,col2rgb(color)["blue",]+delta)),
		maxColorValue=255
	)
}
getDarkColor <- function(color) { changeColorBrightness(color, -100) }
getBrightColor <- function(color) { changeColorBrightness(color, +100) }
darkColor1 <- getDarkColor(color1)
darkColor2 <- getDarkColor(color2)
circosColors <- c(Translocation="#000000",Duplication="#00DB1E", Deletion="#FF0B46", Inversion="#1161FF")

# Convenience functions to add/remove "chr" prefix

addChr <- function(contig) {
	ifelse(contig == "MT", "chrM", paste0("chr", contig))
}
removeChr <- function(contig) {
	sub("^chr", "", sub("^chrM", "MT", contig, perl = TRUE), perl = TRUE)
}

# Convenience function to check if a value is between two others

between <- function(value, start, end) {
	value >= start & value <= end
}

# Read SVs

SVs <- read.table(SVsFile, stringsAsFactors=F, sep="\t", header=T, comment.char="", quote="")
if (colnames(SVs)[1] == "SV_BC") { # Survivor
	colnames(SVs)[colnames(SVs) %in% c("strand1.gene.SV.", "strand2.gene.SV.")] <- c("strand1", "strand2")
	SVs$display_contig1 <- sub(":[^:]*$", "", SVs$Breakpoint1, perl = TRUE)
	SVs$display_contig2 <- sub(":[^:]*$", "", SVs$Breakpoint2, perl = TRUE)
	SVs$contig1 <- removeChr(SVs$display_contig1)
	SVs$contig2 <- removeChr(SVs$display_contig2)
	SVs$Breakpoint1 <- as.numeric(sub(".*:", "", SVs$Breakpoint1, perl = TRUE))
	SVs$Breakpoint2 <- as.numeric(sub(".*:", "", SVs$Breakpoint2, perl = TRUE))
	SVs$split_reads1 <- SVs$split_reads1
	SVs$split_reads2 <- SVs$split_reads2
	SVs$type <- sub(".*(Translocation|Duplication|Deletion|Inversion).*", "\\1", SVs$type, perl = TRUE)
	SVs$SV_transcript <- gsub("[()^$]", "", SVs$SV_transcript)
} else if (colnames(SVs)[1] == "X.SVName") { # STAR-SV
	SVs$gene1 <- sub("\\^.*", "", SVs$LeftGene, perl = TRUE)
	SVs$gene2 <- sub("\\^.*", "", SVs$RightGene, perl = TRUE)
	SVs$strand1 <- sub(".*:(.)$", "\\1/\\1", SVs$LeftBreakpoint, perl = TRUE)
	SVs$strand2 <- sub(".*:(.)$", "\\1/\\1", SVs$RightBreakpoint, perl = TRUE)
	SVs$display_contig1 <- sub(":[^:]*:[^:]*$", "", SVs$LeftBreakpoint, perl = TRUE)
	SVs$display_contig2 <- sub(":[^:]*:[^:]*$", "", SVs$RightBreakpoint, perl = TRUE)
	SVs$contig1 <- removeChr(SVs$display_contig1)
	SVs$contig2 <- removeChr(SVs$display_contig2)
	SVs$Breakpoint1 <- as.numeric(sub(".*:([^:]*):[^:]*$", "\\1", SVs$LeftBreakpoint, perl = TRUE))
	SVs$Breakpoint2 <- as.numeric(sub(".*:([^:]*):[^:]*$", "\\1", SVs$RightBreakpoint, perl = TRUE))
	SVs$direction1 <- ifelse(grepl(":\\+$", SVs$LeftBreakpoint, perl = TRUE), "downstream", "upstream")
	SVs$direction2 <- ifelse(grepl(":\\+$", SVs$RightBreakpoint, perl = TRUE), "upstream", "downstream")
	SVs$gene_id1 <- sub(".*\\^", "", SVs$LeftGene, perl = TRUE)
	SVs$gene_id2 <- sub(".*\\^", "", SVs$RightGene, perl = TRUE)
	SVs$transcript_id1 <- ifelse(rep(!("CDS_LEFT_ID" %in% colnames(SVs)), nrow(SVs)), ".", SVs$CDS_LEFT_ID)
	SVs$transcript_id2 <- ifelse(rep(!("CDS_RIGHT_ID" %in% colnames(SVs)), nrow(SVs)), ".", SVs$CDS_RIGHT_ID)
	SVs$SV_transcript <- ifelse(rep(!("SV_CDS" %in% colnames(SVs)), nrow(SVs)), ".", toupper(sub("([a-z]*)", "\\1|", SVs$SV_CDS, perl = TRUE)))
	SVs$reading_frame <- ifelse(rep(!("PROT_SV_TYPE" %in% colnames(SVs)), nrow(SVs)), ".", ifelse(SVs$PROT_SV_TYPE == "INFRAME", "in-frame", ifelse(SVs$PROT_SV_TYPE == "FRAMESHIFT", "out-of-frame", ".")))
	SVs$split_reads <- SVs$JunctionReadCount
	SVs$discordant_mates <- SVs$SpanningFragCount
	SVs$site1 <- rep("exon", nrow(SVs))
	SVs$site2 <- rep("exon", nrow(SVs))
	SVs$confidence <- rep("High", nrow(SVs))
	SVs$type <- ifelse(SVs$contig1 != SVs$contig2, "Translocation", ifelse(SVs$direction1 == SVs$direction2, "Inversion", ifelse((SVs$direction1 == "downstream") == (SVs$Breakpoint1 < SVs$Breakpoint2), "Deletion", "Duplication")))
} else {
	stop("Unrecognized SV file format")
}

pdf(outputFile, onefile=T, width=pdfWidth, height=pdfHeight, title=ifelse(sampleName != "", sampleName, SVsFile))
par(family=fontFamily)

if (nrow(SVs) == 0) {
	plot(0, 0, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
	text(0, 0, "empty input file")
	warning("empty input file")
	dev.off()
	quit("no")
}

# Read cytoband annotation

cytobands <- NULL
if (cytobandsFile != "") {
	cytobands <- read.table(cytobandsFile, header=T, colClasses=c("character", "numeric", "numeric", "character", "character"))
	cytobands <- cytobands[order(cytobands$contig, cytobands$start, cytobands$end),]
}

# Read exon annotation

message("Loading annotation")
exons <- scan(exonsFile, what=list(contig="",src="",type="",start=0,end=0,score="",strand="",frame="",attributes=""), sep="\t", comment.char="#", quote='"', multi.line=F)
attr(exons, "row.names") <- .set_row_names(length(exons[[1]]))
class(exons) <- "data.frame"
exons <- exons[exons$type %in% c("exon","CDS"),c("contig","type","start","end","strand","attributes")]
exons$contig <- removeChr(exons$contig)
parseGtfAttribute <- function(attribute, gtf) {
	parsed <- sub(paste0(".*", attribute, "[ =]([^;]+).*"), "\\1", gtf$attributes, perl = TRUE)
	failedToParse <- parsed == gtf$attributes
	if (any(failedToParse)) {
		warning(paste0("Failed to parse '", attribute, "' attribute of ", sum(failedToParse), " record(s)."))
		parsed <- ifelse(failedToParse, "", parsed)
	}
	return(parsed)
}
exons$geneID <- parseGtfAttribute("gene_id", exons)
exons$geneName <- parseGtfAttribute("gene_name", exons)
exons$geneName <- ifelse(exons$geneName == "", exons$geneID, exons$geneName)
exons$transcript <- parseGtfAttribute("transcript_id", exons)
exons$exonNumber <- ifelse(rep(printExonLabels, nrow(exons)), parseGtfAttribute("exon_number", exons), "")

# Read protein domain annotation

proteinDomains <- NULL
if (proteinDomainsFile != "") {
	message("Loading protein domains")
	proteinDomains <- scan(proteinDomainsFile, what=list(contig="",src="",type="",start=0,end=0,score="",strand="",frame="",attributes=""), sep="\t", comment.char="", quote="", multi.line=F)
	attr(proteinDomains, "row.names") <- .set_row_names(length(proteinDomains[[1]]))
	class(proteinDomains) <- "data.frame"
	proteinDomains$color <- parseGtfAttribute("color", proteinDomains)
	proteinDomains$proteinDomainName <- sapply(parseGtfAttribute("Name", proteinDomains), URLdecode)
	proteinDomains$proteinDomainID <- parseGtfAttribute("protein_domain_id", proteinDomains)
}

# insert dummy annotations for intergenic breakpoints
if (any(SVs$site1 == "intergenic" | SVs$site2 == "intergenic")) {
	intergenicBreakpoints <- rbind(
		setNames(SVs[SVs$site1 == "intergenic",c("gene1", "strand1", "contig1", "Breakpoint1")], c("gene", "strand", "contig", "breakpoint")),
		setNames(SVs[SVs$site2 == "intergenic",c("gene2", "strand2", "contig2", "Breakpoint2")], c("gene", "strand", "contig", "breakpoint"))
	)
	exons <- rbind(exons, data.frame(
		contig=intergenicBreakpoints$contig,
		type="intergenic",
		start=sapply(intergenicBreakpoints$breakpoint-1000, max, 1),
		end=intergenicBreakpoints$breakpoint+1000,
		strand=".",
		attributes="",
		geneName=intergenicBreakpoints$gene,
		geneID=paste0(intergenicBreakpoints$contig, ":", intergenicBreakpoints$breakpoint),
		transcript=paste0(intergenicBreakpoints$contig, ":", intergenicBreakpoints$breakpoint),
		exonNumber="intergenic"
	))
	SVs[SVs$site1 == "intergenic","gene_id1"] <- paste0(SVs[SVs$site1 == "intergenic","contig1"], ":", SVs[SVs$site1 == "intergenic","Breakpoint1"])
	SVs[SVs$site2 == "intergenic","gene_id2"] <- paste0(SVs[SVs$site2 == "intergenic","contig2"], ":", SVs[SVs$site2 == "intergenic","Breakpoint2"])
}

drawVerticalGradient <- function(left, right, y, color, selection=NULL) {
	# Check if gradient should only be drawn in part of the region
	if (!is.null(selection)) {
		y <- y[selection]
		left <- left[selection]
		right <- right[selection]
	}
	# Draw gradient
	for (i in 1:length(y)) {
		polygon(
			c(left[1:i], right[1:i]),
			c(y[1:i], y[i:1]),
			border=NA,
			col=rgb(col2rgb(color)["red",], col2rgb(color)["green",], col2rgb(color)["blue",], col2rgb(color, alpha=T)["alpha",]*(1/length(y)), max=255)
		)
	}
}
if (!render3dEffect) # nullify function, if no 3D effect should be drawn
	drawVerticalGradient <- function(left, right, y, color, selection=NULL) { }

drawCurlyBrace <- function(left, right, top, bottom, tip) {
	smoothness <- 20
	x <- cumsum(exp(-seq(-2.5, 2.5, len=smoothness)^2))
	x <- x/max(x)
	y <- seq(top, bottom, len=smoothness)
	lines(left+(tip-left)+x*(left-tip), y)
	lines(tip+x*(right-tip), y)
}

drawIdeogram <- function(adjust, left, right, y, cytobands, contig, breakpoint) {
	# Define design of ideogram
	bandColors <- setNames(rgb(100:0, 100:0, 100:0, maxColorValue=100), paste0("gpos", 0:100))
	bandColors <- c(bandColors, gneg="#ffffff", acen="#ec4f4f", stalk="#0000ff")
	cytobands$color <- bandColors[cytobands$giemsa]
	arcSteps <- 30 # defines roundness of arc
	curlyBraceHeight <- 0.03
	ideogramHeight <- 0.04
	ideogramWidth <- 0.4
	# Extract bands of given contig
	bands <- cytobands[cytobands$contig==contig,]
	if (nrow(bands) == 0) {
		warning(paste("Ideogram of contig", contig, "cannot be drawn, because no Giemsa staining information is available."))
		return(NULL)
	}
	# Scale width of ideogram to fit inside given region
	bands$left <- bands$start / max(cytobands$end) * ideogramWidth
	bands$right <- bands$end / max(cytobands$end) * ideogramWidth
	# Seft/right-align cytobands
	offset <- ifelse(adjust=="left", left, right - max(bands$right))
	bands$left <- bands$left + offset
	bands$right <- bands$right + offset
	# Draw curly braces
	tip <- min(bands$left) + (max(bands$right)-min(bands$left)) / (max(bands$end)-min(bands$start)) * breakpoint
	drawCurlyBrace(left, right, y-0.05+curlyBraceHeight, y-0.05, tip)
	# Draw title of chromosome
	text((max(bands$right)+min(bands$left))/2, y+0.07, paste("Chromosome", contig), font=2, cex=fontSize, adj=c(0.5,0))
	# Draw name of band
	bandName <- bands[which(between(breakpoint, bands$start, bands$end)), "name"]
	text(tip, y+0.03, bandName, cex=fontSize, adj=c(0.5,0))
	# Draw start of chromosome
	leftArcX <- bands[1,"left"] + (1+cos(seq(pi/2,1.5*pi,len=arcSteps))) * (bands[1,"right"]-bands[1,"left"])
	leftArcY <- y + sin(seq(pi/2,1.5*pi,len=arcSteps)) * (ideogramHeight/2)
	polygon(leftArcX, leftArcY, col=bands[1,"color"])
	# Draw bands
	centromereStart <- NULL
	centromereEnd <- NULL
	for (band in 2:(nrow(bands)-1)) {
		if (bands[band,"giemsa"] != "acen") {
			rect(bands[band,"left"], y-ideogramHeight/2, bands[band,"right"], y+ideogramHeight/2, col=bands[band,"color"])
		} else { # draw centromere
			if (is.null(centromereStart)) {
				polygon(c(bands[band,"left"], bands[band,"right"], bands[band,"left"]), c(y-ideogramHeight/2, y, y+ideogramHeight/2), col=bands[band,"color"])
				centromereStart <- bands[band,"left"]
			} else {
				polygon(c(bands[band,"right"], bands[band,"left"], bands[band,"right"]), c(y-ideogramHeight/2, y, y+ideogramHeight/2), col=bands[band,"color"])
				centromereEnd <- bands[band,"right"]
			}
		}
	}
	# Draw end of chromosome
	band <- nrow(bands)
	rightArcX <- bands[band,"right"] - (1+cos(seq(1.5*pi,pi/2,len=arcSteps))) * (bands[band,"right"]-bands[band,"left"])
	rightArcY <- y + sin(seq(pi/2,1.5*pi,len=arcSteps)) * ideogramHeight/2
	polygon(rightArcX, rightArcY, col=bands[band,"color"])
	# If there is no centromere, make an artificial one with length zero
	if (is.null(centromereStart) || is.null(centromereEnd)) {
		centromereStart <- bands[1,"right"]
		centromereEnd <- bands[1,"right"]
	}
	# Draw gradients for 3D effect
	drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(0,0,0,0.8), 1:round(arcSteps*0.4)) # black from top on p-arm
	drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(1,1,1,0.7), round(arcSteps*0.4):round(arcSteps*0.1)) # white to top on p-arm
	drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(1,1,1,0.7), round(arcSteps*0.4):round(arcSteps*0.6)) # white to bottom on p-arm
	drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(0,0,0,0.9), arcSteps:round(arcSteps*0.5)) # black from bottom on p-arm
	drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(0,0,0,0.8), 1:round(arcSteps*0.4)) # black from top on q-arm
	drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(1,1,1,0.7), round(arcSteps*0.4):round(arcSteps*0.1)) # white to top on q-arm
	drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(1,1,1,0.7), round(arcSteps*0.4):round(arcSteps*0.6)) # white to bottom on q-arm
	drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(0,0,0,0.9), arcSteps:round(arcSteps*0.5)) # black from bottom on q-arm
}

drawCoverage <- function(left, right, y, coverage, start, end, color) {
	maxResolution <- 5000 # max number of data points to draw coverage
	# Draw coverage as bars
	if (!is.null(coverage)) {
		coverageData <- as.numeric(coverage[IRanges(sapply(start, max, min(start(coverage))), sapply(end, min, max(end(coverage))))])
		# Downsample to maxResolution, if there are too many data points
		coverageData <- aggregate(coverageData, by=list(round(1:length(coverageData) * (right-left) * maxResolution/length(coverageData))), mean)$x
		polygon(c(left, seq(left, right, length.out=length(coverageData)), right), c(y, y+coverageData*0.1, y), col=color, border=NA)
	}
}

drawStrand <- function(left, right, y, color, strand) {
	if (strand %in% c("+", "-")) {
		# Draw strand
		lines(c(left+0.001, right-0.001), c(y, y), col=color, lwd=2)
		lines(c(left+0.001, right-0.001), c(y, y), col=rgb(1,1,1,0.1), lwd=1)
		# Indicate orientation
		if (right - left > 0.01)
			for (i in seq(left+0.005, right-0.005, by=sign(right-left-2*0.005)*0.01)) {
				arrows(i, y, i+0.001*ifelse(strand=="+", 1, -1), y, col=color, length=0.05, lwd=2, angle=60)
				arrows(i, y, i+0.001*ifelse(strand=="+", 1, -1), y, col=rgb(1,1,1,0.1), length=0.05, lwd=1, angle=60)
			}
	}
}

drawExon <- function(left, right, y, color, title, type) {
	gradientSteps <- 10 # defines smoothness of gradient
	exonHeight <- 0.03
	if (type == "CDS") {
		# Draw coding regions as thicker bars
		rect(left, y+exonHeight, right, y+exonHeight/2-0.001, col=color, border=NA)
		rect(left, y-exonHeight, right, y-exonHeight/2+0.001, col=color, border=NA)
		# Draw border
		lines(c(left, left, right, right), c(y+exonHeight/2, y+exonHeight, y+exonHeight, y+exonHeight/2), col=getDarkColor(color), lend=2)
		lines(c(left, left, right, right), c(y-exonHeight/2, y-exonHeight, y-exonHeight, y-exonHeight/2), col=getDarkColor(color), lend=2)
		# Draw gradients for 3D effect
		drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y+0.03, y+0.015, len=gradientSteps), rgb(0,0,0,0.2))
		drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y-0.03, y-0.015, len=gradientSteps), rgb(0,0,0,0.3))
	} else if (type == "exon") {
		rect(left, y+exonHeight/2, right, y-exonHeight/2, col=color, border=getDarkColor(color))
		# Draw gradients for 3D effect
		drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y, y+exonHeight/2, len=gradientSteps), rgb(1,1,1,0.6))
		drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y, y-exonHeight/2, len=gradientSteps), rgb(1,1,1,0.6))
		# Add exon label
		text((left+right)/2, y, title, cex=0.9*fontSize)
	}
}

drawCircos <- function(SV, SVs, cytobands, minConfidenceForCircosPlot, circosColors) {
	# Check if Giemsa staining information is available
	for (contig in unlist(SVs[SV,c("contig1", "contig2")])) {
		if (!any(cytobands$contig==contig)) {
			warning(paste0("Circos plot cannot be drawn, because no Giemsa staining information is available for contig ", contig, "."))
			# Draw empty plots as placeholder
			plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
			plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
			return(NULL)
		}
	}
	# Initialize with empty circos plot
	circos.clear()
	circos.initializeWithIdeogram(cytoband=cytobands, plotType=NULL)
	# Use gene names as labels or <contig>:<position> for intergenic breakpoints
	geneLabels <- data.frame(
		contig=c(SVs[SV,"contig1"], SVs[SV,"contig2"]),
		start=c(SVs[SV,"Breakpoint1"], SVs[SV,"Breakpoint2"])
	)
	geneLabels$end <- geneLabels$start + 1
	geneLabels$gene <- c(SVs[SV,"gene1"], SVs[SV,"gene2"])
	geneLabels$gene <- ifelse(c(SVs[SV,"site1"], SVs[SV,"site2"]) == "intergenic", paste0(c(SVs[SV,"display_contig1"], SVs[SV,"display_contig2"]), ":", geneLabels$start), geneLabels$gene)
	# Draw gene labels
	circos.genomicLabels(geneLabels, labels.column=4, side="outside", cex=fontSize, font=3, labels_height=0.27)
	# Draw chromosome labels in connector plot
	for (contig in unique(cytobands$contig)) {
		set.current.cell(track.index=2, sector.index=contig) # draw in gene label connector track (track.index=2)
		circos.text(CELL_META$xcenter, CELL_META$ycenter, contig, cex=0.85)
	}
	# Draw ideograms
	circos.genomicIdeogram(cytoband=cytobands)
	# Draw arcs
	confidenceRank <- c(Low=0, Medium=1, High=2)
	for (i in c(setdiff(1:nrow(SVs), SV), SV)) { # draw SV of interest last, such that its arc is on top
		f <- SVs[i,]
		if (any(cytobands$contig==f$contig1) && any(cytobands$contig==f$contig2)) # ignore viral contigs, because we have no cytoband information for them
			if (minConfidenceForCircosPlot != "none" && confidenceRank[f$confidence] >= confidenceRank[minConfidenceForCircosPlot] || i==SV)
				circos.link(
					f$contig1, f$Breakpoint1,
					f$contig2, f$Breakpoint2,
					lwd=2, col=ifelse(i==SV, circosColors[f$type], getBrightColor(circosColors[f$type]))
				)
	}
	# Draw legend
	plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", ylab="", xlab="")
	legend(x="top", legend=names(circosColors), col=sapply(circosColors, getBrightColor), lwd=3, ncol=2, box.lty=0)
}

drawProteinDomains <- function(SV, exons1, exons2, proteinDomains, color1, color2, mergeDomainsOverlappingBy, optimizeDomainColors) {

	exonHeight <- 0.2
	exonsY <- 0.5
	geneNamesY <- exonsY - exonHeight/2 - 0.05

	# Find coding exons
	codingExons1 <- exons1[exons1$type == "CDS" & SV$site1 != "intergenic",]
	codingExons2 <- exons2[exons2$type == "CDS" & SV$site2 != "intergenic",]

	# Cut off coding regions beyond breakpoint
	if (SV$direction1 == "upstream") {
		codingExons1 <- codingExons1[codingExons1$end >= SV$Breakpoint1,]
		codingExons1$start <- ifelse(codingExons1$start < SV$Breakpoint1, SV$Breakpoint1, codingExons1$start)
	} else {
		codingExons1 <- codingExons1[codingExons1$start <= SV$Breakpoint1,]
		codingExons1$end <- ifelse(codingExons1$end > SV$Breakpoint1, SV$Breakpoint1, codingExons1$end)
	}
	if (SV$direction2 == "upstream") {
		codingExons2 <- codingExons2[codingExons2$end >= SV$Breakpoint2,]
		codingExons2$start <- ifelse(codingExons2$start < SV$Breakpoint2, SV$Breakpoint2, codingExons2$start)
	} else {
		codingExons2 <- codingExons2[codingExons2$start <= SV$Breakpoint2,]
		codingExons2$end <- ifelse(codingExons2$end > SV$Breakpoint2, SV$Breakpoint2, codingExons2$end)
	}

	# Find overlapping domains
	exonsGRanges1 <- GRanges(codingExons1$contig, IRanges(codingExons1$start, codingExons1$end), strand=codingExons1$strand)
	exonsGRanges2 <- GRanges(codingExons2$contig, IRanges(codingExons2$start, codingExons2$end), strand=codingExons2$strand)
	domainsGRanges <- GRanges(proteinDomains$contig, IRanges(proteinDomains$start, proteinDomains$end), strand=proteinDomains$strand)
	domainsGRanges$proteinDomainName <- proteinDomains$proteinDomainName
	domainsGRanges$proteinDomainID <- proteinDomains$proteinDomainID
	domainsGRanges$color <- proteinDomains$color
	domainsGRanges <- domainsGRanges[suppressWarnings(unique(queryHits(findOverlaps(domainsGRanges, union(exonsGRanges1, exonsGRanges2)))))]

	# Group overlapping domains by domain ID
	domainsGRangesList <- GRangesList(lapply(unique(domainsGRanges$proteinDomainID), function(x) { domainsGRanges[domainsGRanges$proteinDomainID == x] }))

	# Trim protein domains to exon boundaries
	trimDomains <- function(domainsGRangesList, exonsGRanges) {
		do.call(
			"rbind",
			lapply(
				domainsGRangesList,
				function(x) {
					intersected <- as.data.frame(reduce(suppressWarnings(intersect(x, exonsGRanges))))
					if (nrow(intersected) > 0) {
						intersected$proteinDomainName <- head(x$proteinDomainName, 1)
						intersected$proteinDomainID <- head(x$proteinDomainID, 1)
						intersected$color <- head(x$color, 1)
					} else {
						intersected$proteinDomainName <- character()
						intersected$proteinDomainID <- character()
						intersected$color <- character()
					}
					return(intersected)
				}
			)
		)
	}
	retainedDomains1 <- trimDomains(domainsGRangesList, exonsGRanges1)
	retainedDomains2 <- trimDomains(domainsGRangesList, exonsGRanges2)

	# Calculate length of coding exons
	codingExons1$length <- codingExons1$end - codingExons1$start + 1
	codingExons2$length <- codingExons2$end - codingExons2$start + 1

	# Abort, if there are no coding regions
	if (sum(exons1$type == "CDS") + sum(exons2$type == "CDS") == 0) {
		text(0.5, 0.5, "Genes are not protein-coding.")
		return(NULL)
	}
	codingLength1 <- sum(codingExons1$length)
	codingLength2 <- sum(codingExons2$length)
	if (codingLength1 + codingLength2 == 0) {
		text(0.5, 0.5, "No coding regions retained in SV transcript.")
		return(NULL)
	}
	antisenseTranscription1 <- sub("/.*", "", SV$strand1) != sub(".*/", "", SV$strand1)
	antisenseTranscription2 <- sub("/.*", "", SV$strand2) != sub(".*/", "", SV$strand2)
	if ((codingLength1 == 0 || antisenseTranscription1) && (codingLength2 == 0 || antisenseTranscription2)) {
		text(0.5, 0.5, "No coding regions due to antisense transcription.")
		return(NULL)
	}

	# Remove introns from protein domains
	removeIntronsFromProteinDomains <- function(codingExons, retainedDomains) {
		if (nrow(codingExons) == 0) return(NULL)
		cumulativeIntronLength <- 0
		previousExonEnd <- 0
		for (exon in 1:nrow(codingExons)) {
			if (codingExons[exon,"start"] > previousExonEnd)
				cumulativeIntronLength <- cumulativeIntronLength + codingExons[exon,"start"] - previousExonEnd
			domainsInExon <- which(between(retainedDomains$start, codingExons[exon,"start"], codingExons[exon,"end"]))
			retainedDomains[domainsInExon,"start"] <- retainedDomains[domainsInExon,"start"] - cumulativeIntronLength
			domainsInExon <- which(between(retainedDomains$end, codingExons[exon,"start"], codingExons[exon,"end"]))
			retainedDomains[domainsInExon,"end"] <- retainedDomains[domainsInExon,"end"] - cumulativeIntronLength
			previousExonEnd <- codingExons[exon,"end"]
		}
		# Merge adjacent domains
		retainedDomains <- do.call(
			"rbind",
			lapply(
				unique(retainedDomains$proteinDomainID),
				function(x) {
					domain <- retainedDomains[retainedDomains$proteinDomainID == x,]
					merged <- reduce(GRanges(domain$seqnames, IRanges(domain$start, domain$end), strand=domain$strand))
					merged$proteinDomainName <- head(domain$proteinDomainName, 1)
					merged$proteinDomainID <- head(domain$proteinDomainID, 1)
					merged$color <- head(domain$color, 1)
					return(as.data.frame(merged))
				}
			)
		)
		return(retainedDomains)
	}
	retainedDomains1 <- removeIntronsFromProteinDomains(codingExons1, retainedDomains1)
	retainedDomains2 <- removeIntronsFromProteinDomains(codingExons2, retainedDomains2)

	# Abort, if no domains are retained
	if (is.null(retainedDomains1) && is.null(retainedDomains2)) {
		text(0.5, 0.5, "No protein domains retained in SV.")
		return(NULL)
	}

	# Merge domains with similar coordinates
	mergeSimilarDomains <- function(domains, mergeDomainsOverlappingBy) {
		if (is.null(domains)) return(domains)
		merged <- domains[F,] # create empty data frame
		domains <- domains[order(domains$end - domains$start, decreasing=T),] # start with bigger domains => bigger domains are retained
		for (domain in rownames(domains)) {
			if (!any((abs(merged$start - domains[domain,"start"]) + abs(merged$end - domains[domain,"end"])) / (domains[domain,"end"] - domains[domain,"start"] + 1) <= 1-mergeDomainsOverlappingBy))
				merged <- rbind(merged, domains[domain,])
		}
		return(merged)
	}
	retainedDomains1 <- mergeSimilarDomains(retainedDomains1, mergeDomainsOverlappingBy)
	retainedDomains2 <- mergeSimilarDomains(retainedDomains2, mergeDomainsOverlappingBy)

	# If desired, reassign colors to protein domains to maximize contrast
	if (optimizeDomainColors) {
		uniqueDomains <- unique(c(retainedDomains1$proteinDomainID, retainedDomains2$proteinDomainID))
		# Make rainbow of pretty pastell colors
		colors <- rainbow(length(uniqueDomains))
		colors <- apply(col2rgb(colors), 2, function(x) { 0.3 + y/255 * 0.7 }) # make pastell colors
		colors <- apply(colors, 2, function(x) {rgb(x["red"], x["green"], x["blue"])}) # convert back to rgb
		# Reassign colors
		names(colors) <- uniqueDomains
		retainedDomains1$color <- colors[retainedDomains1$proteinDomainID]
		retainedDomains2$color <- colors[retainedDomains2$proteinDomainID]
	}

	# Reverse exons and protein domains, if on the reverse strand
	if (any(codingExons1$strand == "-")) {
		codingExons1$length <- rev(codingExons1$length)
		temp <- retainedDomains1$end
		retainedDomains1$end <- codingLength1 - retainedDomains1$start
		retainedDomains1$start <- codingLength1 - temp
	}
	if (any(codingExons2$strand == "-")) {
		codingExons2$length <- rev(codingExons2$length)
		temp <- retainedDomains2$end
		retainedDomains2$end <- codingLength2 - retainedDomains2$start
		retainedDomains2$start <- codingLength2 - temp
	}

	# Normalize length to 1
	codingExons1$length <- codingExons1$length / (codingLength1 + codingLength2)
	codingExons2$length <- codingExons2$length / (codingLength1 + codingLength2)
	retainedDomains1$start <- retainedDomains1$start / (codingLength1 + codingLength2)
	retainedDomains1$end <- retainedDomains1$end / (codingLength1 + codingLength2)
	retainedDomains2$start <- retainedDomains2$start / (codingLength1 + codingLength2)
	retainedDomains2$end <- retainedDomains2$end / (codingLength1 + codingLength2)

	# Draw coding regions
	rect(0, exonsY-exonHeight/2, sum(codingExons1$length), exonsY+exonHeight/2, col=color1, border=NA)
	rect(sum(codingExons1$length), exonsY-exonHeight/2, sum(codingExons1$length) + sum(codingExons2$length), exonsY+exonHeight/2, col=color2, border=NA)

	# Indicate exon boundaries as dotted lines
	exonBoundaries <- cumsum(c(codingExons1$length, codingExons2$length))
	if (length(exonBoundaries) > 1) {
		exonBoundaries <- exonBoundaries[1:(length(exonBoundaries)-1)]
		for (exonBoundary in exonBoundaries)
			lines(c(exonBoundary, exonBoundary), c(exonsY-exonHeight, exonsY+exonHeight), col="white", lty=3)
	}

	# Find overlapping domains
	# Nest if one is contained in another
	# Stack if they overlap partially
	nestDomains <- function(domains) {
		if (length(unlist(domains)) == 0) return(domains)
		domains <- domains[order(domains$end - domains$start, decreasing=T),]
		rownames(domains) <- 1:nrow(domains)
		# Find nested domains and make tree structure
		domains$parent <- 0
		for (domain in rownames(domains))
			domains[domains$start >= domains[domain,"start"] & domains$end <= domains[domain,"end"] & rownames(domains) != domain,"parent"] <- domain
		# Find partially overlapping domains
		maxOverlappingDomains <- max(1, as.integer(coverage(IRanges(domains$start*10e6, domains$end*10e6))))
		padding <- 1 / maxOverlappingDomains * 0.4
		domains$y <- 0
		domains$height <- 0
		adjustPositionAndHeight <- function(parentDomain, y, height, padding, e) {
			for (domain in which(e$domains$parent == parentDomain)) {
				overlappingDomains <- which((between(e$domains$start, e$domains[domain,"start"], e$domains[domain,"end"]) |
				                             between(e$domains$end  , e$domains[domain,"start"], e$domains[domain,"end"])) &
				                             e$domains$parent == parentDomain)
				e$domains[domain,"height"] <- height/length(overlappingDomains) - padding * (length(overlappingDomains)-1) / length(overlappingDomains)
				e$domains[domain,"y"] <- y + (which(domain==overlappingDomains)-1) * (e$domains[domain,"height"] + padding)
				adjustPositionAndHeight(domain, e$domains[domain,"y"]+padding, e$domains[domain,"height"]-2*padding, padding, e)
			}
		}
		adjustPositionAndHeight(0, 0, 1, padding, environment())
		domains <- domains[order(domains$height, decreasing=T),] # draw nested domains last
		return(domains)
	}
	retainedDomains1 <- nestDomains(retainedDomains1)
	retainedDomains2 <- nestDomains(retainedDomains2)
	retainedDomains1$y <- exonsY - exonHeight/2 + 0.025 + (exonHeight-2*0.025) * retainedDomains1$y
	retainedDomains2$y <- exonsY - exonHeight/2 + 0.025 + (exonHeight-2*0.025) * retainedDomains2$y
	retainedDomains1$height <- retainedDomains1$height * (exonHeight-2*0.025)
	retainedDomains2$height <- retainedDomains2$height * (exonHeight-2*0.025)

	# Draw domains
	drawProteinDomainRect <- function(left, bottom, right, top, color) {
		rect(left, bottom, right, top, col=color, border=getDarkColor(color))
		# Draw gradients for 3D effect
		gradientSteps <- 20
		drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(top, bottom, len=gradientSteps), rgb(1,1,1,0.7))
		drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(bottom, bottom+(top-bottom)*0.4, len=gradientSteps), rgb(0,0,0,0.1))
	}
	if (length(unlist(retainedDomains1)) > 0)
		for (domain in 1:nrow(retainedDomains1))
			drawProteinDomainRect(retainedDomains1[domain,"start"], retainedDomains1[domain,"y"], retainedDomains1[domain,"end"], retainedDomains1[domain,"y"]+retainedDomains1[domain,"height"], retainedDomains1[domain,"color"])
	if (length(unlist(retainedDomains2)) > 0)
		for (domain in 1:nrow(retainedDomains2))
			drawProteinDomainRect(sum(codingExons1$length)+retainedDomains2[domain,"start"], retainedDomains2[domain,"y"], sum(codingExons1$length)+retainedDomains2[domain,"end"], retainedDomains2[domain,"y"]+retainedDomains2[domain,"height"], retainedDomains2[domain,"color"])

	# Draw gene names, if there are coding exons
	if (codingLength1 > 0)
		text(sum(codingExons1$length)/2, geneNamesY, SV$gene1, font=4, cex=fontSize)
	if (codingLength2 > 0)
		text(sum(codingExons1$length)+sum(codingExons2$length)/2, geneNamesY, SV$gene2, font=4, cex=fontSize)

	# Calculate how many non-adjacent unique domains there are
	# We need this info to know where to place labels vertically
	countUniqueDomains <- function(domains) {
		uniqueDomains <- 0
		if (length(unlist(domains)) > 0) {
			uniqueDomains <- 1
			if (nrow(domains) > 1) {
				previousDomain <- domains[1,"proteinDomainID"]
				for (domain in 2:nrow(domains)) {
					if (previousDomain != domains[domain,"proteinDomainID"])
						uniqueDomains <- uniqueDomains + 1
					previousDomain <- domains[domain,"proteinDomainID"]
				}
			}
		}
		return(uniqueDomains)
	}
	if (length(unlist(retainedDomains1)) > 0)
		retainedDomains1 <- retainedDomains1[order(retainedDomains1$start),]
	uniqueDomains1 <- countUniqueDomains(retainedDomains1)
	if (length(unlist(retainedDomains2)) > 0)
		retainedDomains2 <- retainedDomains2[order(retainedDomains2$end, decreasing=T),]
	uniqueDomains2 <- countUniqueDomains(retainedDomains2)

	# draw title of plot
	titleY <- exonsY + exonHeight/2 + (uniqueDomains1 + 2) * 0.05
	text(0.5, titleY+0.01, "RETAINED PROTEIN DOMAINS", adj=c(0.5, 0), font=2, cex=fontSize)
	text(0.5, titleY, ifelse(SV$reading_frame %in% c("in-frame", "out-of-frame"), paste(SV$reading_frame, "SV"), ifelse(SV$reading_frame == "stop-codon", "stop codon before SV junction", "Reading frame unclear")), adj=c(0.5, 1), cex=fontSize)

	# draw domain labels for gene1
	if (length(unlist(retainedDomains1)) > 0) {
		previousConnectorX <- -1
		previousLabelX <- -1
		labelY <- exonsY + exonHeight/2 + uniqueDomains1 * 0.05
		for (domain in 1:nrow(retainedDomains1)) {
			# if possible avoid overlapping lines of labels
			connectorX <- min(retainedDomains1[domain,"start"] + 0.01, (retainedDomains1[domain,"start"] + retainedDomains1[domain,"end"])/2)
			if (connectorX - previousConnectorX < 0.01 && retainedDomains1[domain,"end"] > previousConnectorX + 0.01)
				connectorX <- previousConnectorX + 0.01
			labelX <- max(connectorX, previousLabelX) + 0.02
			# use a signle label for adjacent domains of same type
			adjacentDomainsOfSameType <- domain + 1 <= nrow(retainedDomains1) && retainedDomains1[domain+1,"proteinDomainID"] == retainedDomains1[domain,"proteinDomainID"]
			if (adjacentDomainsOfSameType) {
				labelX <- retainedDomains1[domain+1,"start"] + 0.015
			} else {
				text(labelX, labelY, retainedDomains1[domain,"proteinDomainName"], adj=c(0,0.5), col=getDarkColor(retainedDomains1[domain,"color"]), cex=fontSize)
			}
			lines(c(labelX-0.005, connectorX, connectorX), c(labelY, labelY, retainedDomains1[domain,"y"]+retainedDomains1[domain,"height"]), col=getDarkColor(retainedDomains1[domain,"color"]))
			if (!adjacentDomainsOfSameType)
				labelY <- labelY - 0.05
			previousConnectorX <- connectorX
			previousLabelX <- labelX
		}
	}

	# draw domain labels for gene2
	if (length(unlist(retainedDomains2)) > 0) {
		previousConnectorX <- 100
		previousLabelX <- 100
		labelY <- exonsY - exonHeight/2 - (uniqueDomains2+1) * 0.05
		for (domain in 1:nrow(retainedDomains2)) {
			# if possible avoid overlapping connector lines of labels
			connectorX <- sum(codingExons1$length) + max(retainedDomains2[domain,"end"] - 0.01, (retainedDomains2[domain,"start"] + retainedDomains2[domain,"end"])/2)
			if (previousConnectorX - connectorX < 0.01 && sum(codingExons1$length) + retainedDomains2[domain,"start"] < previousConnectorX - 0.01)
				connectorX <- previousConnectorX - 0.01
			labelX <- min(connectorX, previousLabelX) - 0.02
			# use a signle label for adjacent domains of same type
			adjacentDomainsOfSameType <- domain + 1 <= nrow(retainedDomains2) && retainedDomains2[domain+1,"proteinDomainID"] == retainedDomains2[domain,"proteinDomainID"]
			if (adjacentDomainsOfSameType) {
				labelX <- sum(codingExons1$length) + retainedDomains2[domain+1,"end"] - 0.015
			} else {
				text(labelX, labelY, retainedDomains2[domain,"proteinDomainName"], adj=c(1,0.5), col=getDarkColor(retainedDomains2[domain,"color"]), cex=fontSize)
			}
			lines(c(labelX+0.005, connectorX, connectorX), c(labelY, labelY, retainedDomains2[domain,"y"]), col=getDarkColor(retainedDomains2[domain,"color"]))
			if (!adjacentDomainsOfSameType)
				labelY <- labelY + 0.05
			previousConnectorX <- connectorX
			previousLabelX <- labelX
		}
	}

}

findExons <- function(exons, contig, geneID, direction, breakpoint, coverage, transcriptId, transcriptSelection) {
	# use the provided transcript if desired
	if (transcriptSelection == "provided" && transcriptId != "." && transcriptId != "") {
		candidateExons <- exons[exons$transcript == transcriptId,]
		if (nrow(candidateExons) == 0) {
			warning(paste0("Unknown transcript given in SVs file (", transcriptId, "), selecting a different one"))
		} else {
			return(candidateExons)
		}
	}

	if (transcriptSelection == "canonical") {
		candidateExons <- exons[exons$geneID == geneID & exons$contig == contig,]
	} else {
		# look for exon with breakpoint as splice site
		transcripts <- exons[exons$geneID == geneID & exons$contig == contig & exons$type == "exon" & (direction == "downstream" & abs(exons$end - breakpoint) <= 2 | direction == "upstream" & abs(exons$start - breakpoint) <= 2),"transcript"]
		candidateExons <- exons[exons$transcript %in% transcripts,]
		# if none was found, use all exons of the gene closest to the breakpoint
		if (nrow(candidateExons) == 0)
			candidateExons <- exons[exons$geneID == geneID & exons$contig == contig,]
		# if we have coverage information, use the transcript with the highest coverage if there are multiple hits
		if (!is.null(coverage)) {
			highestCoverage <- -1
			transcriptWithHighestCoverage <- NULL
			lengthOfTranscriptWithHighestCoverage <- 0
			for (transcript in unique(candidateExons$transcript)) {
				exonsOfTranscript <- candidateExons[candidateExons$transcript==transcript,]
				exonsOfTranscript$start <- sapply(exonsOfTranscript$start, max, min(start(coverage)))
				exonsOfTranscript$end <- sapply(exonsOfTranscript$end, min, max(end(coverage)))
				lengthOfTranscript <- sum(exonsOfTranscript$end - exonsOfTranscript$start + 1)
				coverageSum <- sum(as.numeric(coverage[IRanges(exonsOfTranscript$start, exonsOfTranscript$end)]))
				# we prefer shorter transcripts over longer ones, because otherwise there is a bias towards transcripts with long UTRs
				# => a longer transcript must have substantially higher coverage to replace a shorter one
				substantialDifference <- (1 - min(lengthOfTranscript, lengthOfTranscriptWithHighestCoverage) / max(lengthOfTranscript, lengthOfTranscriptWithHighestCoverage)) / 10
				if (lengthOfTranscript > lengthOfTranscriptWithHighestCoverage && coverageSum * (1-substantialDifference) > highestCoverage ||
				    lengthOfTranscript < lengthOfTranscriptWithHighestCoverage && coverageSum > highestCoverage * (1-substantialDifference)) {
					highestCoverage <- coverageSum
					transcriptWithHighestCoverage <- transcript
					lengthOfTranscriptWithHighestCoverage <- lengthOfTranscript
				}
			}
			candidateExons <- candidateExons[candidateExons$transcript==transcriptWithHighestCoverage,]
		}
		# if the gene has multiple transcripts, search for transcripts which encompass the breakpoint
		if (length(unique(candidateExons$transcript)) > 1) {
			transcriptStart <- aggregate(candidateExons$start, by=list(candidateExons$transcript), min)
			rownames(transcriptStart) <- transcriptStart[,1]
			transcriptEnd <- aggregate(candidateExons$end, by=list(candidateExons$transcript), max)
			rownames(transcriptEnd) <- transcriptEnd[,1]
			candidateExons <- candidateExons[between(breakpoint, transcriptStart[candidateExons$transcript,2], transcriptEnd[candidateExons$transcript,2]),]
		}
	}

	# find the consensus transcript, if there are multiple hits
	if (length(unique(candidateExons$transcript)) > 1) {
		consensusTranscript <-
			ifelse(grepl("appris_principal_1", candidateExons$attributes), 12,
			ifelse(grepl("appris_principal_2", candidateExons$attributes), 11,
			ifelse(grepl("appris_principal_3", candidateExons$attributes), 10,
			ifelse(grepl("appris_principal_4", candidateExons$attributes), 9,
			ifelse(grepl("appris_principal_5", candidateExons$attributes), 8,
			ifelse(grepl("appris_principal", candidateExons$attributes), 7,
			ifelse(grepl("appris_candidate_longest", candidateExons$attributes), 6,
			ifelse(grepl("appris_candidate", candidateExons$attributes), 5,
			ifelse(grepl("appris_alternative_1", candidateExons$attributes), 4,
			ifelse(grepl("appris_alternative_2", candidateExons$attributes), 3,
			ifelse(grepl("appris_alternative", candidateExons$attributes), 2,
			ifelse(grepl("CCDS", candidateExons$attributes), 1,
			0
		))))))))))))
		candidateExons <- candidateExons[consensusTranscript == max(consensusTranscript),]
	}
	# use the transcript with the longest coding sequence, if there are still multiple hits
	if (length(unique(candidateExons$transcript)) > 1) {
		codingSequenceLength <- ifelse(candidateExons$type == "CDS", candidateExons$end - candidateExons$start, 0)
		totalCodingSequenceLength <- aggregate(codingSequenceLength, by=list(candidateExons$transcript), sum)
		rownames(totalCodingSequenceLength) <- totalCodingSequenceLength[,1]
		candidateExons <- candidateExons[totalCodingSequenceLength[candidateExons$transcript,2] == max(totalCodingSequenceLength[,2]),]
	}
	# use the transcript with the longest overall sequence, if there are still multiple hits
	if (length(unique(candidateExons$transcript)) > 1) {
		exonLength <- candidateExons$end - candidateExons$start
		totalExonLength <- aggregate(exonLength, by=list(candidateExons$transcript), sum)
		rownames(totalExonLength) <- totalExonLength[,1]
		candidateExons <- candidateExons[totalExonLength[candidateExons$transcript,2] == max(totalExonLength[,2]),]
	}
	# if there are still multiple hits, select the first one
	candidateExons <- unique(candidateExons[candidateExons$transcript == head(unique(candidateExons$transcript), 1),])
	return(candidateExons)
}

findClosestGene <- function(exons, contig, breakpoint, extraConditions) {

	# find exons near breakpoint (extraConditions must define what is considered "near")
	closestExons <- exons[exons$contig == contig & extraConditions,] # find closest exon
	closestExons <- exons[exons$contig == contig & exons$geneID %in% closestExons$geneID,] # select all exons of closest gene

	# when more than one gene found with the given name, use the closest one
	if (length(unique(closestExons$geneID)) > 1) { # more than one gene found with the given name => use the closest one
		distanceToBreakpoint <- aggregate(1:nrow(closestExons), by=list(closestExons$geneID), function(x) { min(abs(closestExons[x,"start"]-breakpoint), abs(closestExons[x,"end"]-breakpoint)) })
		closestGene <- head(distanceToBreakpoint[distanceToBreakpoint[,2] == min(distanceToBreakpoint[,2]),1], 1)
		closestExons <- closestExons[closestExons$geneID == closestGene,]
	}

	# when no gene was found, return default values
	if (nrow(closestExons) == 0) {
		return(IRanges(max(1, breakpoint-1000), breakpoint+1000))
	} else {
		return(IRanges(min(closestExons$start), max(closestExons$end)))
	}
}

# main loop starts here
for (SV in 1:nrow(SVs)) {

	message(paste0("Drawing SV #", SV, ": ", SVs[SV,"gene1"], ":", SVs[SV,"gene2"]))

	# if showIntergenicVicinity is a number, take it as is
	# if it is a keyword (closestGene/closestProteinCodingGene), determine the range dynamically
	showVicinity <- rep(0, 4)
	if (SVs[SV,"site1"] == "intergenic") {
		showVicinity[1] <- ifelse(
			is.numeric(showIntergenicVicinity[[1]]),
			showIntergenicVicinity[[1]],
			SVs[SV,"Breakpoint1"] - start(findClosestGene(exons, SVs[SV,"contig1"], SVs[SV,"Breakpoint1"], exons$end < SVs[SV,"Breakpoint1"] & exons$type == showIntergenicVicinity[[1]]))
		)
		showVicinity[2] <- ifelse(
			is.numeric(showIntergenicVicinity[[2]]),
			showIntergenicVicinity[[2]],
			end(findClosestGene(exons, SVs[SV,"contig1"], SVs[SV,"Breakpoint1"], exons$start > SVs[SV,"Breakpoint1"] & exons$type == showIntergenicVicinity[[2]])) - SVs[SV,"Breakpoint1"]
		)
	}
	if (SVs[SV,"site2"] == "intergenic") {
		showVicinity[3] <- ifelse(
			is.numeric(showIntergenicVicinity[[3]]),
			showIntergenicVicinity[[3]],
			SVs[SV,"Breakpoint2"] - start(findClosestGene(exons, SVs[SV,"contig2"], SVs[SV,"Breakpoint2"], exons$end < SVs[SV,"Breakpoint2"] & exons$type == showIntergenicVicinity[[3]]))
		)
		showVicinity[4] <- ifelse(
			is.numeric(showIntergenicVicinity[[4]]),
			showIntergenicVicinity[[4]],
			end(findClosestGene(exons, SVs[SV,"contig2"], SVs[SV,"Breakpoint2"], exons$start > SVs[SV,"Breakpoint2"] & exons$type == showIntergenicVicinity[[4]])) - SVs[SV,"Breakpoint2"]
		)
	}

	# compute coverage from alignments file
	coverage1 <- NULL
	coverage2 <- NULL
	if (alignmentsFile != "") {
		# determine range in which we need to compute the coverage
		determineCoverageRegion <- function(exons, geneID, contig, breakpoint, showVicinityLeft, showVicinityRight) {
			closestGene <- findClosestGene(exons, contig, breakpoint, exons$geneID == geneID)
			return(IRanges(min(start(closestGene), breakpoint-showVicinityLeft), max(end(closestGene), breakpoint+showVicinityRight)))
		}
		coverageRegion1 <- determineCoverageRegion(exons, SVs[SV,"gene_id1"], SVs[SV,"contig1"], SVs[SV,"Breakpoint1"], showVicinity[1], showVicinity[2])
		coverageRegion2 <- determineCoverageRegion(exons, SVs[SV,"gene_id2"], SVs[SV,"contig2"], SVs[SV,"Breakpoint2"], showVicinity[3], showVicinity[4])
		# function which reads alignments from BAM file with & without "chr" prefix
		readCoverage <- function(alignmentsFile, contig, coverageRegion) {
			coverageData <- tryCatch(
				{
					alignments <- readGAlignments(alignmentsFile, param=ScanBamParam(which=GRanges(contig, coverageRegion)))
					coverage(alignments)[[contig]]
				},
				error=function(e) {
					alignments <- readGAlignments(alignmentsFile, param=ScanBamParam(which=GRanges(addChr(contig), coverageRegion)))
					coverage(alignments)[[addChr(contig)]]
				}
			)
			if (exists("alignments")) rm(alignments)
			return(coverageData)
		}
		# get coverage track
		coverage1 <- readCoverage(alignmentsFile, SVs[SV,"contig1"], coverageRegion1)
		coverage2 <- readCoverage(alignmentsFile, SVs[SV,"contig2"], coverageRegion2)
		# shrink coverage range to chromosome boundaries to avoid subscript out of bounds errors
		coverageRegion1 <- IRanges(max(start(coverageRegion1), min(start(coverage1))), min(end(coverageRegion1), max(end(coverage1))))
		coverageRegion2 <- IRanges(max(start(coverageRegion2), min(start(coverage2))), min(end(coverageRegion2), max(end(coverage2))))
	}

	# find all exons belonging to the fused genes
	exons1 <- findExons(exons, SVs[SV,"contig1"], SVs[SV,"gene_id1"], SVs[SV,"direction1"], SVs[SV,"Breakpoint1"], coverage1, SVs[SV,"transcript_id1"], transcriptSelection)
	if (nrow(exons1) == 0) {
		par(mfrow=c(1,1))
		plot(0, 0, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
		text(0, 0, paste0("exon coordinates of ", SVs[SV,"gene1"], " not found in\n", exonsFile))
		warning(paste("exon coordinates of", SVs[SV,"gene1"], "not found"))
		next
	}
	exons2 <- findExons(exons, SVs[SV,"contig2"], SVs[SV,"gene_id2"], SVs[SV,"direction2"], SVs[SV,"Breakpoint2"], coverage2, SVs[SV,"transcript_id2"], transcriptSelection)
	if (nrow(exons2) == 0) {
		par(mfrow=c(1,1))
		plot(0, 0, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
		text(0, 0, paste0("exon coordinates of ", SVs[SV,"gene2"], " not found in\n", exonsFile))
		warning(paste("exon coordinates of", SVs[SV,"gene2"], "not found"))
		next
	}

	# in case of intergenic breakpoints, show the vicinity
	if (sum(showVicinity) > 0) {
		if (SVs[SV,"site1"] == "intergenic") {
			for (geneID in unique(exons[exons$contig == SVs[SV,"contig1"] & exons$exonNumber != "intergenic" &
			                            (between(exons$end, SVs[SV,"Breakpoint1"]-showVicinity[1], SVs[SV,"Breakpoint1"]+showVicinity[2]) |
			                             between(exons$start, SVs[SV,"Breakpoint1"]-showVicinity[1], SVs[SV,"Breakpoint1"]+showVicinity[2])),"geneID"]))
				exons1 <- rbind(exons1, findExons(exons, SVs[SV,"contig1"], geneID, SVs[SV,"direction1"], SVs[SV,"Breakpoint1"], coverage1, SVs[SV,"transcript_id1"], transcriptSelection))
			# crop genes that are only partially within user-defined vicinity, because the coverage data is incomplete for those
			exons1 <- exons1[exons1$start >= SVs[SV,"Breakpoint1"]-showVicinity[1] & exons1$end <= SVs[SV,"Breakpoint1"]+showVicinity[2] | exons1$exonNumber == "intergenic",]
		}
		if (SVs[SV,"site2"] == "intergenic") {
			for (geneID in unique(exons[exons$contig == SVs[SV,"contig2"] & exons$exonNumber != "intergenic" &
			                            (between(exons$end, SVs[SV,"Breakpoint2"]-showVicinity[3], SVs[SV,"Breakpoint2"]+showVicinity[4]) |
			                             between(exons$start, SVs[SV,"Breakpoint2"]-showVicinity[3], SVs[SV,"Breakpoint2"]+showVicinity[4])),"geneID"]))
				exons2 <- rbind(exons2, findExons(exons, SVs[SV,"contig2"], geneID, SVs[SV,"direction2"], SVs[SV,"Breakpoint2"], coverage2, SVs[SV,"transcript_id2"], transcriptSelection))
			# crop genes that are only partially within user-defined vicinity, because the coverage data is incomplete for those
			exons2 <- exons2[exons2$start >= SVs[SV,"Breakpoint2"]-showVicinity[3] & exons2$end <= SVs[SV,"Breakpoint2"]+showVicinity[4] | exons2$exonNumber == "intergenic",]
		}
	}

	# normalize coverage
	if (alignmentsFile != "") {
		coverageNormalization <- function(coverage, coverageRegion, exons) {
			max(1, ifelse(
				squishIntrons, # => ignore intronic coverage
				max(as.numeric(coverage[IRanges(sapply(exons$start,max,min(start(coverage))),sapply(exons$end,min,max(end(coverage))))])),
				round(quantile(coverage[coverageRegion], 0.9999)) # ignore coverage spikes from read-attracting regions
			))
		}
		coverageNormalization1 <- ifelse(head(coverageRange,1) == 0, coverageNormalization(coverage1, coverageRegion1, exons1), head(coverageRange,1))
		coverageNormalization2 <- ifelse(tail(coverageRange,1) == 0, coverageNormalization(coverage2, coverageRegion2, exons2), tail(coverageRange,1))
		if (length(coverageRange) == 1 && coverageRange == 0) { # harmonize scales of gene1 and gene2
			coverageNormalization1 <- max(coverageNormalization1, coverageNormalization2)
			coverageNormalization2 <- max(coverageNormalization1, coverageNormalization2)
		}
		coverage1 <- coverage1/coverageNormalization1
		coverage2 <- coverage2/coverageNormalization2
		coverage1[coverage1 > 1] <- 1
		coverage2[coverage2 > 1] <- 1
	}

	# sort coding exons last, such that they are drawn over the border of non-coding exons
	exons1 <- exons1[order(exons1$start, -rank(exons1$type)),]
	exons2 <- exons2[order(exons2$start, -rank(exons2$type)),]

	# insert dummy exons, if breakpoints are outside the gene (e.g., in UTRs)
	# this avoids plotting artifacts
	Breakpoint1 <- SVs[SV,"Breakpoint1"]
	Breakpoint2 <- SVs[SV,"Breakpoint2"]
	if (Breakpoint1 < min(exons1$start)) {
		exons1 <- rbind(c(exons1[1,"contig"], "dummy", max(1,Breakpoint1-1000), max(1,Breakpoint1-1000), exons1[1,"strand"], "", "dummy", exons1[1,"geneID"], exons1[1,"transcript"], ""), exons1)
	} else if (Breakpoint1 > max(exons1$end)) {
		exons1 <- rbind(exons1, c(exons1[1,"contig"], "dummy", Breakpoint1+1000, Breakpoint1+1000, exons1[1,"strand"], "", "dummy", exons1[1,"geneID"], exons1[1,"transcript"], ""))
	}
	if (Breakpoint2 < min(exons2$start)) {
		exons2 <- rbind(c(exons2[1,"contig"], "dummy", max(1,Breakpoint2-1000), max(1,Breakpoint2-1000), exons2[1,"strand"], "", "dummy", exons2[1,"geneID"], exons2[1,"transcript"], ""), exons2)
	} else if (Breakpoint2 > max(exons2$end)) {
		exons2 <- rbind(exons2, c(exons2[1,"contig"], "dummy", Breakpoint2+1000, Breakpoint2+1000, exons2[1,"strand"], "", "dummy", exons2[1,"geneID"], exons2[1,"transcript"], ""))
	}
	exons1$start <- as.integer(exons1$start)
	exons1$end <- as.integer(exons1$end)
	exons2$start <- as.integer(exons2$start)
	exons2$end <- as.integer(exons2$end)

	exons1$left <- exons1$start
	exons1$right <- exons1$end
	exons2$left <- exons2$start
	exons2$right <- exons2$end

	squishedIntronSize <- 200
	if (squishIntrons) {
		# hide introns in gene1
		cumulativeIntronLength <- 0
		previousExonEnd <- -squishedIntronSize
		for (exon in 1:nrow(exons1)) {
			if (Breakpoint1 > previousExonEnd+1 && Breakpoint1 < exons1[exon,"left"])
				Breakpoint1 <- (Breakpoint1-previousExonEnd) / (exons1[exon,"left"]-previousExonEnd) * squishedIntronSize + previousExonEnd - cumulativeIntronLength
			if (exons1[exon,"left"] > previousExonEnd) {
				cumulativeIntronLength <- cumulativeIntronLength + exons1[exon,"left"] - previousExonEnd - squishedIntronSize
				previousExonEnd <- exons1[exon,"right"]
			}
			if (Breakpoint1 >= exons1[exon,"left"] && Breakpoint1 <= exons1[exon,"right"]+1)
				Breakpoint1 <- Breakpoint1 - cumulativeIntronLength
			exons1[exon,"left"] <- exons1[exon,"left"] - cumulativeIntronLength
			exons1[exon,"right"] <- exons1[exon,"right"] - cumulativeIntronLength
		}

		# hide introns in gene2
		cumulativeIntronLength <- 0
		previousExonEnd <- -squishedIntronSize
		for (exon in 1:nrow(exons2)) {
			if (Breakpoint2 > previousExonEnd+1 && Breakpoint2 < exons2[exon,"left"])
				Breakpoint2 <- (Breakpoint2-previousExonEnd) / (exons2[exon,"left"]-previousExonEnd) * squishedIntronSize + previousExonEnd - cumulativeIntronLength
			if (exons2[exon,"left"] > previousExonEnd) {
				cumulativeIntronLength <- cumulativeIntronLength + exons2[exon,"left"] - previousExonEnd - squishedIntronSize
				previousExonEnd <- exons2[exon,"right"]
			}
			if (Breakpoint2 >= exons2[exon,"left"] && Breakpoint2 <= exons2[exon,"right"]+1)
				Breakpoint2 <- Breakpoint2 - cumulativeIntronLength
			exons2[exon,"left"] <- exons2[exon,"left"] - cumulativeIntronLength
			exons2[exon,"right"] <- exons2[exon,"right"] - cumulativeIntronLength
		}
	} else { # don't squish introns
		# shift exon coordinates to align the gene to the left border of the plot
		exons1$right <- exons1$right - min(exons1$left)
		Breakpoint1 <- Breakpoint1 - min(exons1$left)
		exons1$left <- exons1$left - min(exons1$left)
		exons2$right <- exons2$right - min(exons2$left)
		Breakpoint2 <- Breakpoint2 - min(exons2$left)
		exons2$left <- exons2$left - min(exons2$left)
	}

	# scale exon sizes to fit on page
	scalingFactor <- max(exons1$right) + max(exons2$right)
	if (fixedScale > 0) {
		if (fixedScale >= scalingFactor) {
			scalingFactor <- fixedScale
		} else {
			warning(paste("fallback to automatic scaling, because value for --fixedScale is too small to fit transcripts on canvas (increase it to", scalingFactor, "to avoid this)"))
		}
	}
	exons1$left <- exons1$left / scalingFactor
	exons1$right <- exons1$right / scalingFactor
	exons2$left <- exons2$left / scalingFactor
	exons2$right <- exons2$right / scalingFactor
	Breakpoint1 <- Breakpoint1 / scalingFactor
	Breakpoint2 <- Breakpoint2 / scalingFactor

	# shift gene2 to the right border of the page
	gene2Offset <- 1 + 0.05 - max(exons2$right)

	# center SV horizontally
	SVOffset1 <- (max(exons1$right)+gene2Offset)/2 - ifelse(SVs[SV,"direction1"] == "downstream", Breakpoint1, max(exons1$right)-Breakpoint1)
	SVOffset2 <- SVOffset1 + ifelse(SVs[SV,"direction1"] == "downstream", Breakpoint1, max(exons1$right)-Breakpoint1)

	# layout: SV on top, circos plot on bottom left, protein domains on bottom center, statistics on bottom right
	layout(matrix(c(1,1,1,2,4,5,3,4,5), 3, 3, byrow=TRUE), widths=c(1.1, 1.2, 0.7), heights=c(1.55, 1.2, 0.25))
	par(mar=c(0, 0, 0, 0))
	plot(0, 0, type="l", xlim=c(-0.12, 1.12), ylim=c(0.4, 1.1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")

	# vertical coordinates of layers
	ySampleName <- 1.04
	yIdeograms <- ifelse(alignmentsFile != "", 0.94, 0.84)
	yBreakpointLabels <- ifelse(alignmentsFile != "", 0.86, 0.76)
	yCoverage <- 0.72
	yExons <- 0.67
	yGeneNames <- 0.58
	ySV <- 0.5
	yTranscript <- 0.45
	yScale <- 0.407
	yTrajectoryBreakpointLabels <- yBreakpointLabels - 0.035
	yTrajectoryExonTop <- yExons + 0.03
	yTrajectoryExonBottom <- yExons - 0.055
	yTrajectorySV <- ySV + 0.03

	# print sample name (title of page)
	text(0.5, ySampleName, sampleName, font=2, cex=fontSize*1.5, adj=c(0.5,0))

	# draw ideograms
	if (!is.null(cytobands)) {
		drawIdeogram("left", min(exons1$left), max(exons1$right), yIdeograms, cytobands, SVs[SV,"contig1"], SVs[SV,"Breakpoint1"])
		drawIdeogram("right", gene2Offset, gene2Offset+max(exons2$right), yIdeograms, cytobands, SVs[SV,"contig2"], SVs[SV,"Breakpoint2"])
	}

	# draw gene & transcript names
	if (SVs[SV,"gene1"] != ".")
		text(max(exons1$right)/2, yGeneNames, SVs[SV,"gene1"], font=4, cex=fontSize, adj=c(0.5,0))
	if (SVs[SV,"site1"] != "intergenic")
		text(max(exons1$right)/2, yGeneNames-0.01, head(exons1$transcript,1), cex=0.9*fontSize, adj=c(0.5,1))
	if (SVs[SV,"gene2"] != ".")
		text(gene2Offset+max(exons2$right)/2, yGeneNames, SVs[SV,"gene2"], font=4, cex=fontSize, adj=c(0.5,0))
	if (SVs[SV,"site2"] != "intergenic")
		text(gene2Offset+max(exons2$right)/2, yGeneNames-0.01, head(exons2$transcript,1), cex=0.9*fontSize, adj=c(0.5,1))

	# if multiple genes in the vicinity are shown, label them
	if (SVs[SV,"site1"] == "intergenic")
		for (gene in unique(exons1$geneName)) {
			exonsOfGene <- exons1[exons1$geneName == gene & exons1$type != "dummy",]
			if (any(exonsOfGene$type == "exon"))
				text(mean(c(min(exonsOfGene$left), max(exonsOfGene$right))), yExons-0.04, gene, cex=0.9*fontSize, adj=c(0.5,1))
		}
	if (SVs[SV,"site2"] == "intergenic")
		for (gene in unique(exons2$geneName)) {
			exonsOfGene <- exons2[exons2$geneName == gene & exons2$type != "dummy",]
			if (any(exonsOfGene$type == "exon"))
				text(gene2Offset+mean(c(min(exonsOfGene$left), max(exonsOfGene$right))), yExons-0.04, gene, cex=0.9*fontSize, adj=c(0.5,1))
		}

	# label breakpoints
	text(Breakpoint1+0.01, yBreakpointLabels-0.03, paste0("Breakpoint1\n", SVs[SV,"display_contig1"], ":", SVs[SV,"Breakpoint1"]), adj=c(1,0), cex=fontSize)
	text(gene2Offset+Breakpoint2-0.01, yBreakpointLabels-0.03, paste0("Breakpoint2\n", SVs[SV,"display_contig2"], ":", SVs[SV,"Breakpoint2"]), adj=c(0,0), cex=fontSize)

	# draw coverage axis
	if (alignmentsFile != "") {
		# left axis (gene1)
		lines(c(-0.02, -0.01, -0.01, -0.02), c(yCoverage, yCoverage, yCoverage+0.1, yCoverage+0.1))
		text(-0.025, yCoverage, "0", adj=c(1,0.5), cex=0.9*fontSize)
		text(-0.025, yCoverage+0.1, coverageNormalization1, adj=c(1,0.5), cex=0.9*fontSize)
		text(-0.05, yCoverage+0.08, "Coverage", srt=90, cex=0.9*fontSize, adj=c(1,0.5))

		# right axis (gene2)
		if (length(coverageRange) == 2) { # separate axes for gene1 and gene2
			rightCoverageAxisX <- gene2Offset+max(exons2$right)
			lines(c(rightCoverageAxisX+0.02, rightCoverageAxisX+0.01, rightCoverageAxisX+0.01, rightCoverageAxisX+0.02), c(yCoverage, yCoverage, yCoverage+0.1, yCoverage+0.1))
			text(rightCoverageAxisX+0.025, yCoverage, "0", adj=c(0,0.5), cex=0.9*fontSize)
			text(rightCoverageAxisX+0.025, yCoverage+0.1, coverageNormalization2, adj=c(0,0.5), cex=0.9*fontSize)
			text(rightCoverageAxisX+0.05, yCoverage+0.08, "Coverage", srt=90, cex=0.9*fontSize, adj=c(1,0.5))
		}

		# plot coverage 1
		rect(min(exons1$left), yCoverage, max(exons1$right), yCoverage+0.1, col="#eeeeee", border=NA)
		if (squishIntrons) {
			for (exon in 1:nrow(exons1))
				if (exons1[exon,"type"] != "CDS") # don't draw coverage twice for coding regions
					drawCoverage(exons1[exon,"left"], exons1[exon,"right"], yCoverage, coverage1, exons1[exon,"start"], exons1[exon,"end"], color1)
		} else {
			drawCoverage(min(exons1$left), max(exons1$right), yCoverage, coverage1, min(exons1$start), max(exons1$end), color1)
		}

		# plot coverage 2
		rect(gene2Offset+min(exons2$left), yCoverage, gene2Offset+max(exons2$right), yCoverage+0.1, col="#eeeeee", border=NA)
		if (squishIntrons) {
			for (exon in 1:nrow(exons2))
				if (exons2[exon,"type"] != "CDS") # don't draw coverage twice for coding regions
					drawCoverage(gene2Offset+exons2[exon,"left"], gene2Offset+exons2[exon,"right"], yCoverage, coverage2, exons2[exon,"start"], exons2[exon,"end"], color2)
		} else {
			drawCoverage(gene2Offset+min(exons2$left), gene2Offset+max(exons2$right), yCoverage, coverage2, min(exons2$start), max(exons2$end), color2)
		}
	}

	# plot gene 1
	lines(c(min(exons1$left), max(exons1$right)), c(yExons, yExons), col=darkColor1)
	for (gene in unique(exons1$geneName))
		drawStrand(min(exons1[exons1$geneName == gene,"left"]), max(exons1[exons1$geneName == gene,"right"]), yExons, darkColor1, head(exons1[exons1$geneName == gene,"strand"],1))
	for (exon in 1:nrow(exons1))
		drawExon(exons1[exon,"left"], exons1[exon,"right"], yExons, color1, exons1[exon,"exonNumber"], exons1[exon,"type"])

	# plot gene 2
	lines(c(gene2Offset, gene2Offset+max(exons2$right)), c(yExons, yExons), col=darkColor2)
	for (gene in unique(exons2$geneName))
		drawStrand(gene2Offset+min(exons2[exons2$geneName == gene,"left"]), gene2Offset+max(exons2[exons2$geneName == gene,"right"]), yExons, darkColor2, head(exons2[exons2$geneName == gene,"strand"],1))
	for (exon in 1:nrow(exons2))
		drawExon(gene2Offset+exons2[exon,"left"], gene2Offset+exons2[exon,"right"], yExons, color2, exons2[exon,"exonNumber"], exons2[exon,"type"])

	# plot gene1 of SV
	if (SVs[SV,"direction1"] == "downstream") {
		# plot strands
		lines(c(SVOffset1, SVOffset1+Breakpoint1), c(ySV, ySV), col=darkColor1)
		for (gene in unique(exons1$geneName)) {
			exonsOfGene <- exons1[exons1$geneName == gene,]
			if (min(exonsOfGene$start) <= SVs[SV,"Breakpoint1"])
				drawStrand(SVOffset1+min(exonsOfGene$left), SVOffset1+min(Breakpoint1, max(exonsOfGene$right)), ySV, col=darkColor1, exonsOfGene$strand[1])
		}
		# plot exons
		for (exon in 1:nrow(exons1))
			if (exons1[exon,"start"] <= SVs[SV,"Breakpoint1"])
				drawExon(SVOffset1+exons1[exon,"left"], SVOffset1+min(Breakpoint1, exons1[exon,"right"]), ySV, color1, exons1[exon,"exonNumber"], exons1[exon,"type"])
		# plot trajectories
		lines(c(0, 0, SVOffset1), c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectorySV), col="red", lty=2)
		lines(c(Breakpoint1, Breakpoint1, SVOffset1+Breakpoint1), c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectorySV), col="red", lty=2)
	} else if (SVs[SV,"direction1"] == "upstream") {
		# plot strands
		lines(c(SVOffset1, SVOffset2), c(ySV, ySV), col=darkColor1)
		for (gene in unique(exons1$geneName)) {
			exonsOfGene <- exons1[exons1$geneName == gene,]
			if (max(exonsOfGene$end+1) >= SVs[SV,"Breakpoint1"])
				drawStrand(SVOffset2-max(exonsOfGene$right)+Breakpoint1, min(SVOffset2, SVOffset2-min(exonsOfGene$left)+Breakpoint1), ySV, col=darkColor1, chartr("+-", "-+", exonsOfGene$strand[1]))
		}
		# plot exons
		for (exon in 1:nrow(exons1))
			if (exons1[exon,"end"]+1 >= SVs[SV,"Breakpoint1"])
				drawExon(SVOffset1+max(exons1$right)-exons1[exon,"right"], min(SVOffset2, SVOffset1+max(exons1$right)-exons1[exon,"left"]), ySV, color1, exons1[exon,"exonNumber"], exons1[exon,"type"])
		# plot trajectories
		lines(c(max(exons1$right), max(exons1$right), SVOffset1), c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectorySV), col="red", lty=2)
		lines(c(Breakpoint1, Breakpoint1, SVOffset1+max(exons1$right)-Breakpoint1), c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectorySV), col="red", lty=2)
	}
	
	# plot gene2 of SV
	if (SVs[SV,"direction2"] == "downstream") {
		# plot strands
		lines(c(SVOffset2, SVOffset2+Breakpoint2), c(ySV, ySV), col=darkColor2)
		for (gene in unique(exons2$geneName)) {
			exonsOfGene <- exons2[exons2$geneName == gene,]
			if (min(exonsOfGene$start) <= SVs[SV,"Breakpoint2"])
				drawStrand(max(SVOffset2, SVOffset2+Breakpoint2-max(exonsOfGene$right)), SVOffset2+Breakpoint2-min(exonsOfGene$left), ySV, col=darkColor2, chartr("+-", "-+", exonsOfGene$strand[1]))
		}
		# plot exons
		for (exon in 1:nrow(exons2))
			if (exons2[exon,"start"] <= SVs[SV,"Breakpoint2"])
				drawExon(max(SVOffset2, SVOffset2+Breakpoint2-exons2[exon,"right"]), SVOffset2+Breakpoint2-exons2[exon,"left"], ySV, color2, exons2[exon,"exonNumber"], exons2[exon,"type"])
		# plot trajectories
		lines(c(gene2Offset, gene2Offset, SVOffset2+Breakpoint2), c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectorySV), col="red", lty=2)
		lines(c(gene2Offset+Breakpoint2, gene2Offset+Breakpoint2, SVOffset2), c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectorySV), col="red", lty=2)
	} else if (SVs[SV,"direction2"] == "upstream") {
		# plot strands
		lines(c(SVOffset2, SVOffset2+max(exons2$right)-Breakpoint2), c(ySV, ySV), col=darkColor2)
		for (gene in unique(exons2$geneName)) {
			exonsOfGene <- exons2[exons2$geneName == gene,]
			if (max(exonsOfGene$end+1) >= SVs[SV,"Breakpoint2"])
				drawStrand(max(SVOffset2, SVOffset2+min(exonsOfGene$left)-Breakpoint2), SVOffset2+max(exonsOfGene$right)-Breakpoint2, ySV, col=darkColor2, exonsOfGene$strand[1])
		}
		# plot exons
		for (exon in 1:nrow(exons2))
			if (exons2[exon,"end"]+1 >= SVs[SV,"Breakpoint2"])
				drawExon(max(SVOffset2, SVOffset2+exons2[exon,"left"]-Breakpoint2), SVOffset2+exons2[exon,"right"]-Breakpoint2, ySV, color2, exons2[exon,"exonNumber"], exons2[exon,"type"])
		# plot trajectories
		lines(c(gene2Offset+max(exons2$right), gene2Offset+max(exons2$right), SVOffset2+max(exons2$right)-Breakpoint2), c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectorySV), col="red", lty=2)
		lines(c(gene2Offset+Breakpoint2, gene2Offset+Breakpoint2, SVOffset2), c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectorySV), col="red", lty=2)
	}
	
	if (SVs[SV,"SV_transcript"] != ".") {
		# print SV transcript colored by gene of origin
		SV_transcript1 <- gsub("\\|.*", "", SVs[SV,"SV_transcript"], perl = TRUE)
		SV_transcript1 <- substr(SV_transcript1, max(1, nchar(SV_transcript1)-30), nchar(SV_transcript1))
		SV_transcript2 <- gsub(".*\\|", "", SVs[SV,"SV_transcript"], perl = TRUE)
		SV_transcript2 <- substr(SV_transcript2, 1, min(nchar(SV_transcript2), 30))
		# check for non-template bases
		non_template_bases <- gsub(".*\\|([^|]*)\\|.*", "\\1", SVs[SV,"SV_transcript"], perl = TRUE)
		if (non_template_bases == SVs[SV,"SV_transcript"]) # no non-template bases found
			non_template_bases <- ""
		# divide non-template bases half-and-half for centered alignment
		non_template_bases1 <- substr(non_template_bases, 1, floor(nchar(non_template_bases)/2))
		non_template_bases2 <- substr(non_template_bases, ceiling(nchar(non_template_bases)/2), nchar(non_template_bases))
		# transcript 1
		text(SVOffset2, yTranscript, bquote(.(SV_transcript1) * phantom(.(non_template_bases1))), col=darkColor1, adj=c(1,0.5), cex=fontSize)
		# transcript 2
		text(SVOffset2, yTranscript, bquote(phantom(.(non_template_bases2)) * .(SV_transcript2)), col=darkColor2, adj=c(0,0.5), cex=fontSize)
		# non-template bases
		text(SVOffset2, yTranscript, non_template_bases1, adj=c(1,0.5), cex=fontSize)
		text(SVOffset2, yTranscript, non_template_bases2, adj=c(0,0.5), cex=fontSize)
	}

	# draw scale
	realScale <- max(exons1$end - exons1$start, exons2$end - exons2$start)
	mapScale <- max(exons1$right - exons1$left, exons2$right - exons2$left)
	# choose scale which is closest to desired scale length
	desiredScaleSize <- 0.2
	realScale <- desiredScaleSize / mapScale * realScale
	mapScale <- desiredScaleSize
	realScaleOptimalFit <- signif(realScale, 1) # round to most significant digit
	mapScaleOptimalFit <- realScaleOptimalFit / realScale * mapScale
	# draw scale line
	lines(c(1-mapScaleOptimalFit, 1), c(yScale, yScale)) # scale line
	lines(c(1-mapScaleOptimalFit, 1-mapScaleOptimalFit), c(yScale-0.007, yScale+0.007)) # left whisker
	lines(c(1, 1), c(yScale-0.007, yScale+0.007)) # right whisker
	# draw units above scale line
	realScaleThousands <- max(0, min(3, floor(log10(realScaleOptimalFit)/3)))
	scaleUnits <- c("bp", "kbp", "Mbp", "Gbp")
	scaleLabel <- paste(realScaleOptimalFit/max(1,1000^realScaleThousands), scaleUnits[realScaleThousands+1])
	text(1-mapScaleOptimalFit/2, yScale+0.005, scaleLabel, adj=c(0.5,0), cex=fontSize*0.9)
	if (squishIntrons)
		text(1-mapScaleOptimalFit/2, yScale-0.005, "Introns not to scale", adj=c(0.5,1), cex=fontSize*0.9, font=3)

	# draw circos plot
	if (is.null(cytobands) || !("circlize" %in% names(sessionInfo()$otherPkgs)) || !("GenomicRanges" %in% names(sessionInfo()$otherPkgs))) {
		plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
		plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
	} else {
		par(mar=c(2,4,0,0), xpd=NA)
		drawCircos(SV, SVs, cytobands, minConfidenceForCircosPlot, circosColors)
		par(mar=c(0,0,0,0), xpd=F)
	}

	# draw protein domains
	plot(0, 0, type="l", xlim=c(-0.1, 1.1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
	par(xpd=NA)
	if (!is.null(proteinDomains))
		drawProteinDomains(SVs[SV,], exons1, exons2, proteinDomains, color1, color2, mergeDomainsOverlappingBy, optimizeDomainColors)
	par(xpd=F)

	# print statistics about supporting alignments
	plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
	text(0, 0.575, "ID & SUPPORTING READ COUNT", font=2, adj=c(0,0), cex=fontSize)
	if ("split_reads" %in% colnames(SVs)) { # STAR-SV reports split reads from both breakpoints combined
		text(0, 0.525, paste0("Split reads = ", SVs[SV,"split_reads"], "\n", "Discordant mates = ", SVs[SV,"discordant_mates"]), adj=c(0,1), cex=fontSize)
	} else { # Arriba reports split reads separately for the two breakpoints
		text(
			0, 0.525,
			paste0(
				"Structural variant ID = ", SVs[SV,"SV_ID"], "\n",
				"Number of reads supporting SV = ", SVs[SV,"coverage1"]
			),
			adj=c(0,1), cex=fontSize
		)
	}

}

devNull <- dev.off()
message("Done")
