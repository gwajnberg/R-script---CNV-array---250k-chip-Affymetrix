library(oligo)
## Obtendo os arquivos
fullFilenames <- list.celfiles(path="/storage/data/elf/MBra/teste",
                               full.names=TRUE)
outputDir     <- file.path(getwd(), "crlmmTest")

## Genotipando
crlmm(fullFilenames, outputDir, verbose=FALSE)
crlmmOut <- getCrlmmSummaries(outputDir)

## Intensidade total por SNP x Strand
theS <- getA(crlmmOut) * 2
grps <- rep(c('normal', 'tumor'), each=7)
colnames(theS) <- paste(grps, 1:7, sep='')

## Acessar anotacao e obter:
##   a) Nome do SNP
##   b) Cromosomo
##   c) Posicao fisica
pkgname <- annotation(crlmmOut)
library(pkgname, character.only=TRUE)
conn <- db(get(pkgname))
sql <- paste('SELECT man_fsetid, chrom, physical_pos',
             'FROM featureSet WHERE man_fsetid LIKE "SNP%"')
snpInfo <- dbGetQuery(conn, sql)
i <- match(snpInfo[['man_fsetid']], rownames(theS))
theS <- theS[i,,]

## Confirme que a order dos SNPs em theS eh a mesma que em theS
stopifnot(identical(snpInfo[['man_fsetid']],
                    rownames(theS)))


logRatioSense <- theS[, 8:14, 2]-theS[, 1:7, 2]
logRatioAntisense <- theS[, 8:14, 1]-theS[, 1:7, 1]

## Parte do DNAcopy
## - Alisamento
## - Segmentacao
library(DNAcopy)
cnObjSense <- CNA(logRatioSense,
                  chrom=snpInfo[["chrom"]],
                  maploc=snpInfo[["physical_pos"]],
                  data.type="logratio")
smoothedCNASense <- smooth.CNA(cnObjSense)
segment.smoothed.CNA.object <- segment(smoothedCNASense)
sdundo.CNA.object <- segment(smoothedCNASense,
                             undo.splits="sdundo",
                             undo.SD=3, verbose=1)

pdf("mainResults.pdf")
plot(segment.smoothed.CNA.object, plot.type='s')
plot(segment.smoothed.CNA.object, plot.type='w')
plot(sdundo.CNA.object, plot.type='s')
dev.off()

madThreshold <- function (xout, threshold = 1){
    if (!inherits(xout, "DNAcopy"))
        stop("First arg must be of class DNAcopy")
    nsample <- ncol(xout$data) - 2
    snames <- names(xout$data)
    status <- matrix(0L, nr=nrow(xout$data), nc=nsample)
    colnames(status) <- snames[-(1:2)]
    for (i in 1:nsample){
        gain <- loss <- rep(FALSE, nr=nrow(xout$data))
        j <- i+2L
        sout <- xout$output[xout$output$ID == snames[j], ]
        xmad <- mad(na.omit(xout$data[, j]) - rep(sout$seg.mean, sout$num.mark))

        genomdat <- xout$data[, j]
        ii <- which(is.finite(genomdat))
        segout <- xout$output[xout$output$ID == snames[j],]
        segmean <- rep(segout$seg.mean, segout$num.mark)
        stat <- (segmean - median(segmean))/xmad
        gain[ii] <- stat > threshold
        loss[ii] <- stat < -threshold
        status[,i] <- gain-loss

    }
    out <- list(chrom=xout$data$chrom,
                maploc=xout$data$maploc,
                status)
    as.data.frame(out)
}

LossGain <- madThreshold(sdundo.CNA.object, 1)
