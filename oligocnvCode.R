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
