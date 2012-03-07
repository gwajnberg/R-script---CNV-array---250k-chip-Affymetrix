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

## CGHcall + DNAcopy
library(CGHcall)
raw <- data.frame(SNP=snpInfo[['man_fsetid']],
                  CHROMOSOME=snpInfo[['chrom']],
                  START_POS=snpInfo[['physical_pos']],
                  END_POS=snpInfo[['physical_pos']]+1L,
                  logRatioSense)
raw$CHROMOSOME <- as.character(raw$CHROMOSOME)
raw <- raw[complete.cases(raw),]
rawAutosomes <- subset(raw, CHROMOSOME %in% as.character(1:22))
rm(raw)
rawAutosomes$CHROMOSOME <- as.integer(rawAutosomes$CHROMOSOME)
idx <- with(rawAutosomes, order(CHROMOSOME, START_POS))
rawAutosomes <- rawAutosomes[idx,]
rm(idx)
rawAutosomes <- cghRaw(rawAutosomes)
prep <- preprocess(rawAutosomes)
norm <- normalize(prep)
segm <- segmentData(norm)
post <- postsegnormalize(segm)
final <- CGHcall(post, nclass=3, prior='all')
final <- ExpandCGHcall(final, segm)
