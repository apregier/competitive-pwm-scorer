#!/usr/bin/env Rscript

library(rtfbs)
library(getopt)

competitiveScoring <- function(pwmFile, germlineFile, somaticFile) {
    pwm <- read.pwm(pwmFile)
    maxScores.germline <- scoreOneFile(pwm, germlineFile)
    maxScores.somatic <- scoreOneFile(pwm, somaticFile)
    combinedMaxScores <- combineMaxScores(maxScores.germline, maxScores.somatic)
    combinedMaxScores
}

scoreOneFile <- function(pwms, fastaFile) {
    ms <- read.ms(fastaFile)
    mm <- build.mm(ms, 3)
    initialNames <- matrix(nrow=length(names.ms(ms)))
    initialScores <- matrix(nrow=length(names.ms(ms)))
    maxScores <- data.frame(pwmName= initialNames, score= initialScores, row.names = names.ms(ms))
    for (i in 1:length(pwms)) {
        scores <-scoreOnePwm(pwms[[i]], ms, mm)
        maxScores <- updateMaxScores(maxScores, scores)
    }
    maxScores
}

scoreOnePwm <- function(pwm, ms, mm) {
    cs <- score.ms(ms, pwm, mm)
    cs
    #v <- simulate.ms(mm, 100000)
    #xs <- score.ms(v, pwm, mm, threshold=-2)
    #fdr <- calc.fdr(ms, cs, v, xs)
}

#maxScores:
#seqName pwmName score
updateMaxScores <- function(maxScores, scores) {
    for (name in scores$seqname) {
        if (!is.na(scores[name, "score"])) {
            if (is.na(maxScores[name, "score"]) || scores[name, "score"] > maxScores[name, "score"]) {
                maxScores[name,"score"] <- scores[name,"score"]
                    maxScores[name,"pwmName"] <- scores[name, "src"]
            }
        }
    }
    maxScores
}

#combinedMaxScores:
#seqName germline.pwmName germline.pwmScore somatic.pwmName somatic.pwmScore tfbsChange scoreChange
combineMaxScores <- function(maxScores.germline, maxScores.somatic) {
    combinedMaxScores <- data.frame(row.names = names(maxScores.germline))
    for (name in combinedMaxScores$name) {
        combinedMaxScores[name, "germline.pwmName"] <- maxScores.germline[name, "pwmName"]
        combinedMaxScores[name, "germline.score"] <- maxScores.germline[name, "score"]
        combinedMaxScores[name, "somatic.pwmName"] <- maxScores.somatic[name, "pwmName"]
        combinedMaxScores[name, "somatic.score"] <- maxScores.somatic[name, "score"]
        if (maxScores.germline[name, "pwmName"] != maxScores.somatic[name, "pwmName"]) {
            combinedMaxScores[name, "tfbsChange"] <- TRUE
        }
        else {
            combinedMaxScores[name, "tfbsChange"] <- FALSE
        }
        combinedMaxScores[name, "scoreChange"] <- combinedMaxScores[name, "somatic.score"] - combinedMaxScores[name, "germline.score"]
    }
    combinedMaxScores
}

optSpec = matrix(c(
    "help", "h", 0, NA,
    "pwmFile", "p", 1, "character",
    "germlineFile", "g", 1, "character",
    "somaticFile", "s", 1, "character"
    ), byrow=TRUE, ncol=4)

opts <- getopt(optSpec)
competition <- competitiveScoring(opts$`pwmFile`, opts$`germlineFile`, opts$`somaticFile`)
competition
