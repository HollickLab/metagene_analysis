# multiple_t_test.R contains the function (multi.t.test) which will perform
# pairwise Welch's t-tests for each window, apply a Holm-Bonferoni correction
# (alpha = 0.05) and plot the sum of the metagene data along with marks of the
# significantly different windows.
# Please see README for full details and examples.
#
# Requires:
#     R (http://cran.us.r-project.org/)
# 
# Joy-El R.B. Talbot Copyright (c) 2014
# 
# The MIT License (MIT)
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

args <- commandArgs(trailingOnly = TRUE)

multi.t.test <- function (file.1.sense,
                          file.1.antisense,
                          file.2.sense, 
                          file.2.antisense,
                          normalization.1, 
                          normalization.2, 
                          output.prefix, 
                          window.size, 
                          window.step,
                          feature.name) {
    options(scipen=9999) # keep things out of scientific notation
    million = 1000000
    
    # load data
    windows.1.sense <- read.table(file.1.sense, header=TRUE, sep=",")
    windows.1.sense$Window <- factor(windows.1.sense$Window)
    total.windows <- length(levels(windows.1.sense$Window))
    
    windows.1.antisense <- read.table(file.1.antisense, header=TRUE, sep=",")
    windows.1.antisense$Window <- factor(windows.1.antisense$Window)
    
    windows.2.sense <- read.table(file.2.sense, header=TRUE, sep=",")
    windows.2.sense$Window <- factor(windows.2.sense$Window)
    
    windows.2.antisense <- read.table(file.2.antisense, header=TRUE, sep=",")
    windows.2.antisense$Window <- factor(windows.2.antisense$Window)

    # create summary file for plots
    windows.sum <- aggregate(windows.1.sense$Abundance, by=list(windows.1.sense$Window), FUN=sum)
    names(windows.sum) <- c("Window", "windows.1.sense")
    
    windows.sum$windows.2.sense <- aggregate(windows.2.sense$Abundance, by=list(windows.2.sense$Window), FUN=sum)$x
    # make the antisense negative to invert the axes...
    windows.sum$windows.1.antisense <- -(aggregate(windows.1.antisense$Abundance, by=list(windows.1.antisense$Window), FUN=sum)$x)
    windows.sum$windows.2.antisense <- -(aggregate(windows.2.antisense$Abundance, by=list(windows.2.antisense$Window), FUN=sum)$x)
    
    # remove factoring from windows
    windows.sum$Window <- 0:(total.windows - 1) # keep in parentheses for order of operations!!
    # add Inclusive_Start
    windows.sum$Inclusive_Start <- windows.1.sense[1:total.windows, "Inclusive_Start"]
    windows.sum$Inclusive_End   <- windows.1.sense[1:total.windows, "Inclusive_End"]
    
    # Normalize windows
    windows.sum[,c("windows.1.sense","windows.1.antisense")] <- (windows.sum[,c("windows.1.sense","windows.1.antisense")] / normalization.1) * million   
    windows.sum[,c("windows.2.sense","windows.2.antisense")] <- (windows.sum[,c("windows.2.sense","windows.2.antisense")] / normalization.2) * million
    
    # Define min and max positions for the plot - scale up to nearest 250
    max.yaxis <- ceiling(max(windows.sum[,c("windows.1.sense","windows.2.sense")]) / 250) * 250
    min.yaxis <- floor(min(windows.sum[,c("windows.1.antisense","windows.2.antisense")]) / 250) * 250
    
    # Normalize the raw data in preparation for the t-tests
    windows.1.sense$Abundance <- (windows.1.sense$Abundance / normalization.1) * million
    windows.1.antisense$Abundance <- (windows.1.antisense$Abundance / normalization.1) * million
    windows.2.sense$Abundance <- (windows.2.sense$Abundance / normalization.2) * million
    windows.2.antisense$Abundance <- (windows.2.antisense$Abundance / normalization.2) * million
    
    # repeat for sense vs antisense
    # t-test calculation for each window position
    # Null Hypothesis = same mean (mu=0)
    # unpaired, unequal variance, two-sided t-test (Welch Two Sample t-test)
    t.test.results.sense <- data.frame(Window = 0:(total.windows-1),
                                       Inclusive_Start = windows.1.sense[1:total.windows, "Inclusive_Start"],
                                       Inclusive_End   = windows.1.sense[1:total.windows, "Inclusive_End"],
        p.value = sapply(0:(total.windows-1), function (x,m,n) t.test(m[m$Window == x, "Abundance"], n[n$Window == x, "Abundance"], paired=FALSE, var.equal = FALSE, conf.level = 0.95, mu = 0, alternative = "two.sided")$p.value, m=windows.1.sense, n=windows.2.sense))

    # Holm-Bonferoni correction for multiple hypotheses
    #(wikipedia.org/wiki/Holm-Bonferroni_method)
    #    1. Order the p-values (and their respective null hypotheses) from smallest to largest H(1)...H(m); P(1)...P(m)
    #    2. Find the minimal index k where P(k) > alpha / (m + 1 - K)
    #    3. Reject the null hypotheses H(1)...H(k-1); fail to reject the hypotheses H(k)...H(m)

    # sort by p.value
    t.test.correction.sense <- t.test.results.sense[order(t.test.results.sense$p.value),]
    t.test.correction.sense$k <- 1:total.windows

    m <- total.windows
    alpha <- 0.05

    t.test.correction.sense$cutoff <- alpha / (m + 1 - t.test.correction.sense$k)

    t.test.correction.sense$reject.null <- t.test.correction.sense$p.value < t.test.correction.sense$cutoff


    t.test.results.antisense <- data.frame(Window = 0:(total.windows-1),
                                           Inclusive_Start = windows.1.antisense[1:total.windows, "Inclusive_Start"],
                                           Inclusive_End   = windows.1.antisense[1:total.windows, "Inclusive_End"],
        p.value = sapply(0:(total.windows-1), function (x,m,n) t.test(m[m$Window == x, "Abundance"], n[n$Window == x, "Abundance"], paired=FALSE, var.equal = FALSE, conf.level = 0.95, mu = 0, alternative = "two.sided")$p.value, m=windows.1.antisense, n=windows.2.antisense))

    # Holm-Bonferoni correction for multiple hypotheses
    #(wikipedia.org/wiki/Holm-Bonferroni_method)
    #    1. Order the p-values (and their respective null hypotheses) from smallest to largest H(1)...H(m); P(1)...P(m)
    #    2. Find the minimal index k where P(k) > alpha / (m + 1 - K)
    #    3. Reject the null hypotheses H(1)...H(k-1); fail to reject the hypotheses H(k)...H(m)

    # sort by p.value
    t.test.correction.antisense <- t.test.results.antisense[order(t.test.results.antisense$p.value),]
    t.test.correction.antisense$k <- 1:total.windows

    m <- total.windows
    alpha <- 0.05

    t.test.correction.antisense$cutoff <- alpha / (m + 1 - t.test.correction.antisense$k)

    t.test.correction.antisense$reject.null <- t.test.correction.antisense$p.value < t.test.correction.antisense$cutoff

    
    # outputs
    write.table(t.test.correction.sense, 
                file = paste(output.prefix,".sense_pvalues.csv",sep=""),
                quote = FALSE,
                col.names = TRUE,
                row.names = FALSE,
                sep = "\t")
    write.table(t.test.correction.antisense, 
                file = paste(output.prefix,".antisense_pvalues.csv",sep=""),
                quote = FALSE,
                col.names = TRUE,
                row.names = FALSE,
                sep = "\t")            
    
    feature.range = windows.sum[windows.sum$Inclusive_Start <= 0 & windows.sum$Inclusive_End >= 0, c("Inclusive_Start", "Inclusive_End")]
                     
    pdf(paste(output.prefix, ".plot.pdf", sep=""))
    plot(windows.sum$Inclusive_Start, windows.sum$windows.1.sense, type="n", ylim=c(min.yaxis, max.yaxis), xlab=paste("Start position of ",window.size,"bp window relative to ",feature.name," (nt)", sep=""), ylab="Window Coverage in Reads per Million Mapped", xaxs="i", yaxs="i")
    #draw feature setup (only works now with start position...
    abline(h=0, col="darkgrey")
    rect(0, -0.01*(max.yaxis-min.yaxis), max(windows.sum$Inclusive_End), 0.01*(max.yaxis-min.yaxis), col="darkgrey", border="darkgrey")
    # draw region including the start feature
    rect(min(feature.range$Inclusive_Start), min.yaxis, max(feature.range$Inclusive_End), max.yaxis, col="darkgrey", border="darkgrey")
    abline(v=0, col="black")
    
    points(windows.sum$Inclusive_Start, windows.sum$windows.1.sense, type="l", col="blue", lwd=2)
    points(windows.sum$Inclusive_Start, windows.sum$windows.1.antisense, type="l", col="blue", lwd=2)
    points(windows.sum$Inclusive_Start, windows.sum$windows.2.sense, type="l", col="green", lwd=2)
    points(windows.sum$Inclusive_Start, windows.sum$windows.2.antisense, type="l", col="green", lwd=2)

    # add segments for statistically signficant regions
    if (sum(t.test.correction.sense$reject.null, na.rm=T) > 0) {
        segments(x0=(c(t.test.correction.sense[t.test.correction.sense$reject.null,"Window"]) * window.step) - 1000, 
                y0=rep(0.99 * max.yaxis, sum(t.test.correction.sense$reject.null)), 
                x1=(c(t.test.correction.sense[t.test.correction.sense$reject.null,"Window"]) * window.step) - 1000 + window.size, 
                col="red", lwd=3)
    } 
    
    if (sum(t.test.correction.sense$reject.null, na.rm=T) > 0) {
    segments(x0=(c(t.test.correction.antisense[t.test.correction.antisense$reject.null,"Window"]) * window.step) - 1000, 
             y0=rep(0.99 * min.yaxis, sum(t.test.correction.antisense$reject.null)), 
             x1=(c(t.test.correction.antisense[t.test.correction.antisense$reject.null,"Window"]) * window.step) - 1000 + window.size, 
             col="red", lwd=3)
    }
    dev.off()
}


multi.t.test(args[1], #file.1.sense
             args[2], #file.1.antisense
             args[3], #file.2.sense 
             args[4], #file.2.antisense
             as.numeric(args[5]), #normalization.1, 
             as.numeric(args[6]), #normalization.2, 
             args[7], #output.prefix, 
             as.numeric(args[8]), #window.size, 
             as.numeric(args[9]), #window.step,
             args[10]) #feature.name)
