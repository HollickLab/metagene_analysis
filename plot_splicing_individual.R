# plot_splicing_individual.R contains the function (plot.splicing.individual) 
# make bar-plot of abundance of gapped vs ungapped reads
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
options(scipen=9999)

plot.splicing.individual <- function (data, title, gene.length) {
    max.yaxis <- (1.1 * max(c(data$Gapped, data$Ungapped)))
    min.yaxis <- (-0.1 * max(c(data$Gapped, data$Ungapped)))
    plot(0,0, type="n",
         yaxs="i",
         xaxs="i",
         ylim=c(min.yaxis, max.yaxis),
         xlim=c(-1000, gene.length + 1000),
         xlab="Interval Position Relative to Gene Start (nt)",
         ylab="Abundance (Reads per Million Mapped to Interval)",
         main=title)
    rect(0, 0.2*min.yaxis, gene.length-1, -0.2*min.yaxis, col="darkgrey", border="darkgrey")
    abline(h=0, col="black")
    points(data$Inclusive_Start, data$Ungapped, type="h", col="blue")
    points(data$Inclusive_Start, data$Gapped, type="h", col="red")
    text(-400, 0.9*max.yaxis, "Unspliced", col="blue")
    text(-400, 0.8*max.yaxis, "Spliced", col="red")
}

read.data <- function (file) {
    return (read.table(file, header=TRUE, sep=",", skip=1))
}

normalize.data <- function (data, total.reads, base.unit) {
    return ((data / total.reads) * base.unit)
}


gapped.file <- args[1]
ungapped.file <- args[2]
normalization <- as.numeric(args[3])
gene.length <- as.numeric(args[4])
output.prefix <- args[5]
base.unit <- 1000000

data <- read.data(gapped.file)
gene <- as.character(levels(data$Gene))
data <- data[,c("Window","Inclusive_Start","Inclusive_End","Abundance")]
names(data) <- c("Window","Inclusive_Start","Inclusive_End","Gapped")

ungapped <- read.data(ungapped.file)
data$Ungapped <- ungapped$Abundance

data$Gapped <- normalize.data(data$Gapped, normalization, base.unit)
data$Ungapped <- normalize.data(data$Ungapped, normalization, base.unit)

pdf(paste(output.prefix, gene, "plot.pdf", sep="."))
plot.splicing.individual(data, paste(output.prefix, "Reads on Gene:", gene, sep=" "), gene.length)
dev.off()


