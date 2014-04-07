# plot_only.R contains the function (plot.only) which will plot several pairs
# of sense/antisense or ungapped/gapped data on the same plot.
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

load.data <- function (file){
    return (read.table(file, header=TRUE, sep=",", skip=1))
}

normalize.data <- function (data, total.reads, base.unit) {
    return ((data / total.reads) * base.unit)
}

make.window.sum <- function (is.first, data) {
    data$Window <- factor(data$Window)
    total.windows <- length(levels(data$Window))
    
    sum <- aggregate(data$Abundance, by=list(data$Window), FUN=sum)
    names(sum) <- c("Window", "Sum")
    
    if (is.first) {
        sum$Window <- 0:(total.windows - 1)
        sum$Inclusive_Start <- data[1:total.windows, "Inclusive_Start"]
        sum$Inclusive_End <- data[1:total.windows, "Inclusive_End"]
        return (sum[,c("Window","Inclusive_Start","Inclusive_End","Sum")])  
    } else {
        return (sum$Sum)
    }
}

# for testing:
#args <- c("output.prefix",
#          "TSS",
#          0,
#          0,
#          -1000,
#          1000,
#          "sense",
#          "antisense",
#          2,
#          "140403_header_wt.10bpX10bp.sense_allreads.csv", "140403_header_wt.10bpX10bp.antisense_allreads.csv", 12239069, "blue", "WT", "140403_header_rpd1.10bpX10bp.sense_allreads.csv", "140403_header_rpd1.10bpX10bp.antisense_allreads.csv", 11739444, "magenta", "rpd1")


# (file.top, file.bottom, normalization, color, name)
output.prefix <- args[1]
feature.name <- args[2]
window.size <- args[3]
interval.range <- c(as.numeric(args[4]), as.numeric(args[5]))
total.range <- c(as.numeric(args[6]), as.numeric(args[7]))
top <- args[8]
bottom <- args[9]
sets <- as.numeric(args[10])
data <- data.frame(file.1 = sapply(0:(sets-1), function(i,c,m) as.character(m[c + 5*i]), m=args, c=11),
                   file.2 = sapply(0:(sets-1), function(i,c,m) as.character(m[c + 5*i + 1]), m=args, c=11),
                   normalization = sapply(0:(sets-1), function(i,c,m) as.numeric(m[c + 5*i + 2]), m=args, c=11),
                   color = sapply(0:(sets-1), function(i,c,m) as.character(m[c + 5*i + 3]), m=args, c=11),
                   name = sapply(0:(sets-1), function(i,c,m) as.character(m[c + 5*i + 4]), m=args, c=11))
                   
# build windows.sum file
total.names <- c("Windows",
                 "Inclusive_Start",
                 "Inclusive_End", 
                 as.character(sapply(1:sets, function(i,m) paste(m[i, "name"], top, sep="."), m=data)),
                 as.character(sapply(1:sets, function(i,m) paste(m[i, "name"], bottom, sep="."), m=data)))

plotting.data <- cbind(make.window.sum(TRUE, load.data(as.character(data[1, "file.1"]))),
                       sapply(2:sets, function(i,m) make.window.sum(FALSE, load.data(as.character(m[i, "file.1"]))), m=data),
                       sapply(1:sets, function(i,m) make.window.sum(FALSE, load.data(as.character(m[i, "file.2"]))), m=data))

names(plotting.data) <- total.names

# normalize
normalized.data <- cbind(plotting.data[,1:3],
                         sapply(1:sets, function(i,m,n) normalize.data(n[,as.character(paste(m[i,"name"], top, sep="."))], m[i,"normalization"], 1000000), m=data, n=plotting.data),
                         sapply(1:sets, function(i,m,n) normalize.data(n[,as.character(paste(m[i,"name"], bottom, sep="."))], m[i,"normalization"], 1000000), m=data, n=plotting.data))

names(normalized.data) <- total.names

# set parameters for plotting
max.yaxis <- 1.1*max(normalized.data[,4:(4+sets-1)])
min.yaxis <- -1.1*max(normalized.data[,(4+sets):(4+2*sets-1)])

pdf(paste(output.prefix, "plot_only.pdf", sep="."))
plot(0,0, type="n",
     xaxs="i",
     yaxs="i",
     xlim=total.range,
     ylim=c(min.yaxis,max.yaxis),
     xlab=paste("Start position of ",window.size,"bp window relative to ",feature.name," (nt)", sep=""), 
     ylab="Window Coverage in Reads per Million Mapped",
     main=paste(output.prefix, "plot_only.pdf", sep="."))

#draw gene model
gene.width <- c(-0.01*(max.yaxis-min.yaxis), 0.01*(max.yaxis-min.yaxis))
if (feature.name == "TSS" || feature.name == "Start" || feature.name == "start") {
        rect(interval.range[1], gene.width[1], total.range[2], gene.width[2], col="darkgrey", border="darkgrey")
    } else if (feature.name == "end" || feature.name == "End") {
        rect(total.range[1], gene.width[1], interval.range[2], gene.width[2], col="darkgrey", border="darkgrey")
    } else {
        rect(interval.range[1], gene.width[1], interval.range[2], gene.width[2], col="darkgrey", border="darkgrey")
    }
abline(h=0, col="black")

# draw region including the start feature
# may be matrices!
start.range = normalized.data[normalized.data$Inclusive_Start <= interval.range[1] & normalized.data$Inclusive_End >= interval.range[1], c("Inclusive_Start", "Inclusive_End")]
end.range = normalized.data[normalized.data$Inclusive_Start <= interval.range[2] & normalized.data$Inclusive_End >= interval.range[1], c("Inclusive_Start", "Inclusive_End")]

rect(min(start.range$Inclusive_Start), min.yaxis, max(start.range$Inclusive_End), max.yaxis, col="darkgrey", border="darkgrey")
rect(min(end.range$Inclusive_Start), min.yaxis, max(end.range$Inclusive_End), max.yaxis, col="darkgrey", border="darkgrey")
abline(v=interval.range[1], col="black")
abline(v=interval.range[2], col="black")

# plot data
sapply(1:sets, function(i,m,n) points(m$Inclusive_Start, m[,as.character(paste(n[i,"name"], top, sep="."))], type="l", col=as.character(n[i,"color"]), lwd=2), m=normalized.data, n=data)

sapply(1:sets, function(i,m,n) points(m$Inclusive_Start, -m[,as.character(paste(n[i,"name"], bottom, sep="."))], type="l", col=as.character(n[i,"color"]), lwd=2), m=normalized.data, n=data)

# add legend text (may not position correctly, but fixable in Illustrator or Inkscape)
sapply(1:sets, function(i,x,y.start,y.interval,n) text(x,(y.start + i*y.interval), as.character(n[i, "name"]), col=as.character(n[i, "color"]), adj=c(0,0)), x=total.range[1] + 100, y.start=max.yaxis, y.interval=(-0.1*max.yaxis), n=data)
    
dev.off()
