setwd("/Users/tobias/GitHub/seqcap_processor/data/processed/stacey_analysis/run_on_spacemule/")
workdir<-getwd()
workdir
dir()
x<-read.table("species_da_results_1e-3.txt", sep="", header=TRUE)
x
#x$Florisuga <- NULL
y<-names(x)
y
renames <- matrix(c(
"X1166_allele1", "1166_1",
"X1166_allele0", "1166_0",
"X1140_allele1", "1140_1",
"X1140_allele0", "1140_0",
"X1087_allele1", "1087_1",
"X1087_allele0", "1087_0",
"X1086_allele1", "1086_1",
"X1086_allele0", "1086_0",
"X1085_allele1", "1085_1",
"X1085_allele0", "1085_0",
"X1083_allele1", "1083_1",
"X1083_allele0", "1083_0",
"X1082_allele1", "1082_1",
"X1082_allele0", "1082_0",
"X1080_allele1", "1080_1",
"X1080_allele0", "1080_0",
"X1079_allele1", "1079_1",
"X1079_allele0", "1079_0",
"X1074_allele1", "1074_1",
"X1074_allele0", "1074_0",
"X1073_allele1", "1073_1",
"X1073_allele0", "1073_0",
"X1070_allele1", "1070_1",
"X1070_allele0", "1070_0",
"X1068_allele1", "1068_1",
"X1068_allele0", "1068_0",
"X1065_allele1", "1065_1",
"X1065_allele0", "1065_0",
"X1064_allele1", "1064_1",
"X1064_allele0", "1064_0",
"X1063_allele1", "1063_1",
"X1063_allele0", "1063_0",
"X1061_allele1", "1061_1",
"X1061_allele0", "1061_0"),
nrow=34, ncol=2, byrow=TRUE)
renames
# define the columns that should be pursued (remove the 'count' 'fraction' 'similarity' and 'nclusters' column)
mincl.names<-colnames(x)[-(1:4)]
for (i in 1:length(mincl.names)) {
  stopifnot(mincl.names[i] == renames[i,1])
}
mincl.names[1]
renames[1,1]
#make similarity matrix
displaynames <- renames[,2]
nmincls <- length(displaynames)
sim <- matrix(0, ncol=nmincls, nrow=nmincls, dimnames=list(displaynames, displaynames))
for (i in 1:nmincls) {
  for (j in 1:nmincls) {
    coli <- x[,mincl.names[i]]
    colj <- x[,mincl.names[j]]
    w <- coli == colj
    sim[i,j] <- sum(x[w,"fraction"])
  }
}
sim <- pmin(sim,1)
neworder <- c(34,33,26,25,22,21,18,17,20,19,32,31,30,29,8,7,28,27,24,23,10,9,12,11,16,15,4,3,2,1,14,13,6,5)
dividers<-c(0,34)
plot.rectangle <- function(v1,v2,...)
{
polygon(c(v1[1],v2[1],v2[1],v1[1]), c(v1[2],v1[2],v2[2],v2[2]), ...)
}
plot.simmatrix <- function() {
  par(mar= c(0,5,5,0)+.1)
  plot(NULL, xlim=c(0,nmincls), ylim=c(nmincls,0), axes=FALSE, ylab="", xlab="")
  axis(3, at=(1:nmincls)-.5, displaynames[neworder], tick=FALSE, las=2, line=-1)
  axis(2, at=(1:nmincls)-.5, displaynames[neworder], tick=FALSE, las=2, line=-1)
  for (i in 1:nmincls) {
    for (j in 1:nmincls) {
      d <- 1 - sim[neworder[i],neworder[j]]
      plot.rectangle(c(i-1,j-1), c(i,j), col=rgb(d,d,d), border="white")
    }
  }
  for (b in dividers) {
    lines(x=c(-.5,nmincls), y=c(b,b))
    lines(x=c(b,b), y=c(-.5,nmincls))
  }
  #legend(nmincls-4,0,legend = c("0%","25%","50%","75%","100%"),col=c(rgb(1,1,1),rgb(.75,.75,.75),rgb(.5,.5,.5),rgb(.25,.25,.25),rgb(0,0,0)),bg="white", lwd=15, cex=2, box.lty = 1)
}
print(sim[neworder,neworder], digits=2)
#plot.simmatrix()
pdf(file=paste('/Users/tobias/GitHub/seqcap_processor/data/processed/stacey_analysis/run_on_spacemule/', "simmatrix_geonoma_allele_1e3.pdf", sep=""))
plot.simmatrix()
dev.off()

