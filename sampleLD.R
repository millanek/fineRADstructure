#!/usr/bin/env Rscript
## Author: Daniel Lawson (dan.lawson@bristol.ac.uk)
## Licence: GPL v3

options(warn=-1)
args <- commandArgs(trailingOnly = TRUE)
if(length(args)<2){
    cat("Usage: Rscript sampleLD.R <options> <input chunks file> <output chunks file>
Estimate the LD between all pairs of snps, and reorder them to make high LD snps close.
options:
-s : seed. Set this to get reproducible results.
-n : number of times to sample rad tag variants for correlate. Default: 500
")
    stop("Must provide 2 arguments")
}
tempargs=args

seed=0
nsamples=500
i=1
verbose=FALSE
while(i<length(tempargs)){
    if(tempargs[i] == "-s"){
        seed=as.numeric(tempargs[i+1] )
        tempargs=tempargs[- (i+(0:1)) ]
    }else if(tempargs[i] == "-v"){
        verbose=TRUE
        tempargs=tempargs[- i]
    }else if(tempargs[i]=="-n"){
        nsamples=as.numeric(tempargs[i+1])
        tempargs=tempargs[- (i+(0:1)) ]
    }else{
        i=i+1
    }
}
print(paste(tempargs,collapse=" : "))
if(seed>0) set.seed(seed)
if(length(tempargs)!=2) stop("Must provide input and output files")
infile=tempargs[1]
outfile=tempargs[2]

sampleLD=function(test,nsamples=10,verbose=TRUE){
    tunique=apply(test,1,function(x)unique(x))
    snpcor=test
    snpcor[]=0
    for(sample in 1:nsamples){
        if(verbose)print(paste("Performing sample",sample,"of",nsamples))
        tsample=sapply(tunique,function(x){
            sample(x,1)
        })
        ttest=test
        tmat=sapply(1:dim(test)[1],function(i){
            test[i,]==tsample[i]
        })
        tret=t(tmat)%*%tmat
        snpcor=snpcor + tret
    }
    snpcor/nsamples
}
orderFromSim<-function(x,scale=I){
    x[]=scale(x)
    tdend=as.dendrogram(hclust(dist(x)))
    Rowv <- rowMeans(x)
    tdend=reorder(tdend,Rowv)
    labels(tdend)
}

print(paste("Reading data from",infile))
test=read.table(infile,header=T,sep="\t",as.is=T)
print(paste("Computing LD with",nsamples,"samples"))
testld=sampleLD(test,nsamples)
print(paste("Reordering data"))
testorder=orderFromSim(testld)

print(paste("Writing data to",outfile))
write.table(test[testorder,],file=outfile,quote=FALSE,sep="\t",row.names=FALSE)
