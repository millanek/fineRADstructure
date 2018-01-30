##################################################################
## A simple R example for plotting fineRADstructure output
## Author: Milan Malinsky (millanek@gmail.com), adapted from a Finestructure R Example by Daniel Lawson (dan.lawson@bristol.ac.uk) and using his library of R functions
## Date: 04/04/2016
## Notes:
##    These functions are provided for help working with fineSTRUCTURE output files
## but are not a fully fledged R package for a reason: they are not robust
## and may be expected to work only in some specific cases - often they may require 
## at least minor modifications! USE WITH CAUTION!
## SEE FinestrictureLibrary.R FOR DETAILS OF THE FUNCTIONS
##
## Licence: GPL V3
## 
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.

##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.

##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.


### 1) EDIT THE FOLLOWING THREE LINES TO PROVIDE PATHS TO THE fineRADstructure OUTPUT 
setwd("/PATH/TO/fineRADstructure/OUTPUT/") ## The directory where the files are located
chunkfile<-"INPUT_FILE_chunks.out" ## RADpainter output file
mcmcfile<-"INPUT_FILE_chunks.mcmc.xml" ## finestructure mcmc file
treefile<-"INPUT_FILE_chunks.mcmcTree.xml" ## finestructure tree file
### 2) EDIT THIS PATH TO WHERE YOU WANT THE PLOTS:
plotsFolder <- "/PATH/TO/WHERE/I/WANT/PLOTS/"
### 3) SET VALUES FOR THESE VARIABLES: "analysisName" will be included in output plots
analysisName <- "INPUT_FILE";  maxIndv <- 10000; maxPop<-10000


### 4) EDIT THE PATH TO YOUR COPY of FinestructureLibrary.R
source("/PATH/TO/FinestructureLibrary.R", chdir = TRUE) # read in the R functions, which also calls the needed packages

### 5) EXECUTE THE CODE ABOVE AND THE REST OF THE CODE BELOW
## make some colours
some.colors<-MakeColorYRP() # these are yellow-red-purple
some.colorsEnd<-MakeColorYRP(final=c(0.2,0.2,0.2)) # as above, but with a dark grey final for capped values
###### READ IN THE CHUNKCOUNT FILE
dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=1)) # read in the pairwise coincidence 
###### READ IN THE MCMC FILES
mcmcxml<-xmlTreeParse(mcmcfile) ## read into xml format
mcmcdata<-as.data.frame.myres(mcmcxml) ## convert this into a data frame
###### READ IN THE TREE FILES
treexml<-xmlTreeParse(treefile) ## read the tree as xml format
ttree<-extractTree(treexml) ## extract the tree into ape's phylo format

## Reduce the amount of significant digits printed in the posteror assignment probabilities (numbers shown in the tree):
ttree$node.label[ttree$node.label!=""] <-format(as.numeric(ttree$node.label[ttree$node.label!=""]),digits=2)
 # convert to dendrogram format
tdend<-myapetodend(ttree,factor=1)
## Now we work on the MAP state
mapstate<-extractValue(treexml,"Pop") # map state as a finestructure clustering
mapstatelist<-popAsList(mapstate) # .. and as a list of individuals in populations
popnames<-lapply(mapstatelist,NameSummary) # population names IN A REVERSIBLE FORMAT (I.E LOSSLESS)
## NOTE: if your population labels don't correspond to the format we used (NAME<number>) YOU MAY HAVE TROUBLE HERE. YOU MAY NEED TO RENAME THEM INTO THIS FORM AND DEFINE YOUR POPULATION NAMES IN popnamesplot BELOW
popnamesplot<-lapply(mapstatelist,NameMoreSummary) # a nicer summary of the populations
names(popnames)<-popnamesplot # for nicety only
names(popnamesplot)<-popnamesplot # for nicety only
popdend<-makemydend(tdend,mapstatelist) # use NameSummary to make popdend
popdend<-fixMidpointMembers(popdend) # needed for obscure dendrogram reasons
popdendclear<-makemydend(tdend,mapstatelist,"NameMoreSummary")# use NameMoreSummary to make popdend
popdendclear<-fixMidpointMembers(popdendclear) # needed for obscure dendrogram reasons

	
########################
## Plot 1: COANCESTRY MATRIX
fullorder<-labels(tdend) # the order according to the tree
datamatrix<-dataraw[fullorder,fullorder] # reorder the data matrix

tmpmat<-datamatrix 
tmpmat[tmpmat>maxIndv]<-maxIndv #  # cap the heatmap
pdf(file=paste(plotsFolder,analysisName,"-SimpleCoancestry.pdf",sep=""),height=25,width=25)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=1.1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=1.2))
dev.off()

########################
## Plot 2: POPULATIONS AND COANCESTRY AVERAGES
popmeanmatrix<-getPopMeanMatrix(datamatrix,mapstatelist)

tmpmat<-popmeanmatrix
tmpmat[tmpmat>maxPop]<-maxPop # cap the heatmap
pdf(file=paste(plotsFolder,analysisName,"-PopAveragedCoancestry.pdf",sep=""),height=20,width=20)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=1.1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=1.2))
dev.off()


########################
## Plot 3: POPULATIONS AND COANCESTRY AVERAGES WITH PERHAPS MORE INFORMATIVE LABELS
mappopcorrectorder<-NameExpand(labels(popdend))
mappopsizes<-sapply(mappopcorrectorder,length)
labellocs<-PopCenters(mappopsizes)
xcrt=0
ycrt=45

pdf(file=paste(plotsFolder,analysisName,"-PopAveragedCoancestry2.pdf",sep=""),height=25,width=25)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],labelsx=labels(popdendclear),labelsatx=labellocs,xcrt=xcrt,cols=some.colorsEnd,ycrt=ycrt,dend=tdend,cex.axis=1.1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=1.2),hmmar=c(3,0,0,1))
dev.off()
