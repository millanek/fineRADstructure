##################################################################
## Finestructure R Library
## Author: Daniel Lawson (dan.lawson@bristol.ac.uk)
## For more details see www.paintmychromosomes.com ("R Library" page)
## Date: 04/12/2016
## Notes:
##    These functions are provided for help working with fineSTRUCTURE output files
## but are not a fully fledged R package for a reason: they are not robust
## and may be expected to work only in some very specific cases! USE WITH CAUTION!
## SEE FinestrictureExample.R FOR USAGE
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

## THESE ARE REQUIRED LIBRARIES
require(ape)
require(XML)

################################################
## Beginning of functions

########################
## XML reading

simplelabels<-function(x,labmode="mcmc",oneminus=TRUE){
## Utility function for tree extraction
		if(x=="")return("");
		if(labmode=="mcmc"){
			x<-strsplit(x,"-")[[1]][1];
		}else if (labmode=="pop"){
			x<-strsplit(x,"-")[[1]];
			if(length(x)>1){ x<-x[2] }else {return("")}
		}else {return(x);}
		if(oneminus){return(1-as.numeric(x))}
		return(as.numeric(x)) 
}

extractTree<-function(txml,labmode="mcmc",oneminus=TRUE,hidecertain=FALSE){
## extracts the tree from the MAP/Tree file xml
	tmptree<-txml$doc$children$outputFile[[length(txml$doc$children$outputFile)]]
	t<-read.tree(text=xmlValue(tmptree))
	t$node.label<-as.vector(sapply(t$node.label,simplelabels,labmode=labmode,oneminus=oneminus))
	if(hidecertain) t$node.label[which(t$node.label=="1")]<-""
	return(t)
}

extractPop<-function(txml){# ONLY for iterationless outputfiles, such as the MAP/Tree file
## (for extracting states from iterations, use extractValue)
	as.character(xmlChildren(txml$doc$children$Pop)$text)[6]
}
getV<-function(it,v){
## extract element v from iteration it
	return(xmlValue(xmlChildren(it)[[v]]))
}

extractValue<-function(txml,v,getNames=FALSE){
## Important utility function for extracting element v from the xml
	tmpits<-txml$doc$children$outputFile[which(names(txml$doc$children$outputFile)=="Iteration")]
	if(getNames) num<-sapply(tmpits,getV,v="Number")
	res<-sapply(tmpits,getV,v=v)
	if(!getNames) names(res)<-NULL
	else names(res)<-num
	return(res)
}

as.data.frame.myres<-function(txml){
## Converts our xml file format into a matrix, one row per iteration
	tmpits<-txml$doc$children$outputFile[which(names(txml$doc$children$outputFile)=="Iteration")]
#	xmlChildren(tmpits[[1]])
	cnames<-names(tmpits[[1]])
	res<-as.data.frame(matrix(nrow=length(tmpits),ncol=length(cnames)))
	for(i in 1:length(cnames)) {
		res[,i]<-extractValue(txml,cnames[i])
	}
	colnames(res)<-cnames
	res[,which(!names(res) %in% c("Pop","P","Q"))]  <-apply(res[,which(!names(res) %in% c("Pop","P","Q"))],2,as.numeric)
	res
}

popAsList<-function(s) {
## convert bracketed population string into population list
	s<-gsub(")","",s)
	s<-strsplit(s,c("(",")"),fixed=TRUE)[[1]]
	s<-s[sapply(s,nchar)>0]
	sapply(s,strsplit,",",fixed=TRUE)
}

popCor<-function(pop1,pop2){
#	((sum(sapply(pop1,function(x){x %in% pop2})) + sum(sapply(pop2,function(x){x %in% pop1})))/(length(pop1)+length(pop2)))^2
	m<-sum(sapply(pop1,function(x){x %in% pop2}))
	(m/(length(pop1)+length(pop2)-m))^2
}

popCorMat<-function(state1,state2){
## correlation matrix of the populations between two states
	r<-sapply(state1,function(x){sapply(state2,popCor,x)})	
	dimnames(r)<-NULL
	r
}

popCorMatDiag<-function(state1,state2){
    ## diagonalise the population correlation matrix by finding the best population for each
                                        #rmat<-popCorMat(state2,state1) # do it this way around to allow state2 to be a "reference state"
    rmat<-popCorMat(state1,state2)
	rmat2<-(rmat)
	if(length(state1)==1) {return(rmat)
	}else bestmatch<-apply(rmat,1,function(x){y<-which(x==max(x)); if(length(y)>1) y<-sample(y,1); y})
	for(i in 1:length(bestmatch)) rmat2[,i]<-rmat[,bestmatch[i]]
	rmat2
}

stateCor<-function(state1,state2,fn="popCorMat"){
## Correlation between states
	rmat<-get(fn)(state2,state1) # do it this way around to allow state2 to be a "reference state"
	if(length(state1)==1) return(rmat[1])
	else bestmatch<-apply(rmat,1,function(x){y<-which(x==max(x)); if(length(y)>1) y<-sample(y,1); y})
	corvec<-rep(0,length(bestmatch))
	for(i in 1:length(bestmatch)) corvec[i]<-rmat[i,bestmatch[i]]
	lens<-as.vector(sapply(state1,length))
	sum(corvec*lens/sum(lens))
}

groupingAsMatrix<-function(pg,res=NULL)
{
        tmp<-unlist(pg)
        names(tmp)<-NULL
        if(is.null(res)){
                res<-matrix(0,ncol=length(tmp),nrow=length(tmp))
                colnames(res)<-rownames(res)<-tmp
        }else res[]<-0
        for(i in 1:length(pg)) {
                tdat<-expand.grid(pg[[i]],pg[[i]], KEEP.OUT.ATTRS = FALSE)
                for(j in 1:dim(tdat)[1]) res[as.character(tdat[j,1]),as.character(tdat[j,2])]<-res[as.character(tdat[j,1]),as.character(tdat[j,2])]+1
        }
        res
}


matrixAsPopList<-function(mat)
{
	if(any((mat!=0) * (mat!=1)>0)) stop("Must be binary matrix")
	namlist<-dimnames(mat)[[1]]
	res<-list(namlist[length(namlist)])
	namlist<-namlist[-length(namlist)]
	while(1){
		if(length(namlist)<1) break;
		i<-tail(namlist,1)
		found<-FALSE
		for(j in 1:length(res)) {
			if(any(mat[i,res[[j]]]==1)) {
			  res[[j]]<-c(res[[j]],i)
			  found<-TRUE
			  break
			}
		}
		if(!found) res<-c(res,i)
		namlist<-namlist[-length(namlist)]
	}
	res
}

aveMat<-function(pophist) {
## Compute pairwise coincidence from the population history vector
	phl<-lapply(pophist,popAsList)
	for(i in 1:length(pophist)){cat("-")};cat("\n");
	res<-groupingMatrix(pophist[1])
	cat(".")
	for(i in 2:length(pophist)){
		tp<-groupingMatrix(pophist[i])
		res<-res+tp
		cat(".")
	}
	cat("\n")
	res<-res/length(pophist)
	res
}

matColMeans<-function(mat,poplist)
## Sums over columns in a matrix, by grouping all columns listed in poplist
## e.g. if mat is M*N matrix and poplist is length K, returns a M*K matrix
## the names of poplist are used to assign names to the returned matrix
{
	res<-matrix(0,nrow=dim(mat)[1],ncol=length(poplist))
	colindex<-lapply(poplist,function(x){which(colnames(mat)%in%x)})
	res<-t(apply(mat,1,function(x){
		sapply(colindex,function(y){mean(x[y])})
	}))
	colnames(res)<-names(poplist)
	rownames(res)<-rownames(mat)
	res
}

matColSums<-function(mat,poplist)
## Sums over columns in a matrix, by grouping all columns listed in poplist
## e.g. if mat is M*N matrix and poplist is length K, returns a M*K matrix
## the names of poplist are used to assign names to the returned matrix
{
	res<-matrix(0,nrow=dim(mat)[1],ncol=length(poplist))
	colindex<-lapply(poplist,function(x){which(colnames(mat)%in%x)})
	res<-t(apply(mat,1,function(x){
		sapply(colindex,function(y){sum(x[y])})
	}))
        if(dim(res)[1]==1) res<-t(res)
	colnames(res)<-names(poplist)
	rownames(res)<-rownames(mat)
	res
}


avgDist<-function(avgmat,poplist)
{
	res<-vector("numeric",length(poplist))
	for(i in 1:length(poplist)){cat("-")};cat("\n");
	for(i in 1:length(poplist)){
		tp<-groupingMatrix(poplist[i])
		res[i]<-sum(abs(tp-avgmat))
		cat(".")
	}
	cat("\n")
	res
}

gelman.diag<-function (x, confidence = 0.95, transform = FALSE, autoburnin = TRUE)
{# MODIFY THIS (from coda) TO WORK WITH MY DATA?
    x <- as.mcmc.list(x)
    if (nchain(x) < 2)
        stop("You need at least two chains")
    if (autoburnin && start(x) < end(x)/2)
        x <- window(x, start = end(x)/2 + 1)
    Niter <- niter(x)
    Nchain <- nchain(x)
    Nvar <- nvar(x)
    xnames <- varnames(x)
    if (transform)
        x <- gelman.transform(x)
    x <- lapply(x, as.matrix)
    S2 <- array(sapply(x, var, simplify = TRUE), dim = c(Nvar,
        Nvar, Nchain))
    W <- apply(S2, c(1, 2), mean)
    xbar <- matrix(sapply(x, apply, 2, mean, simplify = TRUE),
        nrow = Nvar, ncol = Nchain)
    B <- Niter * var(t(xbar))
    if (Nvar > 1) {
        if (is.R()) {
            CW <- chol(W)
            emax <- eigen(backsolve(CW, t(backsolve(CW, B, transpose = TRUE)),
                transpose = TRUE), symmetric = TRUE, only.values = TRUE)$values[1]
        }
        else {
            emax <- eigen(qr.solve(W, B), symmetric = FALSE,
                only.values = TRUE)$values
        }
        mpsrf <- sqrt((1 - 1/Niter) + (1 + 1/Nvar) * emax/Niter)
    }
    else mpsrf <- NULL
    w <- diag(W)
    b <- diag(B)
    s2 <- matrix(apply(S2, 3, diag), nrow = Nvar, ncol = Nchain)
    muhat <- apply(xbar, 1, mean)
    var.w <- apply(s2, 1, var)/Nchain
    var.b <- (2 * b^2)/(Nchain - 1)
    cov.wb <- (Niter/Nchain) * diag(var(t(s2), t(xbar^2)) - 2 *
        muhat * var(t(s2), t(xbar)))
    V <- (Niter - 1) * w/Niter + (1 + 1/Nchain) * b/Niter
    var.V <- ((Niter - 1)^2 * var.w + (1 + 1/Nchain)^2 * var.b +
        2 * (Niter - 1) * (1 + 1/Nchain) * cov.wb)/Niter^2
    df.V <- (2 * V^2)/var.V
    df.adj <- (df.V + 3)/(df.V + 1)
    B.df <- Nchain - 1
    W.df <- (2 * w^2)/var.w
    R2.fixed <- (Niter - 1)/Niter
    R2.random <- (1 + 1/Nchain) * (1/Niter) * (b/w)
    R2.estimate <- R2.fixed + R2.random
    R2.upper <- R2.fixed + qf((1 + confidence)/2, B.df, W.df) *
        R2.random
    psrf <- cbind(sqrt(df.adj * R2.estimate), sqrt(df.adj * R2.upper))
    dimnames(psrf) <- list(xnames, c("Point est.", paste(50 *
        (1 + confidence), "% quantile", sep = "")))
    out <- list(psrf = psrf, mpsrf = mpsrf)
    class(out) <- "gelman.diag"
    out
}


matrixAsPopList<-function(mat)
{
	if(any((mat!=0) * (mat!=1)>0)) stop("Must be binary matrix")
	namlist<-dimnames(mat)[[1]]
	res<-list(namlist[length(namlist)])
	namlist<-namlist[-length(namlist)]
	while(1){
		if(length(namlist)<1) break;
		i<-tail(namlist,1)
		found<-FALSE
		for(j in 1:length(res)) {
			if(any(mat[i,res[[j]]]==1)) {
			  res[[j]]<-c(res[[j]],i)
			  found<-TRUE
			  break
			}
		}
		if(!found) res<-c(res,i)
		namlist<-namlist[-length(namlist)]
	}
	res
}

popListAsPop<-function(pl)
{
	tf<-function(x){paste("(",paste(x,collapse=","),")",sep="")}
	pl2<-sapply(pl,tf)
	paste(pl2,collapse="")
}

writePopFile<-function(pop,file)
{
	cat(paste("<Pop>",pop,"</Pop>\n",sep=""),file=file)
	invisible(NULL)
}

pairwiseCoincidence<-function(pl)
{
	pl2<-popAsList(pl)
	pl2<-sapply(pl2,function(x){gsub("[0-9]","",x)})
	pcft<-function(x){sapply(x,function(y){sum(x==y)/length(x)})}
	tmp<-unlist(sapply(pl2,pcft))
	mean(tmp)
}## Penalises merges of different labels

pairwiseEfficiency<-function(pl){pairwiseCoincidence(pl)}

pairwiseSplit<-function(pl)
{
	pl2<-popAsList(pl)
	pl2<-sapply(pl2,function(x){paste("a",x,sep="")})# cope with numbers only
	pl2<-sapply(pl2,function(x){gsub("[0-9]","",x)})
	plistall<-unlist(pl2)
	pcounts<-sapply(unique(plistall),function(x){sum(x==plistall)})
	pcft<-function(x){sapply(x,function(y){sum(y==x)/pcounts[y]})}
	tmp<-unlist(sapply(pl2,pcft),use.names = FALSE)
	mean(tmp)
}## Penalises splitting of labels, allows merges of different labels


pairwiseAccuracy<-function(pl)
{
	pl2<-popAsList(pl)
	pl2<-sapply(pl2,function(x){paste("a",x,sep="")})# cope with numbers only
	pl2<-sapply(pl2,function(x){gsub("[0-9]","",x)})
	plistall<-unlist(pl2)
	pcounts<-sapply(unique(plistall),function(x){sum(x==plistall)})
	pcft<-function(x){sapply(x,function(y){sum(y!=x)/length(x)*(1-sum(x==y)/pcounts[y])})}
	tmp<-unlist(sapply(pl2,pcft),use.names = FALSE)
	1-mean(tmp)
}## For each non-pop indiv in your group, get a penalty of the proportion of your pop NOT in this group
## i.e. penalise only when correct labelling as a population becomes impossible
# i.e. no penalty for mixing two populations, no penalty for splitting a population
# Get a penalty when you share a group with another pop and are also present in another group


getPop<-function(str,half=TRUE,n=100)
{
txml<-xmlTreeParse(str)
tres<-as.data.frame.myres(txml)
if(half) tres<-tres[(dim(tres)[1]/2):dim(tres)[1],]
if(n>0) mys<-trunc(seq(floor(length(tres$Pop)/n),length(tres$Pop),length.out=n))
else mys<-1:length(tres$Pop)
tres$Pop[mys]
}

testMCMC<-function(poplist)
{
	r1<-sapply(poplist,pairwiseCoincidence)
	r2<-sapply(poplist,pairwiseAccuracy)
	
	ret<-c(mean(r1),mean(r2),quantile(r1,c(0.025,0.975)),quantile(r2,c(0.025,0.975)))
	names(ret)<-c("meanPc","meanPa","2.5%Pc","97.5%Pc","2.5%Pa","97.5%Pa")
	ret
}

###########################
## Obtaining the population-averaged chunkcount matrix

getPopCountMatrix<-function(datamatrix,mapstatelist,remdiag=TRUE){
  ## Takes the chunkcount matrix and the list of which individuals are in which population
  ## Returns a matrix of the number of individuals in each population pair
  ## ASSUMES THAT WE EXCLUDE THE DIAGONAL, AND NOTHING ELSE, FROM THE COMPARISONS
  popnmatrix<-datamatrix
  popnmatrix[]<-1
  if(remdiag) diag(popnmatrix)<-0
  popcountmatrix<-matColSums(popnmatrix,mapstatelist)
  popcountmatrix<-matColSums(t(popcountmatrix),mapstatelist)
  popcountmatrix
}

getPopMatrix<-function(datamatrix,mapstatelist,correction=1,remdiag=TRUE){
  ## Takes the chunkcount matrix and the list of which individuals are in which population
  ## Returns a matrix of the chunkcount average for each population pair
  datamatrix[is.nan(datamatrix[])]<-0
  popcountmatrix<-getPopCountMatrix(datamatrix,mapstatelist,remdiag=remdiag)
  popmatrix<-matColSums(datamatrix,mapstatelist)
  popmatrix<-t(matColSums(t(popmatrix),mapstatelist))
  popmatrix<-popmatrix/popcountmatrix
  popmatrix[is.nan(popmatrix[])]<-0
  (popmatrix/correction)
}

getPopIndices<-function(indnames,mapstatelist){
  # takes a list of individual names, for example as rownames(datamatrix), and the list of individuals in each population
  # returns a list of which population each individual is in
  popindices<-sapply(indnames,function(x){
    which(sapply(mapstatelist,function(y){x %in%y}))
  })
  popindices
}

getOrderedPopMatrix<-function(datamatrix,mapstatelist,tlabels=NULL){
  # takes the chunkcount matrix and the list of individuals in each population
  # returns a matrix of the dimension of each population pair
  if(is.null(tlabels)) tlabels<-names(mapstatelist)
  popmeanmatrix<-datamatrix
  popindices<-getPopIndices(rownames(datamatrix),mapstatelist)
  popmatrix<-getPopMatrix(datamatrix,mapstatelist)
  popmatrix<-popmatrix[unique(popindices),unique(popindices)]
  rownames(popmatrix)<-colnames(popmatrix)<-tlabels
  t(popmatrix)
}

getPopMeanMatrix<-function(datamatrix,mapstatelist){
  # takes the chunkcount matrix and the list of individuals in each population
  # returns a matrix of the same dimension, where x[i,j] has been replaced by the average x[i,j] for each population pair
  popmeanmatrix<-datamatrix
  popindices<-getPopIndices(rownames(datamatrix),mapstatelist)
  popmatrix<-getPopMatrix(datamatrix,mapstatelist)
  for(i in 1:dim(popmeanmatrix)[1]) for(j in 1:dim(popmeanmatrix)[2]){
    popmeanmatrix[i,j]<-popmatrix[popindices[i],popindices[j]]
  }
  popmeanmatrix
}


###########################
## Start of dendrogram functions
fixMidpointMembers<-function(x)
{
    ## Fixes the midpoint of a dendrogram node (dendrapply and removing midpoints on tips needed)
attr(x, "x.member")<-NULL
attr(x, "value")<-NULL
if(is.leaf(x)) {
	attr(x, "label")<-attr(x, "label")[[1]]
	attr(x, "midpoint")<-NULL
	attr(x, "members")<-1
}
if(!is.leaf(x)) {
    attr(x, "members")<-length(labels(x))
}
x
}

cutdend<- function(tdend,height,summary="NameSummary") {
    ## Cuts the dendrogram at a given height and returns the corresponding dendrogram.  This is invertible with summary="NameSummary" or otherwise simplified in summary="NameMoreSummary"
    testj<-cut(tdend,height)
    newdend<-testj$upper
    newdend<-dendrapply(newdend,function(x){
        if(is.leaf(x)) {
            attr(x,"height")<-0
            tbranch<-as.numeric(gsub("Branch","",attr(x,"label")))
            attr(x,"label")<-get(summary)(labels(testj$lower[[tbranch]]))
        }
        x
    })
    ## Fixing the mid points
    newdend <- dendrapply(newdend, fixMidpointMembers)
##        tdend <- fix_members_attr.dendrogram(tdend)   
    newdend <- suppressWarnings(midcache.dendrogram(newdend)) # fixing the middle point thing
    newdend
}

relabelDend<-function(tdend,namemap=NULL){
    ## Relabel a dendrogram (which must have unique labels)
    ## Either (default) generate the relabelling from NameMoreSummary,
    ## Or use a map provide (in named list format)
    if(is.null(namemap)) {
        namemap<-sapply(NameExpand(labels(tdend)),NameMoreSummary)
    }
    dendrapply(tdend,function(x){
        if(!is.null(attr(x,"label")) && attr(x,"label")%in%names(namemap)){
            attr(x,"label")<-namemap[[attr(x,"label")]]
        }
        x
    })
}


makemydend<- function(tdend,lablist,summary="NameSummary") {
	for(i in 1:2) {
	j<-i+1;	if(j==3)j<-1
	test<-which(sapply(lablist,function(x){all(labels(tdend[[i]]) %in% x)}))
	if(length(test)>0) {
		if(all(lablist[[test]]%in%labels(tdend[[i]]))){
		testj<-cut(tdend,attributes(tdend)$height)$upper
		attr(testj[[i]],"height")<-0
                attr(testj[[i]],"label")<-get(summary)(lablist[[test]])
		testj[[j]]<-tdend[[j]]
		tdend<-testj
		}
	}else if(length(labels(tdend[[i]])>1)) tdend[[i]]<-makemydend(tdend[[i]],lablist,summary)
    }
        ## Fixing the mid points
        tdend <- dendrapply(tdend, fixMidpointMembers)
##        tdend <- fix_members_attr.dendrogram(tdend)   
        tdend <- suppressWarnings(midcache.dendrogram(tdend)) # fixing the middle point thing
        tdend

}

popIn<-function(x){
	tlab2<-gsub("[0-9]","",x)
	unique(tlab2)
}

popIn2<-function(x){
	tlab2<-gsub("[0-9,-]","",x)
	unique(tlab2)
}

popUnder<-function(x){
	tlab2<-gsub("[0-9]","",labels(x))
	unique(tlab2)
}

getTip1<-function(x){
	while(!is.leaf(x)) x<-getTip1(x[[1]])
	x
}

relabel<-function(x,l1,l2){
	if(!is.null(attributes(x)$label)) attributes(x)$label<-l2[which(l1==attributes(x)$label)]
	x
}

dend.relabel<-function(x,l1,l2)
{
	dendrapply(x,relabel,l1=l1,l2=l2)
}

## Reversing the above (as much as possible)
## is only possible for dendrograms with *binary* splits


.memberDend <- function(x) {
    r <- attr(x,"x.member")
    if(is.null(r)) {
	r <- attr(x,"members")
	if(is.null(r)) r <- 1L
    }
    r
}

.midDend <- function(x)
    if(is.null(mp <- attr(x, "midpoint"))) 0 else mp


midcache.dendrogram <- function (x, type = "hclust", quiet=FALSE)
{
    ## Recompute "midpoint" attributes of a dendrogram, e.g. after reorder().

    type <- match.arg(type) ## currently only "hclust"
    stopifnot( inherits(x, "dendrogram") )
    setmid <- function(d, type) {
	if(is.leaf(d))# no "midpoint"
	    return(d)
	k <- length(d)
	if(k < 1)
	    stop("dendrogram node with non-positive #{branches}")
	r <- d # incl. attributes!
	midS <- 0
	for(j in 1L:k) {
	    r[[j]] <- unclass(setmid(d[[j]], type))
	    midS <- midS + .midDend(r[[j]])
	}
	if(!quiet && type == "hclust" && k != 2)
	    warning("midcache() of non-binary dendrograms only partly implemented")
	## compatible to as.dendrogram.hclust() {MM: doubtful if k > 2}
	attr(r, "midpoint") <- (.memberDend(d[[1L]]) + midS) / 2
	r
    }
    setmid(x, type=type)
}
    
as.hclust.dendrogram.orig <- function(x, ...)
{
    stopifnot(is.list(x), length(x) == 2)
    n <- length(ord <- unlist(x))
    stopifnot(n == attr(x, "members"))
    n.h <- n - 1L
    ## labels: not sure, if we'll use this; there should be a faster way!
    labsu <- unlist(labels(x))
    labs <- labsu[sort.list(ord)]
    x <- .add.dendrInd(x)

    SIMP <- function(d) {
	if(is.leaf(d)) {
	    - as.vector(d)# dropping attributes
	} else {
	    j <<- j + 1L
	    height[j] <<- attr(d, "height")
	    inds[[j]] <<- attr(d, ".indx.")
	    attributes(d) <- NULL # drop all, incl. class
	    ## recursively apply to components:
	    d[] <- lapply(d, SIMP)
	    d
	}
    }

    height <- numeric(n.h);  inds <- vector("list",n.h);  j <- 0L
    xS <- SIMP(x)
    ii <- sort.list(height)

    merge <- matrix(NA_integer_, 2L, n.h)
    for(k in seq_len(n.h)) {
	if(k < n.h) { in.k <- inds[[ ii[k] ]] ; s <- xS[[in.k]] } else s <- xS
	##cat(sprintf("ii[k=%2d]=%2d -> s=xS[[in.k]]=", k, ii[k])); str(s)
	s<-lapply(s, as.integer)
	stopifnot(length(s) == 2L, all( vapply(s, is.integer, NA) ))# checking..
	merge[,k] <- unlist(s)
	if(k < n.h)
	    xS[[in.k]] <- + k
    }

    r <- list(merge = t(merge),
	      height = height[ii],
	      order = ord,
	      labels = labs,
	      call = match.call(),
	      method = NA_character_,
	      dist.method = NA_character_)
    class(r) <- "hclust"
    r
}

##' add the c(i1,i2,..) list indices to each non-leaf of a dendrogram
##' --> allowing "random access" into the dendrogram
.add.dendrInd <- function(x)
{
    add.I <- function(x, ind) {
        if(!is.leaf(x)) {
            for(i in seq_along(x))
                x[[i]] <- add.I(x[[i]], c(ind, i))
            attr(x, ".indx.") <- ind
        }
        x
    }
    ## apply recursively:
    add.I(x, integer(0))
}


plot.node.labels<-function(x,m=0)
{
	if(!is.leaf(x)) {
		text(attr(x,"height")+0.5,m+attr(x,"midpoint"),attr(x,"edgetext"))
		plot.node.labels(x[[1]],m+attr(x,"midpoint"))
		plot.node.labels(x[[2]],m)
	}
}

simpledend<-function(x) {
	print(x)
	if(length(popUnder(x[[1]]))>1){
#		print(paste("PopUnder1:",popUnder(x[[1]])))
#		print(x[[1]])
		r1<-simpledend(x[[1]])
		attributes(r1)$members<-length(labels(r1))
	}else {
		r1<-getTip1(x[[1]])
		attributes(r1)$label<- paste(popUnder(x[[1]]),length(labels(x[[1]])),sep="")
#		attributes(r1)$height<-attributes(x[[1]])$height
	}
	if(length(popUnder(x[[2]]))>1){
#		print(paste("PopUnder2:",popUnder(x[[2]])))
#		print(x[[2]])
		r2<-simpledend(x[[2]])
		attributes(r2)$members<-length(labels(r2))
	}else {
		r2<-getTip1(x[[2]])
		attributes(r2)$label<- paste(popUnder(x[[2]]),length(labels(x[[2]])),sep="")
#		attributes(r2)$height<-attributes(x[[2]])$height
	}
	x[[1]]<-r1
	x[[2]]<-r2
	attributes(x)$members <- length(labels(x))
	x<-midcache.dendrogram(x)
	return(x)
}


my.as.hclust.phylo<-function (x,tol=0.01, ...)
{
    if (!is.ultrametric(x,tol))
        stop("the tree is not ultrametric")
    if (!is.binary.tree(x))
        stop("the tree is not binary")
    n <- length(x$tip.label)
    bt <- rev(branching.times(x))
    N <- length(bt)
    nm <- as.numeric(names(bt))
    merge <- matrix(NA, N, 2)
    for (i in 1:N) {
        ind <- which(x$edge[, 1] == nm[i])
        for (k in 1:2) merge[i, k] <- if (x$edge[ind[k], 2] <=
            n)
            -x$edge[ind[k], 2]
        else which(nm == x$edge[ind[k], 2])
    }
    names(bt) <- NULL
    obj <- list(merge = merge, height = bt, order = 1:(N + 1),
        labels = x$tip.label, call = match.call(), method = "unknown")
    class(obj) <- "hclust"
    obj
}# Adds tolerance to the function

positiveheights<-function(x){
	if(attributes(x)$height<0){attributes(x)$height<-0.0001;}
x}

flattenheights<-function(x,factor=0.25){
	if(attributes(x)$height>0){attributes(x)$height<-attributes(x)$height^factor;}
x}

myapetodend<-function(ttree,simplify=FALSE,tol=0.1,factor=0.25){
	nodelab<-ttree$node.label
	ttree$node.label<-NULL
	htree<-my.as.hclust.phylo(ttree,tol)
	dend<-as.dendrogram(htree)
#	if(!is.null(nodelab))dend<-setNodeLabels(dend,nodelab)
	dend<-dendrapply(dend,positiveheights)
	if(simplify) dend<-simpledend(dend)
	dend<-dendrapply(dend,flattenheights,factor=factor)
	if(!is.null(nodelab))dend<-setNodeLabels(dend,nodelab)
	dend
}

findsplits<-function(tdend,tol=0.75){
  i<-0
  tlist<-list()
  popExt <<- function(n,tol=0.75) {
      if(!is.leaf(n)) {
	myval<-as.numeric(attr(n, "edgetext"))
        if(length(myval)>0) {if(myval>tol) {
		i<<-i+1
		tlist[[i]]<<-labels(n)
	}}
      }
      n
  }
  tmp<-dendrapply(tdend, popExt,tol=tol)
  tlist
}

findsplitsNoText<-function(tdend){
  i<-0
  tlist<-list()
  popExt <<- function(n) {
      if(!is.leaf(n)) {
		i<<-i+1
		tlist[[i]]<<-labels(n)
      }
      n
  }
  tmp<-dendrapply(tdend, popExt)
  tlist
}


checkPops<-function(l1,l2) {
	res<-rep(FALSE,length(l1))
	for( pon in 1:length(l1)) {
		p<-l1[[pon]]
		plist<-lapply(p,function(x){grep(x,l2)})
		tvec<-rep(FALSE,length(plist[[1]]))
		if(length(plist[[1]])>0){
		for(pop in 1:length(plist[[1]])) tvec[pop]<-all(sapply(plist,function(x){plist[[1]][pop] %in% x}))
		}
		for(pon2 in which(tvec)) {
			p2<-l2[[plist[[1]][pon2]]]
			if(all(sapply(p2,function(x){x %in% p}))) res[[pon]]<-TRUE 
		}
#		res[[pon]]<-any(tvec)
	}
	return(sum(res)/length(res))
}

###################################
## DENDROGRAM MANIPULATION - probably worth ignoring this

setCols<-function(n,lab.cex=0.5,cex=0,cols=1:5,collist){
	if(is.leaf(n)){
            a <- attributes(n)
	    label<-gsub("[0-9]","",labels(n))
	    area<-which(sapply(lapply(collist,function(x,pattern){which(x==pattern)},pattern=label),length)>0)
	    thiscol=cols[area]
             attr(n, "nodePar") <-
                 c(a$nodePar, list(cex=cex,lab.cex=lab.cex,lab.col = thiscol))
	}
	n
}

setCols2<-function(n,lab.cex=0.5,cex=0,cols=1:5,collist){
	if(is.leaf(n)){
            a <- attributes(n)
	    labeltmp<-gsub("[0-9,-]","",labels(n))
	    label<-strsplit(labeltmp,";")[[1]]
	    area<-which(sapply(lapply(collist,function(x,pattern){which(x%in%pattern)},pattern=label),length)>0)
	print(paste(labels(n),":::",area))
	    if(length(area)==1) {thiscol=cols[area]}else {thiscol="yellow"}
             attr(n, "nodePar") <-
                 c(a$nodePar, list(cex=cex,lab.cex=lab.cex,lab.col = thiscol))
	}
	n
}


getCols<-function(n){
	if(is.leaf(n)) return(attr(n,"nodePar")$lab.col)
}

setLabels<-function(din,mylabs){
 i <- 0
  edgeLab <<- function(n,lablist) {
      if(is.leaf(n)) {
        a <- attributes(n)
        i <<- i+1
        attr(n, "label") <-lablist[i]
      }
      n
  }
  dendrapply(din, edgeLab,mylabs)
}	

setDend<-function(tdend,labrep,tdendc){
  replacetip <<- function(n) {
      if(is.leaf(n)) {
        if(attr(n, "label")==labrep) n<-tdendc
      }
      n
  }
  dendrapply(tdend,replacetip)
}

expandDendElements<-function(tdend,labrep,newsize)
{
tdend<-dendrapply(tdend,function(x){
if(!is.leaf(x)){attr(x,"members")<-attr(x,"members")+2*sum(labels(x)%in% labrep);}else{ if(attr(x,"label")%in%labrep) attr(x,"members")<-newsize}; x
})
tdend
}

removeDend<-function(tdendc,rem)
{
  while(any(labels(tdendc)%in%rem)){
  if(all(labels(tdendc[[1]])%in%rem)) {tdendc<-tdendc[[2]]
  }else if(all(labels(tdendc[[2]])%in%rem)){ tdendc<-tdendc[[1]]}
  }
  tdendc
}

setNodeLabels<-function(din,mylabs){
 i <- 0
  edgeLab <<- function(n,lablist) {
      if(!is.leaf(n)) {
        a <- attributes(n)
        i <<- i+1
#        attr(n, "label") <-lablist[i]
	if(nchar(lablist[i])>0 &&i>1){
        #attr(n, "edgePar") <- list(p.col="plum")
        attr(n, "edgetext") <-lablist[i]# paste(attr(n,"members"),"members")
	}
      }
      n
  }
  dendrapply(din, edgeLab,mylabs)
}	

setEdgePar<-function(din,mylabs){
 i <- 0
  setAtt<<-function(n,lab){
        a <- attributes(n)
#        attr(n, "label") <-lablist[i]
	if(nchar(lab)>0 &&i>1){
        #attr(n, "edgePar") <- list(p.col="plum")
        attr(n, "edgePar") <-lab# paste(attr(n,"members"),"members")
	}
	n
  }
  edgeLab <<- function(n,lablist) {
      if(!is.leaf(n)) {
        i <<- i+1
	n<-setAtt(n,lablist[[i]])
	#n[[1]]<-setAtt(n[[1]],lablist[[i]])
	#n[[2]]<-setAtt(n[[2]],lablist[[i]])
     }
      n
  }
  dendrapply(din, edgeLab,mylabs)
}

fixEdges<-function(din){
	res<-dendrapply(din, function(x){
		if(is.leaf(x)) {attributes(x)<-list(members=attr(x,"members"),height=attr(x,"height"),leaf=attr(x,"leaf"),label=attr(x,"label"),class=attr(x,"class"));x}
	})
	class(res)<-"dendrogram"
	res
}

applyTreeTextPart<-function(x1,x2){
	if(labels(x1[[1]])==labels(x2[[1]])){
	  attributes(x1[[1]])$edgetext<-attributes(x2[[1]])$edgetext
	  attributes(x1[[2]])$edgetext<-attributes(x2[[2]])$edgetext
	}else if (labels(x1[[1]])==labels(x2[[2]])){
	  attributes(x1[[1]])$edgetext<-attributes(x2[[2]])$edgetext
	  attributes(x1[[2]])$edgetext<-attributes(x2[[1]])$edgetext
	}else{
	  stop("Trees don't match!")
	}
	x1
}

dend.relabel<-function(x,l1,l2)
{
	dendrapply(x,relabel,l1=l1,l2=l2)
}

edgelablist<-function(tdend){
  i<-0
  tlist<-list()
  ell<-vector("character",0)
  popExt <<- function(n) {
      if(!is.leaf(n)) {
	myval<-as.numeric(attr(n, "edgetext"))
        if(length(myval)>0) {
		i<<-i+1
		tlist[[i]]<<-labels(n)
		ell<<-c(ell,myval)
	}
      }
      n
  }
  tmp<-dendrapply(tdend, popExt)
  list(edgelab=ell,inds=tlist)
}

matchlabels<-function(lab,tlablist){
  r1<-sapply(tlablist,function(x){all(x %in% lab)})
  if(length(r1)==0) return(r1)
  tmp<-which(r1)
  r2<-sapply(tmp,function(x){all(lab %in% tlablist[[x]])})
  if(length(r2)==0) return(r2)
  tmp[which(r2)]
}

applyTreeText<-function(tdend,tdendA){
  lablist<-edgelablist(tdendA)
  dendrapply(tdend,function(x){
      tlab<-matchlabels(labels(x),lablist$inds)
      if(length(tlab)>0){
	attr(x,"edgetext")<-lablist$edgelab[tlab]
      }else attr(x,"edgetext")<-NULL
      x
  })
}

#############################

## FIXME: need larger par("mar")[1L] or [4L] for longish labels !
## {probably don't change, just print a warning ..}
fs.plot.dendrogram <-
    function (x, type = c("rectangle", "triangle"), center = FALSE,
	      edge.root = is.leaf(x) || !is.null(attr(x, "edgetext")),
	      nodePar = NULL, edgePar = list(),
	      leaflab = c("perpendicular", "textlike", "none"), dLeaf = NULL,
	      xlab = "", ylab = "", xaxt="n", yaxt="s",
	      horiz = FALSE, frame.plot = FALSE, xlim, ylim,height=0, ...)
{
    if(height>0)x<-dendrapply(x,function(x,h){attr(x, "height")<-attr(x, "height")+h;x},height)
    type <- match.arg(type)
    leaflab <- match.arg(leaflab)
    hgt <- attr(x, "height")
    if (edge.root && is.logical(edge.root))
	edge.root <- 0.0625 * if(is.leaf(x)) 1 else hgt
    mem.x <- .memberDend(x)
    yTop <- hgt + edge.root
    if(center) { x1 <- 0.5 ; x2 <- mem.x + 0.5 }
    else       { x1 <- 1   ; x2 <- mem.x }
    xl. <- c(x1 - 1/2, x2 + 1/2)
    yl. <- c(0, yTop)
    if (horiz) {## swap and reverse direction on `x':
	tmp <- xl.; xl. <- rev(yl.); yl. <- tmp
	tmp <- xaxt; xaxt <- yaxt; yaxt <- tmp
    }
    if(missing(xlim) || is.null(xlim)) xlim <- xl.
    if(missing(ylim) || is.null(ylim)) ylim <- yl.
    plot(0, xlim = xlim, ylim = ylim, type = "n", xlab = xlab, ylab = ylab,
	 xaxt = xaxt, yaxt = yaxt, frame.plot = frame.plot, ...)
    if(is.null(dLeaf))
        dLeaf <- .75*(if(horiz) strwidth("w") else strheight("x"))

    if (edge.root) {
### FIXME: the first edge + edgetext is drawn here, all others in plotNode()
### -----  maybe use trick with adding a single parent node to the top ?
	x0 <- fs.plotNodeLimit(x1, x2, x, center)$x
	if (horiz)
	    segments(hgt, x0, yTop, x0)
	else segments(x0, hgt, x0, yTop)
	if (!is.null(et <- attr(x, "edgetext"))) {
	    my <- mean(hgt, yTop)
	    if (horiz)
		text(my, x0, et)
	    else text(x0, my, et)
	}
    }
    fs.plotNode(x1, x2, x, type = type, center = center, leaflab = leaflab,
             dLeaf = dLeaf, nodePar = nodePar, edgePar = edgePar, horiz = horiz)
}

### the work horse: plot node (if pch) and lines to all children
fs.plotNode <-
    function(x1, x2, subtree, type, center, leaflab, dLeaf,
	     nodePar, edgePar, horiz = FALSE)
{
    inner <- !is.leaf(subtree) && x1 != x2
    yTop <- attr(subtree, "height")
    bx <- fs.plotNodeLimit(x1, x2, subtree, center)
    xTop <- bx$x

    ## handle node specific parameters in "nodePar":
    hasP <- !is.null(nPar <- attr(subtree, "nodePar"))
    if(!hasP) nPar <- nodePar

    if(getOption("verbose")) {
	cat(if(inner)"inner node" else "leaf", ":")
	if(!is.null(nPar)) { cat(" with node pars\n"); str(nPar) }
	cat(if(inner)paste(" height", formatC(yTop),"; "),
	    "(x1,x2)= (",formatC(x1,width=4),",",formatC(x2,width=4),")",
	    "--> xTop=", formatC(xTop, width=8),"\n", sep="")
    }

    Xtract <- function(nam, L, default, indx)
	rep(if(nam %in% names(L)) L[[nam]] else default,
	    length.out = indx)[indx]
    asTxt <- function(x) # to allow 'plotmath' labels:
	if(is.character(x) || is.expression(x) || is.null(x)) x else as.character(x)

    i <- if(inner || hasP) 1 else 2 # only 1 node specific par

    if(!is.null(nPar)) { ## draw this node
	pch <- Xtract("pch", nPar, default = 1L:2,	 i)
	cex <- Xtract("cex", nPar, default = c(1,1),	 i)
	col <- Xtract("col", nPar, default = par("col"), i)
	bg <- Xtract("bg", nPar, default = par("bg"), i)
	points(if (horiz) cbind(yTop, xTop) else cbind(xTop, yTop),
	       pch = pch, bg = bg, col = col, cex = cex)
    }

    if(leaflab == "textlike")
        p.col <- Xtract("p.col", nPar, default = "white", i)
    lab.col <- Xtract("lab.col", nPar, default = par("col"), i)
    lab.cex <- Xtract("lab.cex", nPar, default = c(1,1), i)
    lab.font <- Xtract("lab.font", nPar, default = par("font"), i)
    lab.crt <- Xtract("lab.crt", nPar, default = par("crt"), i)
    if (is.leaf(subtree)) {
	## label leaf
	if (leaflab == "perpendicular") { # somewhat like plot.hclust
	    if(horiz) {
                X <- yTop + dLeaf * lab.cex 
                Y <- xTop; srt <- 0; adj <- c(0, 0.5)
	    }
	    else {
                Y <- yTop - dLeaf * lab.cex 
                X <- xTop; srt <- 90; adj <- 1
	    }
            nodeText <- asTxt(attr(subtree,"label"))
	    text(X, Y, nodeText, xpd = TRUE, srt = srt, adj = adj,
                 cex = lab.cex, col = lab.col, font = lab.font,crt=lab.crt)
	}
    }
    else if (inner) {
	segmentsHV <- function(x0, y0, x1, y1) {
	    if (horiz)
		segments(y0, x0, y1, x1, col = col, lty = lty, lwd = lwd)
	    else segments(x0, y0, x1, y1, col = col, lty = lty, lwd = lwd)
	}
	for (k in seq_along(subtree)) {
	    child <- subtree[[k]]
	    ## draw lines to the children and draw them recursively
	    yBot <- attr(child, "height")
	    if (getOption("verbose")) cat("ch.", k, "@ h=", yBot, "; ")
	    if (is.null(yBot))
		yBot <- 0
	    xBot <-
		if (center) mean(bx$limit[k:(k + 1)])
		else bx$limit[k] + .midDend(child)

	    hasE <- !is.null(ePar <- attr(child, "edgePar"))
	    if (!hasE)
		ePar <- edgePar
	    i <- if (!is.leaf(child) || hasE) 1 else 2
	    ## define line attributes for segmentsHV():
	    col <- Xtract("col", ePar, default = par("col"), i)
	    lty <- Xtract("lty", ePar, default = par("lty"), i)
	    lwd <- Xtract("lwd", ePar, default = par("lwd"), i)
	    if (type == "triangle") {
		segmentsHV(xTop, yTop, xBot, yBot)
	    }
	    else { # rectangle
		segmentsHV(xTop,yTop, xBot,yTop)# h
		segmentsHV(xBot,yTop, xBot,yBot)# v
	    }
	    vln <- NULL
	    if (is.leaf(child) && leaflab == "textlike") {
		nodeText <- asTxt(attr(child,"label"))
		if(getOption("verbose"))
		    cat('-- with "label"',format(nodeText))
		hln <- 0.6 * strwidth(nodeText, cex = lab.cex)/2
		vln <- 1.5 * strheight(nodeText, cex = lab.cex)/2
		rect(xBot - hln, yBot,
		     xBot + hln, yBot + 2 * vln, col = p.col)
		text(xBot, yBot + vln, nodeText, xpd = TRUE,
                     cex = lab.cex, col = lab.col, font = lab.font)
	    }
	    if (!is.null(attr(child, "edgetext"))) {
		edgeText <- asTxt(attr(child, "edgetext"))
		if(getOption("verbose"))
		    cat('-- with "edgetext"',format(edgeText))
		if (!is.null(vln)) {
		    mx <-
			if(type == "triangle")
			    (xTop+ xBot+ ((xTop - xBot)/(yTop - yBot)) * vln)/2
			else xBot
		    my <- (yTop + yBot + 2 * vln)/2
		}
		else {
		    mx <- if(type == "triangle") (xTop + xBot)/2 else xBot
		    my <- (yTop + yBot)/2
		}
		## Both for "triangle" and "rectangle" : Diamond + Text

                p.col <- Xtract("p.col", ePar, default = "white", i)
                p.border <- Xtract("p.border", ePar, default = par("fg"), i)
                ## edge label pars: defaults from the segments pars
                p.lwd <- Xtract("p.lwd", ePar, default = lwd, i)
                p.lty <- Xtract("p.lty", ePar, default = lty, i)
                t.col <- Xtract("t.col", ePar, default = col, i)
                t.cex <- Xtract("t.cex", ePar, default =  1,  i)
                t.font<- Xtract("t.font",ePar, default= par("font"), i)
                t.srt<- Xtract("t.srt",ePar, default= 0, i)
				toffdefault<-0
				if(t.srt==90 || t.srt== -90) toffdefault<-0.15
				t.off<- Xtract("t.off",ePar, default= toffdefault, i)

				mx<-mx + t.off
				vlm <- strheight(c(edgeText,"h"), cex = t.cex)/2
		hlm <- strwidth (c(edgeText,"m"), cex = t.cex)/2
		hl3 <- c(hlm[1L], hlm[1L] + hlm[2L], hlm[1L])
                if(horiz) {
                    if(p.lwd>0) polygon(my+ c(-hl3, hl3), mx + sum(vlm)*c(-1L:1L, 1L:-1L),
                            col = p.col, border= p.border,
                            lty = p.lty, lwd = p.lwd)
                    text(my, mx, edgeText, cex = t.cex, col = t.col,
                         font = t.font,srt=t.srt)
                } else {
                     if(t.srt==90 || t.srt== -90) {
			mx <- mx -min( hlm)
		    }else{if(p.lwd>0) polygon(mx+ c(-hl3, hl3), my + sum(vlm)*c(-1L:1L, 1L:-1L),
                            col = p.col, border= p.border,
                            lty = p.lty, lwd = p.lwd)}
                    text(mx, my, edgeText, cex = t.cex, col = t.col,
                         font = t.font,srt=t.srt)
                }
	    }
	    fs.plotNode(bx$limit[k], bx$limit[k + 1], subtree = child,
		     type, center, leaflab, dLeaf, nodePar, edgePar, horiz)
	}
    }
    invisible()
}

fs.plotNodeLimit <- function(x1, x2, subtree, center)
{
    ## get the left borders limit[k] of all children k=1..K, and
    ## the handle point `x' for the edge connecting to the parent.
    inner <- !is.leaf(subtree) && x1 != x2
    if(inner) {
	K <- length(subtree)
	mTop <- .memberDend(subtree)
	limit <- integer(K)
	xx1 <- x1
	for(k in 1L:K) {
	    m <- .memberDend(subtree[[k]])
	    ##if(is.null(m)) m <- 1
	    xx1 <- xx1 + (if(center) (x2-x1) * m/mTop else m)
	    limit[k] <- xx1
	}
	limit <- c(x1, limit)
    } else { ## leaf
	limit <- c(x1, x2)
    }
    mid <- attr(subtree, "midpoint")
    center <- center || (inner && !is.numeric(mid))
    x <- if(center) mean(c(x1,x2)) else x1 + (if(inner) mid else 0)
    list(x = x, limit = limit)
}

#############################


################################
## Name manipulation functions - these may cause problems if your naming system isn't NAME<number> e.g. French12

NumToEnd<-function(x){
## moves numbers to the end of a name - careful as it might break things
	ymain<-gsub("[0-9]","",x)
	y<-gsub("[a-z]|[.,/_]","",x,ignore.case =TRUE)
	ret<-ymain
	ret[nchar(y)>0]<-paste(ymain,abs(as.numeric(y)),sep="")[nchar(y)>0]
	return(ret)
}


popIn<-function(x){
## which population labels are in a list of individuals
	tlab2<-gsub("[0-9]","",x)
	unique(tlab2)
}

popNos<-function(x,name) {
## part of the name expansion , removes name from the label and returns the numbers that this must represent
	x<-x[which(substr(x,1,nchar(name))==name)]
	x<-x[which(nchar(gsub("[0-9]","",substr(x,nchar(name)+1,nchar(name)+1)))==0)]
	tnames<-gsub("[0-9,]","",x)
	tmp1<-strsplit(tnames,"-")[[1]]
	tmp1<-tmp1[nchar(tmp1)>0]
	tnames<-paste(tmp1,collapse=".")
	x<-x[tnames==name]
	x<-gsub(name,"",x)
	x<-sort(as.numeric(x))
	x
}

NameCollapse<-function(x){
## Collapses individual names down to a list of numbers
	r1<-1
	ret<-NULL
	r2<-2
	if(r2>length(x)) ret<-x
	while(r2<=length(x)){
	  while(x[r2]-x[r1]==r2-r1) {r2<-r2+1;if(r2>length(x)) break;}
	  if(r2==r1+1){
		ret<-c(ret,x[r1])
	  }else{ ret<-c(ret,paste(x[r1],x[r2-1],sep="-"))}
	  r1<-r2
	}
	paste(ret,collapse=",")
}

getIndivsFromSummary2<-function(x){
## computes individuals within a single label from their name summary
	ret<-NULL
	tnames<-gsub("[0-9,]","",x)
	if(x==tnames) return(x)
	tmp1<-strsplit(tnames,"-")[[1]]
	tmp1<-tmp1[nchar(tmp1)>0]
	tmppop<-paste(tmp1,collapse=".")
	tmpdat<-substring(x,nchar(tmppop)+1)
	tmp1<-strsplit(gsub("[a-zA-Z._]","",x),",")[[1]]
	for(i in tmp1) {
		j<-strsplit(i,"-")[[1]]
		j<-j[nchar(j)>0]
		if(length(j)==1) ret<-c(ret,paste(tmppop,j,sep=""))
		else ret<-c(ret,paste(tmppop,seq(as.numeric(j[1]),as.numeric(j[2])),sep=""))
	}
	ret
}

getIndivsFromSummary<-function(x) {
## computes individuals in a population from their name summary
	tmp1<-strsplit(x,";")[[1]]
	as.vector(unlist(sapply(tmp1,getIndivsFromSummary2)))
}

NameExpand<-function(x){
	ret<-list()
	for(i in x){
		tlist<-NULL
		j<-strsplit(i,";")[[1]]
		for(k in j) {
			tlist<-c(tlist,getIndivsFromSummary(k))
		}
		ret[[i]]<-tlist
	}
	ret
} # takes a summary list of populations and returns the full set of individuals

NameExpandList<-function(fullpops){
    ## Run NameExpand on every population in a list (removing labels)
    clusterlist<-lapply(fullpops,function(x){
        ret<-unlist(NameExpand(x))
        names(ret)<-NULL
        ret
    })
    names(clusterlist)<-sapply(clusterlist,NameMoreSummary)
    clusterlist
}

NameSummary<-function(x){#x: a poplist element
	x<-sapply(x,NumToEnd)
	pin<-popIn(x)
	pnames<-lapply(pin,popNos,x=x)
	ret<-NULL
	for(i in 1:length(pin)) {
		ret<-c(ret,paste(pin[i],NameCollapse(pnames[[i]]),sep=""))
	}
	paste(ret,collapse=";")
}# returns summaries that completely represent the data

NameMoreSummary<-function(x){#x: a poplist element
	x<-sapply(x,NumToEnd)
	pin<-popIn(x)
	pnames<-lapply(pin,popNos,x=x)
	ret<-NULL
	for(i in 1:length(pin)) {
		if(length(pnames[[i]])>0){ret<-c(ret,paste(length(pnames[[i]]),pin[i],sep=""))
		}else ret<-c(ret,paste("1",pin[i],sep=""))
	}
	paste(ret,collapse=";")
} # returns a summary of the couints of each label only

summarizeNames<-function(l,omitbelow=0,fracbelow=0,addcount=FALSE){
    ## Heavily summarize population names from a list of names
    ## Specifically, making each only the names of the labels within
    ## Ordered by number of individuals from each label
    ## And then with A-Z added when duplicate names are created
    temppoplabels<-sapply(l,function(x){
        x<-gsub("[0-9,-]","",x)
        ttab<-sort(table(x),decreasing=TRUE)
        tcumsum<-c(1,1-cumsum(ttab)/sum(ttab))[1:length(ttab)]
        tkeep<-ttab>omitbelow & tcumsum>fracbelow
        tkeep[1]<-TRUE
        ttab<-ttab[tkeep]
        ret<-paste(names(ttab),collapse="_")
        if(any(!tkeep)) ret<-paste0(ret,"_etal")
        if(addcount) ret<-paste0(ret,sum(ttab))
        ret
    },USE.NAMES = FALSE)
    names(temppoplabels)<-NULL

    sapply(1:length(temppoplabels),function(i){
        x<-temppoplabels[i]
        if(sum(temppoplabels==x)>1){
            tw<-which(temppoplabels==x)
            return(paste(x,LETTERS[which(tw==i)],sep="_"))
        }
        return(x)
    },USE.NAMES = FALSE)
}

NameLessSummary<-function(x){
  paste(x,collapse=";")
}

FixNames<-function(x) {
## Some helpful modifications to make our names better behaved, do this to everything as you read it in
	x<-NumToEnd(x)
	x<-gsub("-",".",x)
	x
}

PopCenters<-function(popsizes){
## Gets the distance (on the 1:N scale) of the center of each population
	csum<-ret<-rep(0,length(popsizes))
	csum[1]<-popsizes[1]
	ret[1]<-(popsizes[1]+1)/2
	for(i in 2:length(popsizes)) {
		csum[i]<-csum[i-1]+popsizes[i]
		ret[i]<-csum[i-1] + (popsizes[i]+1)/2
	}
	ret
}

######################################
## PLOTTING FUNCTIONS
MakeColorYRP<-function(colby=0.05,final=NULL){
## makes yellow/red/purple colour scheme, adding rgb(final) if not null
	tmp<-c(rgb(1,seq(1,0,-colby),0),rgb(1,0,seq(colby,1,colby)),rgb(seq(1-colby,0,-colby),0,1.0))
	if(is.null(final)) return(tmp)
	c(tmp,rgb(final[1],final[2],final[3]))
}

## MAIN PLOTTING FUNCTION FOR HEATMAPS.  This function can be modified to get most desired results.
## You can add things to the main heatmap image after calling this function

writePopFile<-function(pl,file){ ## Create a "fake" mcmc file from a specified popualtion list pl, from which we can generate a tree as normal
  K<-length(pl)
  cat(file=file,
  "<?xml version = '1.0' encoding = 'UTF-8'?>
<outputFile>
<header>
<runtype>MCMC</runtype>
<inflation>0.0123730763</inflation>
<burnin>1</burnin>
<datafilename>NA</datafilename>
<copymodel>1</copymodel>
</header>
<comment>Created fom FinestructureLibrary.R</comment>
<Iteration>
<Pop>")
  cat(file=file,popListAsPop(pl),append = T)
  cat(file = file,append = T,
      paste0("</Pop>
<Posterior>-1</Posterior>
<K>",K,"</K>
<alpha>1</alpha>
<beta>1</beta>
<delta>1</delta>
<F>1</F>
<AccSAMS>0</AccSAMS>
<AccMS>0</AccMS>
<AccIndiv>0</AccIndiv>
<AccHyper>0</AccHyper>
<Number>0</Number>
</Iteration>
</outputFile>
"))
}

plotFinestructure<-function(tmpmat,# this is the heatmap to be plotted
		labelsx,# the x,y labels to be plotted
		labelsy=NULL,# assumed equal to labelsx unless otherwise specified
		labelsatx=NULL, # location of the X labels - note that the heatmap X ranges 1:dim(tmpmat)[1] and good locations are (1:dim(tmpmat)[1]) -0.5
		labelsaty=NULL, # assumed equal to labelsatx unless otherwise specified
		cols=NULL,# colour scale. Will be generated by MakeColors if NULL
		dend=NULL, # optional dendrogram to be placed at the top
		optpts=NULL, # optional points to add, e.g. MAP state
		labmargin=8, # margin for x,y labels
		layoutd=0.2, # proportion of the plot to be used by the dendrogram
		layoutf=0.1,# proportion of the plot to be used by the scale
		cex.axis=0.5, # cex.axis applies to the x,y labels
		xcrt=0, # character rotation of x,y labels (x is rotated the opposite direction)
		ycrt=0,
		colscale=NULL, # specify a colour scale, meaning the range of the z-values for the image
		text.col=NULL, # colour of the labels (can be a vector, is recycled
		ignorebelow=0, # if text.col=0, then we white out everything below ignorebelow
                ignoreabove=Inf, # Ignore values above this
		nodePar=list(cex=0,lab.cex=0.5,las=2),# nodePar as seen by plot.dendrogram
		edgePar=list(p.lwd=0,t.srt=0,t.off=0.2),# edgePar as seen by plot.dendrogram
		dendmar=c(0,0,2,1), # additional modification to the margins of the dendrogram
		scalemar=c(2,5,2,1), # additional modification to the margins of the scale
		hmmar=c(0,0,0,1), # additional modification to the margins of the heatmap
		cex.scale=1, # size of the characters in the scale
		scalelocs=NULL, # optional firced positions of the scale
		scalenum=10, # if scalelocs=NULL, the number of points to label on the scale
		scalesignif=3, # number of significant digits on scale
		scalelabel="", # label for the scale
		optcex=1.0, # size of the optional points
		optpch=20, # pch for the optional points
		optcol="black", # colour for the optional points
		labelsoff=c(1,1), # "margin" like, distance from colour/rotated labels and the axis
		tickmarks=1,# size of tickmarks when using colour/rotated labels
                startpt="bottomleft", # whether the first individual vs first individual is at the "topleft", or the "bottomleft"
                dolayout=TRUE,
                            main="", # dendrogram title
                            cex.main=1, # title size
                            adj=0 # title adj
                            ) 
{
if(dolayout) layout(matrix(c(2,1,4,3), 2, 2, byrow=TRUE),widths=c(1-layoutf,layoutf),height=c(layoutd))
if(is.null(cols)) cols<-MakeColorYRP()
## TOP RIGHT
if(dolayout || !is.null(dend))
  {
    par(mar=c(0,0,0,0))
    plot(c(0,1),c(0,1),type="n",axes=F,xlab="",ylab="")# null plot for top right
    
## TOP LEFT: DENDROGRAM
    par(mar=c(0,labmargin,0,0)+dendmar)
  }
if(dolayout && is.null(dend)){
	plot(c(0,1),c(0,1),type="n",axes=F,xlab="",ylab="",main=main,cex.main=cex.main,adj=adj)# null plot for top right
}else if(!is.null(dend)) {
  fs.plot.dendrogram(dend,horiz=FALSE,axes=FALSE,xaxs = "i", leaflab = "none",nodePar=nodePar,edgePar=edgePar,main=main,cex.main=cex.main,adj=adj)
}

## BOTTOM RIGHT: SCALE
if(is.null(colscale)) {
  ttmpmat<-na.omit(as.numeric(tmpmat))
  ttmpmat<-ttmpmat[ttmpmat>ignorebelow]
  ttmpmat<-ttmpmat[ttmpmat<ignoreabove]
  colscale<-c(max(ignorebelow,min(ttmpmat)),min(ignoreabove,max(ttmpmat)))
}
par(mar=c(labmargin,0,0,0)+scalemar)
colindex<-t(matrix(seq(min(colscale),max(colscale),length.out=100),ncol=1,nrow=100)) # colour scale
image(1,1:100,colindex,xaxt="n",yaxt="n",xlab=scalelabel,ylab="",col=cols,zlim=range(colindex))
if(is.null(scalelocs)){ # fill the range with scale labels
	scalelocs <- min(colindex)+(max(colindex)-min(colindex))*seq(0,1,length.out=scalenum)
	scalephysicalpos <- seq(1,100,length.out=scalenum)
}else{
	stop("NOT IMPLEMENTED...")
}
axis(2,at=scalephysicalpos,labels=signif(scalelocs,scalesignif),las=2,cex.axis=cex.scale)

## BOTTOM LEFT: MAIN HEATMAP
par(mar=c(labmargin,labmargin,0,0)+hmmar)
if(is.null(labelsatx)) { # labels x positions
	labelsatx<-seq(1:dim(tmpmat)[2]) 
}
if(is.null(labelsaty)) { # labels y positions
	labelsaty<-seq(1:dim(tmpmat)[1]) 
}
if(is.null(labelsy)) { # y labels
	labelsy<-labelsx
}
if(startpt!="bottomleft"){
  if(startpt!="topleft") stop("Invalid \"startpt\" option: valid choices are \"topleft\" and \"bottomleft\".")
  tmpmat<-tmpmat[dim(tmpmat)[1]:1,,drop=FALSE]
  if(!is.null(optpts)) {
    optpts<-optpts[dim(optpts)[1]:1,,drop=FALSE]
  }
  labelsaty<-dim(tmpmat)[1] - labelsaty + 1
}

#plot the heatmap
image(1:dim(tmpmat)[2],1:dim(tmpmat)[1],t(tmpmat),xaxt="n",yaxt="n",xlab="",ylab="",col=cols,zlim=colscale)
# draw the axes
if(xcrt==0 && ycrt==0 && is.null(text.col)) {
axis(1,at=labelsatx,labels=labelsx,las=2,cex.axis=cex.axis)
axis(2,at=labelsaty,labels=labelsy,las=2,cex.axis=cex.axis)
}else{
    text(labelsatx, -labelsoff[1], srt = xcrt, adj = 1,
          labels = labelsx, xpd = TRUE,cex=cex.axis,col=text.col)
    text(-labelsoff[2], labelsaty, srt = ycrt, adj = 1,
          labels = labelsy, xpd = TRUE,cex=cex.axis,col=text.col)
	if(tickmarks>0) {
		axis(1,at=labelsatx,labels=rep("",length(labelsatx)),las=2,cex.axis=cex.axis)
		axis(2,at=labelsaty,labels=rep("",length(labelsaty)),las=2,cex.axis=cex.axis)
	}

}
## add optional points
if(!is.null(optpts))tmp<-sapply(1:dim(optpts)[1],function(x){
  plist<-which(optpts[,x]==1);
  xpos<-rep(x,length(plist))
  ypos<-plist
#  if(yrev) ypos<-dim(tmpmat)[1] - plist + 1
  points(xpos,ypos,pch=optpch,cex=optcex,col=optcol)
  invisible(NULL)
})

## DONE

}

#############################################
## Things to do with PCA
mypcanorm<-function(x,zeromean=T,fdiag=T){
	## Normalise for PCA by setting the diagonal to be the average of the rows, then zero meaning by substracting the row mean (so the diagonal ends up as zero)
	if(fdiag){
		diag(x)<-rowSums(x)/(dim(x)[1]-1)
	}
	if(zeromean){
          return(x-rowMeans(x))
	}else{
		return(x)
	}
}
mypca<-function(x,zeromean=T,fdiag=T){
	## do PCA on x %*% t(x) , normalising first if desired
	xN<-mypcanorm(x,zeromean,fdiag)
	tcov<-xN %*% t(xN)
	eigen(tcov)	
}
optimalHorn<-function(x){
	require(paran)
	# Horns Parallel Analysis criterion
	tcov<-mypcanorm(x)
	return(paran(mypcanorm(tcov),quiet=T)$Retained)
}
optimalMap<-function(x){
	require(psych)
	# Velicer MAP criterion
	tcov<-mypcanorm(x)
		tcov<-tcov %*% t(tcov)
	tmap1<-VSS(tcov,n=20,rotate="none",fm="pc",n.obs=1000,plot=F)
	return(which(tmap1$map==min(tmap1$map)))
}


###############################################
## CHROMOPAINTER FUNCTIONS (currently limited)


plotdensitybetweensnps<-function(xstart,xend,dens,collist,bottom=0,rangemax=1,eps=1){
## This plots the density between two snps (at xstart and xend) as a box
# dens is the vector of densities
# collist is the vector of colours (should be the same length as dens)
# bottom is the where the vector will start (ignore unless using flexible plotting options)
# rangemax is the distance up from bottom that is plotted (ignore unless using flexible plotting options)
# eps is the "overlap" in x that is drawn (default is fine unless you use a distance that isn't SNP position)
      rangemax<-1
      for(ohap in 1:length(dens)) {
			tcopyprob<-dens[ohap]
			if(tcopyprob>0) rect(xstart-eps,bottom,xend+eps,bottom+rangemax*tcopyprob,density= -1,col=collist[ohap],border=NA)
			bottom<-bottom + rangemax*tcopyprob
      }
}

cpdensityplotInternal<-function(tlocs,copyprobs,collist,xmin=NULL, xmax=NULL,eps=1){
## Plots the chromopainter SNP probabilities, WITHOUT calling plot
# tlocs: the location of the SNPs, in base-pairs (makes most sense trerated as an integer in genome location)
# copyprobs: for each SNP, the copyprob output from chromopainter (an L by N matrix)
# collist is the vector of colours (should be the same length as dens)
# xmin, xmax: the beginning/end of the data range (assumed to be the SNP range if omitted)
# eps: the "overlap" between SNPs; for integers tlocs, use 1
	if(is.null(xmin))xmin<-min(tlocs)
	if(is.null(xmax))xmax<-max(tlocs)
	if(xmin < min(tlocs)) {
		plotdensitybetweensnps(xmin,tlocs[1],copyprobs[1,],collist=collist,eps=eps)
	}else{xmin<-min(tlocs)}
	for(snpon in 2:length(tlocs)){
		plotdensitybetweensnps(tlocs[snpon-1],tlocs[snpon],copyprobs[snpon-1,],collist=collist,eps=eps)
	}
	if(xmax > max(tlocs)) {
		plotdensitybetweensnps(tlocs[length(tlocs)],xmax,copyprobs[dim(copyprobs)[1],],collist=collist,eps=eps)
	}
}

cpdensityplot<-function(tlocs,copyprobs,collist=NULL,xmin=NULL, xmax=NULL,eps=1,xlab="SNP location",ylab="Density",...){
## Plots the chromopainter SNP probabilities, calling plot and drawing axes
# tlocs: the location of the SNPs, in base-pairs (makes most sense trerated as an integer in genome location)
# copyprobs: for each SNP, the copyprob output from chromopainter (an L by N matrix)
# xmin, xmax: the beginning/end of the data range (assumed to be the SNP range if omitted)
# eps: the "overlap" between SNPs; for integers tlocs, use 1
# xlab, ylab: plot labels
# ...: extra arguments passed to plot
	if(is.null(collist)) collist<-collist<-c(rgb(0,0,0),rgb(0.66,0,0),rgb(1,0,0),rgb(1,0,1),rgb(1,0.5,0),
	   rgb(0,0.5,0),rgb(0.25,0.8,0.25),rgb(0,0.5,0.75),rgb(0,0,1),rgb(0,0,0.5))
	if(length(collist)<dim(copyprobs)[2]) collist<-rep(collist,ceiling(dim(copyprobs)[2]/length(collist)))
	if(is.null(xmin))xmin<-min(tlocs)
	if(is.null(xmax))xmax<-max(tlocs)
	if(xmin>min(tlocs)) {
		twhich<-which(tlocs>xmin)
		tlocs<-tlocs[twhich]
		copyprobs<-copyprobs[twhich,]
	}
	if(xmax<max(tlocs)) {
		twhich<-which(tlocs<xmax)
		tlocs<-tlocs[twhich]
		copyprobs<-copyprobs[twhich,]
	}
	plot(c(xmin,xmax),c(0,1),type="n",xlab=xlab,ylab=ylab,axes=F,...)
	axis(1)
	axis(2)
	cpdensityplotInternal(tlocs,copyprobs,xmin,xmax,eps=eps,collist=collist)
}

numLines<-function(filename){
## Counts the number of lines in a file to pre-allocate data vector (not that fast!)
# filename : name of the file (can be g/bzipped)
	testconn <- file(filename, open="r") 
	csize <- 10000
	nolines <- 0
	while((readnlines <- length(readLines(testconn,csize))) >0 ) nolines<- nolines+readnlines
	close(testconn)
	nolines
}

getHap<-function(haplotype,filename,nlines=NULL,gzip=NULL,verbose=FALSE){
## Reads the copyprobsperlocus of a haplotype number from a file called filename (assumed .gz file if gzip=T, plain if gzip=F, detected from line ending if gzip=NULL)
# haplotype: the number of the desired haplotype. You must read each individually at present, it is best to create files with only the desired individuals in as they can get very large
# filename: filename of the chromopainter ".copyprobsperlocus.out.gz" file (can optionally be uncompressed first to ".copyprobsperlocus.out"
# nlines: the number of SNPs. if specified, this will speed the file read (by preallocating data)
# gzip: whether the file is gzipped. The default is to detect this from the filename (.gz ending)
# verbose: whether to print out which line is being processed (reasuring if the file is big, but some R versions may buffer this output)
	if(is.null(nlines))nlines<-numLines(filename)
	print(paste("Reading",nlines,"lines"))
	if(is.null(gzip)){
		if(length(grep(".gz$",filename)>0)) { gzip<-TRUE
		}else{ gzip<-FALSE}
	}
	if(gzip){zz <- gzfile(filename, "r")  # compressed file
	}else {zz <- file(filename, "r")}  # uncompressed file
	hapline<-paste("HAP",haplotype)
	if(verbose) lineon<-0
	if(verbose) {lineon<-lineon+1; print(paste("Line",lineon))}
	tres<-suppressWarnings(readLines(zz,1))
	while(tres!=hapline){
		tres<-suppressWarnings(readLines(zz,1,warn=FALSE))
#		if(verbose) {lineon<-lineon+1; print(paste("Skipping Line",lineon, "length", length(tres)))}
		if(length(tres)==0 || tres=="") {
			stop("File I/O error: haplotype not found in file!")
		}
	}
	## We are now at the start of a haplotype
	snpvec<-numeric(nlines) # maximum size
	snpon<-1
	while(TRUE){
		tline<-suppressWarnings(readLines(zz,1))
		if(length(tline)<1) break;
		tres<-strsplit(tline," ")[[1]]
		if(verbose) {lineon<-lineon+1; print(paste("Processing Line",lineon, "length", length(tres)))}
		if(length(tres)<=2) break;
		tres<-as.numeric(tres)
		snpvec[snpon]<-tres[1]
		if(snpon==1){
			probs<-matrix(0,nrow=nlines,ncol=length(tres)-1)
		}
		probs[snpon,]<-tres[-1]
		snpon<-snpon+1
	}
	close(zz)
	## remove unused lines
	keep<-which(snpvec>0)
	snpvec<-snpvec[keep]
	probs<-probs[keep,,drop=F]
	## reverse SNPS to get increasing order, and decide whether we need to add the 0 copy vector for self (depends on how the data was obtained)
        haveself<-any(colSums(probs)==0)
        if(haveself) return(list(snps=rev(snpvec),probs=probs[dim(probs)[1]:1,]))
	rprobs<-cbind(probs,rep(0,dim(probs)[1]))[dim(probs)[1]:1,]
	torder<-1:dim(probs)[2]
	torder[torder>=haplotype]<-torder[torder>=haplotype]+1
	torder<-c(torder,haplotype)
	rprobs<-rprobs[,order(torder)]
	rsnps<-rev(snpvec)
	return(list(snps=rsnps,probs=rprobs))
}
