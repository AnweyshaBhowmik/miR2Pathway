#!/usr/bin/env  
args = commandArgs(trailingOnly=TRUE)

## Program: miR2Pathway Official
## Date: 9/19/2018

##Rscript miR2Pathway.r inputs... 


# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Please provide at least one argument (input file).", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[7] = "miR2Pathway_Results.txt"
}


miR2Pathway <- function(mydata.gene,mydata.miR,genelist,name.genelist,miRlist,miRlist.full,N.miR,N.gene,N.path,Num.sample.normal,Num.sample.case,Pathway.database,cor.cutoff,N.parallel){
  
  require(compiler) 
  require(foreach) 
  require(iterators)
  require(parallel)
  require(doParallel)
  require(AnnotationDbi)
  require(stats4)
  require(BiocGenerics)
  require(miRNAtap) 
  require(miRNAtap.db) 
  require(graphite)
  require(BiocGenerics)
  require(graph)
  require(igraph)
  require(org.Hs.eg.db)
  require(topGO)
  require(dplyr)
  options(warn=-1)
  cl <- makeCluster(N.parallel)
  registerDoParallel(cl)
  getDoParWorkers()
  
  amiRPath<-lapply(rep(0,1),function(x) lapply(rep(0,N.miR),function(x) graph.empty(x,directed=FALSE)))
  
  amiRPath.normal1<-amiRPath
  amiRPath.tumor1<-amiRPath
  amiRPath.normal<-amiRPath
  amiRPath.tumor<-amiRPath
  length.pathway<-c()
  title.pathway<-c()
  PR.path.sum<-c()
  fromResults3 <- c()
  humanReactome <- pathways("hsapiens", "reactome")
  humanBiocarta <- pathways("hsapiens", "biocarta")
  humanKEGG <- pathways("hsapiens", "kegg")
  
  miRtarnetwork<-lapply(rep(0,N.miR),function(x) graph.empty(x,directed=FALSE))
  
	##Order miRNA list
	
	miRlist.1 = as.data.frame(miRlist)
	miRlist.2 = as.character(miRlist.1[order(miRlist.1$miRlist),])
	print("miR2Pathway is now beginning the analysis")
	
	
	##Use lapply with the getPredictedTargets function and save output to the targetgenes variable

	targets <- function(w,x,y,z) {
				  results = getPredictedTargets(w, species = x, method = y, min_src = z)
				  return(results)
	}
	
	targetgenes.1 <- lapply(miRlist.2, targets, 'hsa', 'geom', 2)

	##Extract the rownames

	targetgenes.2 <- lapply(targetgenes.1, rownames)

	##Extract each dataframe from the targetgenes.2 object and concatentate into one numeric matrix

	tg.name.1 = as.matrix(do.call(rbind, targetgenes.2))
	columns = ncol(tg.name.1)
	tg.name.2 = c(t(as.matrix(do.call(rbind, targetgenes.2))))
	

	##Replicate the sorted miRNAs so that their dimentions match those of the target genes 
        miRlist.full.1 = as.data.frame(miRlist.full)
	miRlist.full.2 = as.character(miRlist.full.1[order(miRlist.full.1$miRlist.full),])
	
	
	miRlist.full.3 = rep(miRlist.full.2, columns)
	miRlist.full.4 = as.data.frame(miRlist.full.3)
	miRlist.full.5 = as.character(miRlist.full.4[order(miRlist.full.4$miRlist.full.3),])
	
	edgelist.1 <- as.matrix(miRlist.full.5)
	edgelist.2 <- as.matrix(tg.name.2)

	edgelist <- cbind(edgelist.1, edgelist.2)

	##Create the edges 

	miRtarnetwork <- graph_from_edgelist(edgelist,direct=FALSE)
	print("The edges for the miRNAs have been created")
	print("The results are ready, please use ls() to find the results object")
 
   
  
  cor.normal1<-matrix(0,1,N.gene)
 
  cor.normal<-foreach (i=1:N.miR, .combine="c") %dopar% {
    
    for (j in 1:N.gene){
      
      cor1<-cor.test(as.numeric(mydata.gene[j,1:Num.sample.normal]),as.numeric(mydata.miR[i,1:Num.sample.normal]))[cor.test(as.numeric(mydata.gene[j,1:Num.sample.normal]),as.numeric(mydata.miR[i,1:Num.sample.normal]))[[4]]<(cor.cutoff)]$p.value
      
      if (length(as.numeric(cor1))==0){
        cor.normal1[1,j]<- 1
      }else{
        cor.normal1[1,j]<- cor1
      }
      
    }
    return(cor.normal1)
  
  dim(cor.normal)=c(N.miR,N.gene)
  
  cor.normal.adj<-as.matrix(cor.normal)         #convert data.frame to matrix
  
  ##convert a matrix into adjacency matrix for each pathway.
  
  result1<-foreach (y=1:N.path, .combine="c") %dopar% {
    require(graphite)
    require(igraph)
    print(y)
    
    amiRPath.normal<-amiRPath
    
    KEGG.p0 <-Pathway.database[[y]]      
    KEGG.p <- convertIdentifiers(KEGG.p0, "entrez") 
    KEGG.g<-pathwayGraph(KEGG.p) 
    nodelist<-nodes(KEGG.g)
    KEGG.g1<-igraph.from.graphNEL(KEGG.g) 
    KEGG.g2<-graph_from_edgelist(get.edgelist(KEGG.g1),direct=FALSE)
    print(y)
    length.pathway[y]<-length(nodelist)
    
    loc <- c()
    for (i in 1:length(nodes(KEGG.p0))){
      ##probably can't find KEGG.p0 genes in genelist.
      if (is.element(nodes(KEGG.p0)[i],unlist(genelist))==FALSE){
        next
      }else{
        loc[i] <-which(as.vector(genelist)[1]==nodes(KEGG.p0)[i])
        
      }
    }
    
    loc<-loc[!is.na(loc)]
    loc<-loc[!is.nan(loc)]
    
    cor.normal.adj.path<-matrix(2,N.miR,length(loc))
    
    cor.normal.adj.path<-cor.normal.adj[,loc] 
    
    cor.normal.adjacency.path<-matrix(2,(N.miR+length(loc)),(N.miR+length(loc)))
    
    colnames(cor.normal.adjacency.path)[(N.miR+1):(N.miR+length(loc))]<- names(name.genelist[loc])       
    rownames(cor.normal.adjacency.path)[(N.miR+1):(N.miR+length(loc))]<- names(name.genelist[loc]) 
    colnames(cor.normal.adjacency.path)[1:N.miR]<-miRlist.full[1:N.miR]  
    rownames(cor.normal.adjacency.path)[1:N.miR]<-miRlist.full[1:N.miR] 
    
    for (j in 1:N.miR){
      cor.normal.adjacency.path[j,(N.miR+1):(N.miR+length(loc))]<- cor.normal.adj.path[j,]       ##length(loc) take over of length(nodes(KEGG.p0))
      
    }
    
    cor.normal.adjacency.path[cor.normal.adjacency.path>0.05]<-0
    cor.normal.adjacency.path[ is.na(cor.normal.adjacency.path)] <- 0
    cor.normal.adjacency.path[ is.nan(cor.normal.adjacency.path)] <- 0
    diag(cor.normal.adjacency.path) <- 0
 
    path.interaction.normal<-graph.empty()
    
    path.interaction.normal<-graph.adjacency(cor.normal.adjacency.path,weighted=TRUE,mode="undirected")
    
    miRPath.normal.stat<-graph.union(path.interaction.normal,KEGG.g2) 
    
    miRPath.normal.stat<-delete_vertices(miRPath.normal.stat, which(degree(miRPath.normal.stat)==0))   ##remove isolated nodes
    
    
    miRPath.normal.stat.miRlist<-delete_vertices(miRPath.normal.stat, nodes(KEGG.p0))   ##remove pathway's genes only leave miRNA nodes
    
    miR.left<-V(miRPath.normal.stat.miRlist)
    
    length.miR<-length(miR.left)
    
    if (length.miR==0){
      
      print(c("This pathway is not targeted by any miRNAs"))
      return(amiRPath)
      next
      
    }
    else{
      
      ##Find the locations of miR.left in name.miR(1046)
      loc.miR<-c()
      miRPath.normal.stat1<-get.adjacency(miRPath.normal.stat)
      miRPath.normal.stat<-graph.adjacency(miRPath.normal.stat1,weighted=TRUE,mode="undirected")
      
      for (i in 1:length.miR){
        loc.miR[i]<-which(miRlist.full==names(miR.left[i]))
      }
      
      for (x in c(loc.miR)){ ##The numebr of left miRNAs
        
        amiRPath.normal[[1]][[x]]<- induced.subgraph(miRPath.normal.stat,c(miRlist.full[x],names(V(KEGG.g2)))) ##names(miR.left) is full name.
        
        graph.intersect.normal<-graph.intersection(miRtarnetwork[[x]],amiRPath.normal[[1]][[x]])
        
        graph.intersect.normal<-delete_vertices(graph.intersect.normal, which(degree(graph.intersect.normal)==0)) 
        
        amiRPath.normal[[1]][[x]]<-graph.union(graph.intersect.normal,KEGG.g2) 
        
        if (length(V(amiRPath.normal[[1]][[x]]))==length(V(KEGG.g2))){
          print(c("Pathway",y,"miR",x,"has no connection"))
          amiRPath.normal[[1]][[x]]<-graph.empty()
        } 
      }
    }
    return(amiRPath.normal)
  } 
  } 
  ##Construct a miR-gene network for tumor
  
  cor.tumor1<-matrix(0,1,N.gene)
  
  cor.tumor<-foreach (i=1:N.miR, .combine="c") %dopar% {
    
    for (j in 1:N.gene){
      
      cor1<-cor.test(as.numeric(mydata.gene[j,(Num.sample.normal+1):(Num.sample.normal+Num.sample.case)]),as.numeric(mydata.miR[i,(Num.sample.normal+1):(Num.sample.normal+Num.sample.case)]))[cor.test(as.numeric(mydata.gene[j,(Num.sample.normal+1):(Num.sample.normal+Num.sample.case)]),as.numeric(mydata.miR[i,(Num.sample.normal+1):(Num.sample.normal+Num.sample.case)]))[[4]]<(cor.cutoff)]$p.value
      
      if (length(as.numeric(cor1))==0){
        cor.tumor1[1,j]<- 1
      }else{
        cor.tumor1[1,j]<- cor1
      }
      
    }
    return(cor.tumor1)
  
  dim(cor.tumor)=c(N.miR,N.gene)
  
  cor.tumor.adj<-as.matrix(cor.tumor)        #convert data.frame to matrix
  
  ##convert a matrix into adjacency matrix for each pathway.
  
  result2 <- foreach (y=1:N.path, .combine="c") %dopar% {
    require(graphite)
    require(igraph)
    print(y)
    
    amiRPath.tumor<-amiRPath
    
    KEGG.p0 <-Pathway.database[[y]]      
    KEGG.p <- convertIdentifiers(KEGG.p0, "entrez") 
    KEGG.g<-pathwayGraph(KEGG.p) 
	
    KEGG.g1<-igraph.from.graphNEL(KEGG.g) 
    KEGG.g2<-graph_from_edgelist(get.edgelist(KEGG.g1),direct=FALSE)
    print(y)
    
    loc <- c()
    for (i in 1:length(nodes(KEGG.p0))){
      ##probably can't find KEGG.p0 genes in genelist.
      if (is.element(nodes(KEGG.p0)[i],unlist(genelist))==FALSE){
        next
      }else{
        loc[i] <-which(as.vector(genelist)[1]==nodes(KEGG.p0)[i])
        
      }
    }
    
    loc<-loc[!is.na(loc)]
    loc<-loc[!is.nan(loc)]
    
    cor.tumor.adj.path<-matrix(2,N.miR,length(loc))
    
    cor.tumor.adj.path<-cor.tumor.adj[,loc] 
    
    cor.tumor.adjacency.path <- matrix(2,(N.miR+length(loc)),(N.miR+length(loc)))
	
    colnames(cor.tumor.adjacency.path)[(N.miR+1):(N.miR+length(loc))] <- names(name.genelist[loc])       
    rownames(cor.tumor.adjacency.path)[(N.miR+1):(N.miR+length(loc))] <- names(name.genelist[loc]) 
    colnames(cor.tumor.adjacency.path)[1:N.miR] <- miRlist.full[1:N.miR]  
    rownames(cor.tumor.adjacency.path)[1:N.miR] <- miRlist.full[1:N.miR] 
     
	for (j in 1:N.miR){
      cor.tumor.adjacency.path[j,(N.miR+1):(N.miR+length(loc))]<- cor.tumor.adj.path[j,]       ##length(loc) take over of length(nodes(KEGG.p0))
    }
    
    cor.tumor.adjacency.path[cor.tumor.adjacency.path>0.05]<-0
    cor.tumor.adjacency.path[ is.na(cor.tumor.adjacency.path)] <- 0
    cor.tumor.adjacency.path[ is.nan(cor.tumor.adjacency.path)] <- 0
    diag(cor.tumor.adjacency.path) <- 0
    
    path.interaction.tumor<-graph.empty()
    
    path.interaction.tumor<-graph.adjacency(cor.tumor.adjacency.path,weighted=TRUE,mode="undirected")
    
    miRPath.tumor.stat<-graph.union(path.interaction.tumor,KEGG.g2) 
    
    miRPath.tumor.stat<-delete_vertices(miRPath.tumor.stat, which(degree(miRPath.tumor.stat)==0))   ##remove isolated nodes
    
    miRPath.tumor.stat.miRlist<-delete_vertices(miRPath.tumor.stat, nodes(KEGG.p0))   ##remove pathway's genes only leave miRNA nodes
    
    miR.left<-V(miRPath.tumor.stat.miRlist)
    
    length.miR<-length(miR.left)
    
    if (length.miR==0){
      
      print(c("This pathway is not targeted by any miRNAs"))
      return(amiRPath)
      next
      
    }else{
      
      ##Find the locations of miR.left in name.miR(1046)
      loc.miR<-c()
      miRPath.tumor.stat1<-get.adjacency(miRPath.tumor.stat)
      miRPath.tumor.stat<-graph.adjacency(miRPath.tumor.stat1,weighted=TRUE,mode="undirected")
      
      for (i in 1:length.miR){
        loc.miR[i]<-which(miRlist.full==names(miR.left[i]))
      }
      
      for (x in c(loc.miR)){ ##The number of left miRNAs
        
        amiRPath.tumor[[1]][[x]]<- induced.subgraph(miRPath.tumor.stat,c(miRlist.full[x],names(V(KEGG.g2)))) ##names(miR.left) is full name.
        
        graph.intersect.tumor<-graph.intersection(miRtarnetwork[[x]],amiRPath.tumor[[1]][[x]])
        
        graph.intersect.tumor<-delete_vertices(graph.intersect.tumor, which(degree(graph.intersect.tumor)==0)) 
        
        amiRPath.tumor[[1]][[x]]<-graph.union(graph.intersect.tumor,KEGG.g2) 
        
        if (length(V(amiRPath.tumor[[1]][[x]]))==length(V(KEGG.g2))){
          print(c("Pathway",y,"miR",x,"has no connection"))
          amiRPath.tumor[[1]][[x]]<-graph.empty()
        } 
      }
    }
    return(amiRPath.tumor)
  } 
  
  result3<-foreach (y=1:N.path, .combine="c") %dopar% {
    require(graphite)
    require(igraph)
    print(y)
    PR.miR.sum<-matrix(0,N.miR,1)
    
    KEGG.p0 <-Pathway.database[[y]]      
    KEGG.p <- convertIdentifiers(KEGG.p0, "entrez") 
    KEGG.g<-pathwayGraph(KEGG.p) 
    KEGG.g1<-igraph.from.graphNEL(KEGG.g) 
    KEGG.g2<-graph_from_edgelist(get.edgelist(KEGG.g1),direct=FALSE)
    
    for (x in 1:N.miR){
      mat1<-get.adjacency(result1[[y]][[x]])
      mat2<-get.adjacency(result2[[y]][[x]])
      if (dim(mat1)==dim(mat2)){
        if (dim(mat1)==0){
          print("diffenetwork is empty")
          PR.miR.sum[x,1]<-0
          next
        }
        mat3<-mat1-mat2
        
        diffnetwork <- graph.adjacency(mat3,weighted=TRUE,mode="undirected")
        
        neigh<-neighbors(diffnetwork,v=miRlist.full[[x]])
        PR.miR.sum[x,1]<-sum(page.rank(KEGG.g2)$vector[neigh])
        
      }else{
        if (dim(mat1)==0){    
          diffnetwork<-graph.adjacency(mat2,weighted=TRUE,mode="undirected")
          
          neigh<-neighbors(diffnetwork,v=miRlist.full[[x]])
          PR.miR.sum[x,1]<-sum(page.rank(KEGG.g2)$vector[neigh])
          
        }else{
          diffnetwork<-graph.adjacency(mat1,weighted=TRUE,mode="undirected")
          
          neigh<-neighbors(diffnetwork,v=miRlist.full[[x]])
          PR.miR.sum[x,1]<-sum(page.rank(KEGG.g2)$vector[neigh])
        }
      }
    }
    return(PR.miR.sum)
  } ## end of the results3 loop
 
  ## Extract the T scores
  ## Note: results3 is not NULL. The following line of code would not work if results3 were NULL
  dim(results3)=c(N.miR,N.path)
  fromResults3 <<- results3
} ## end of the cor loop
  
  ## Extract the pathway names
  for (xt in 1:length(Pathway.database)) {
    p0 <-Pathway.database[[xt]]
    title.pathway[xt] <- p0@title
  } 

  ## Extract the T scores
  for (y in 1:length(Pathway.database))  {
    PR.path.sum[y]<-sum(fromResults3[,y])
  } 
  
  
#write.table(fromResults3, "fromResults3.txt", sep=",", quote=F)
#write.table(PR.path.sum, "T_scores.txt", sep=",", quote=F)
#write.table(sprintf("%.05f",cor.tumor), "P_values.txt", sep=",", quote=F)
#write.table(title.pathway, "Pathways.txt", sep=",", quote=F)
Full_Results <- list(Pathways=title.pathway, P.Values=sprintf("%.05f",cor.tumor), T.scores= PR.path.sum)
}


library(compiler) 
library(foreach) 
library(iterators)
library(parallel)
library(doParallel)
library(AnnotationDbi)
library(stats4)
library(BiocGenerics)
library(miRNAtap) 
library(miRNAtap.db) 
library(graphite)
library(BiocGenerics)
library(graph)
library(igraph)
library(org.Hs.eg.db)
library(topGO)
library(dplyr)
####### INPUT ARGS KEY ##############
# mydata.gene = mydata.gene2             args[1]
# mydata.miR  = mydata.miR2              args[2]
# genelist = genelist                    args[3]
# name.genelist = names.genelist2        args[4]
# miRlist = miRlist                      args[5]
# miRlist.full = miRlist.full            args[6]
# N.miR = nrow(miRlist)
# N.gene = nrow(genelist)
# N.path = length(humanKEGG)
# Num.sample.normal =  ncol(select(mydata.gene, contains("normal")))
# Num.sample.case = ncol(select(mydata.gene, contains("tumor")))
# humanKEGG <- pathways("hsapiens", "kegg")
# Pathway.database = humanKEGG
# cor.cutoff=-0.4
# N.parallel=4

humanKEGG <- pathways("hsapiens", "kegg")
 
mydata.gene = read.table(args[1], sep="\t", header=TRUE)
mydata.miR  = read.table(args[2], sep="\t", header=TRUE)
genelist = read.table(args[3], sep="\t", header=TRUE)
name.genelist = read.table(args[4], sep="\t", header=TRUE)
miRlist = read.table(args[5], sep="\t", header=TRUE)
miRlist.full = read.table(args[6], sep="\t", header=TRUE)


N.miR = nrow(miRlist)
N.gene = nrow(genelist)
N.path = length(humanKEGG)
Num.sample.normal =  ncol(select(mydata.gene, contains("normal")))
Num.sample.case = ncol(select(mydata.gene, contains("tumor")))
Pathway.database = humanKEGG
cor.cutoff=-0.4
N.parallel=4

results<-miR2Pathway(mydata.gene=mydata.gene,mydata.miR=mydata.miR,genelist=genelist,name.genelist=name.genelist,miRlist=miRlist,miRlist.full=miRlist.full,N.miR=N.miR,N.gene=N.gene,N.path=N.path,Num.sample.normal=Num.sample.normal,Num.sample.case=Num.sample.case,Pathway.database=humanKEGG,cor.cutoff=cor.cutoff,N.parallel=N.parallel)

#### Post-Processing


Pathways <- as.data.frame(unlist(results[1]))
names(Pathways) <- c("Pathways")
rownames(Pathways) <- NULL
PathM <- as.matrix(Pathways)

P.values <- as.data.frame(unlist(results[2]))
names(P.values) = c("P.values")
rownames(P.values) <- NULL
P.valM <- as.matrix(P.values)

max.len <- max(length(PathM), length(P.valM))

P.values <- as.data.frame(c(P.valM, rep(1, max.len - length(P.valM))))
names(P.values) = c("P.values")

T.scores = as.data.frame(results[3])
Rank = as.data.frame(c(1:301))
names(Rank)  <- c("Rank")

df1 <- data.frame(Pathways, P.values, T.scores)
df2 <- df1[order(df1$P.values),]
df3 <- data.frame(df2,Rank)
df4 <- df3[order(df3$Pathways),]

df_sub <- as.data.frame(gsub(" ", "_", df4$Pathways)) ##args[]
names(df_sub) <- c("Pathways")
results$Pathways = df_sub$Pathways

write.table(results, args[7], sep="\t", quote=FALSE, row.names=FALSE)
