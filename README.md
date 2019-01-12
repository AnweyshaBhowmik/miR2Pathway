# miR2Pathway
miR2Pathway is a PageRank-based method that can be used to rank disease risk of miRNA-mediated biological pathways. miR2Pathway can help explore how much miRNAs differentially influence the activity of biological pathways between two classes of phenotypes.

<b> For information about the method, please read: </b>

Chaoxing Li and Valentin Dinu. “miR2Pathway: A Novel Analytical Method to Discover MicroRNA-mediateding Dysregulated Pathways Involved in Disease”, 2017.

Main website: http://dinulab.org/tools/mir2pathway/

<b> R code, sample data and sample scripts: </b>

The R script to perform miR2Pathway: <samp>   <br> miR2Pathway.R</samp>

Data used in miR2Pathway_Example_Dataset.RData: <samp> <br>   miR2Pathway_Example_Dataset.RData</samp>



# Usage

<b>******* Please use miR2Pathway with R x64 3.5.0 Pre-release *******</b>
  
<b> Example Usage (to be used with example data) </b>  

library(compiler) <br>
library(foreach) <br>
library(iterators) <br>
library(parallel) <br>
library(doParallel) <br>
library(AnnotationDbi) <br>
library(stats4) <br>
library(BiocGenerics) <br>
library(miRNAtap)  <br>
library(miRNAtap.db)  <br>
library(graphite) <br>
library(BiocGenerics) <br>
library(graph) <br>
library(igraph) <br>

<code>humanKEGG <- pathways("hsapiens", "kegg")</code>

<b>If you wish, you may download the example dataset files individually and run the code below to pre-process them <br> </b>

<code>genelist <- read.table("genelist.txt", header = TRUE, sep = "\t")</code><br>
<code>miRlist <- read.table("miRlist.txt", header = TRUE, sep = "\t")</code><br>
<code>miRlist.full <- read.table("miRlist.full", header = TRUE, sep = "\t")</code><br>

<code>mydata.miR <- read.table("mydata.miR", header = TRUE, sep = "\t")</code><br>
<code>names.genelist <- read.table("names.genelist", header = TRUE, sep = "\t")</code><br>
<code>mydata.gene <- read.table("mydata.gene", header = TRUE, sep = "\t")</code><br>

<code>mydata.gene2 <- mydata.gene[,-1]</code><br>
<code>rownames(mydata.gene2) <- mydata.gene[,1]</code><br>

<code>mydata.miR2 <- mydata.miR[,-1]</code><br>
<code>rownames(mydata.miR2) <- mydata.miR[,1]</code><br>

<code>names.genelist2 <- names.genelist[,-1]</code><br>
<code>rownames(names.genelist2) <- names.genelist[,1]</code><br>
 
<b>Or you may simply load the workspace file: <code>load("miR2Pathway_Example_Dataset.RData")</code> <br></b>
<b>Then run the following code: </b><br>

<code>results<-miR2Pathway(mydata.gene=mydata.gene2,mydata.miR=mydata.miR2,genelist=genelist,name.genelist=names.genelist2,miRlist=miRlist,miRlist.full=miRlist.full,N.miR=4,N.gene=4,N.path=4,Num.sample.normal=4,Num.sample.case=4,Pathway.database=humanKEGG,cor.cutoff=-0.4,N.parallel=4)</code> 
   
   

# Arguments

<code> miR2Pathway <- function(mydata.gene,mydata.miR,genelist,name.genelist,miRlist,miRlist.full,N.miR,N.gene,N.path,Num.sample.normal,Num.sample.case,Pathway.database,cor.cutoff,N.parallel) </code>

<html>

<body>

<table>
    <tr>
    <th>Object</th>
     <th>Description</th>
   </tr>
  <tr>
   <tr>
    <td>mydata.gene</td>
    <td>a gene expression dataset (matrix). Rows represent genes, and columns represent samples (from control to case).</td>
  </tr>
  <tr>
    <td>mydata.miR</td>
    <td>a miRNA expression dataset (matrix). Rows represent miRNAs, and columns represent samples (from control to case).</td>
  </tr>
  <tr>
    <td>genelist</td>
    <td>a list of gene names (e.g., gene symobol).</td>
  </tr>
  <tr>
    <td>name.genelist</td>
    <td>a matrix of gene expression dataset (The order of gene names must be consistent with genelist). This matrix is 1*N, where N is the number of genes. The value of each element could be a ramdom value. See an example in Examples Section. </td>
  </tr>
  <tr>
    <td>miRlist</td>
    <td>miRNA in a standard format</td>
  </tr>
  <tr>
    <td>miRlist.full</td>
    <td>the full name of all miRNAs.</td>
  </tr>
</table>

</body>
</html>

  

# Values

<html>
<body>

<table>
   <tr>
    <th>Object</th>
    <th>Description</th>
   </tr>
  <tr>
  <td>P.Value</td>
  <td>the p-value that resulted from the correlation test</td>
  </tr>
   <tr>
  <td>Pathways</td>
  <td>the reported pathways</td>
  </tr>
  <tr>
  <td>T.score</td>
  <td>the total differential influence of all the miRNAs on the activity of a single pathway between control and case</td>
  </tr>
  <tr>
   <td>LengthOfPathway</td>
   <td>gene counts of pathways</td>
  </tr>
</table>

</body>
</html>      

# Examples

See miR2Pathway_Example_Dataset.RData <br>
An example for the matrix of name.genelist:

<html>
<body>

<table>
   <tr>
    <th>samples</th>
    <th>KIF3B</th>
    <th>SP4</th>
    <th>DEF8</th>
    <th>USP9X</th>
   </tr>
   <tr>
    <td>sample1</td>
    <td>50 </td>
      <td>78</td>
      <td>87</td>
      <td>80</td>
    </tr>
    <tr>
    <td>sample2</td>
    <td>70</td>
      <td>87</td>
      <td>64</td>
      <td>55</td>
    </tr>
    <tr>
    <td>sample3</td>
    <td>90</td>
      <td>96</td>
      <td>41</td>
      <td>30</td>
    </tr>
   <tr>
    <td>sample4</td>  	 	 	 
    <td>110</td>
      <td>105</td>
      <td>18</td>
      <td>5</td>
    </tr>
    <tr>
    <td>sample5</td>  	 	 	 
    <td>130</td>
      <td>114</td>
      <td>4</td>
      <td>32</td>
    </tr>
    <tr>
    <td>sample6</td>  	 	 	 
    <td>150</td>
      <td>123</td>
      <td>124</td>
      <td>74</td>
    </tr>
   <tr>
    <td>sample7</td>  	 	 	 
    <td>170</td>
      <td>132</td>
      <td>54</td>
      <td>33</td>
    </tr>
   <tr>
    <td>sample8</td>  	 	 	 
    <td>190</td>
      <td>141</td>
      <td>53</td>
      <td>30</td>
    </tr>
</table>

</body>
</html>

 
# Credits

<b> Authors: </b> Chaoxing Li and Valentin Dinu

<b> If you use or modify the code, please cite: </b>

“miR2Pathway: A Novel Analytical Method to Discover MicroRNA-mediateding Dysregulated Pathways Involved in Disease”, 2017.

<b> Issues? </b>

Please email Chaoxing Li <chaoxing@asu.edu> or Valentin Dinu <Valentin.Dinu@asu.edu> if you have any questions, concerns or feedback.
