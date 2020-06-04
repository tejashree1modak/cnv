#!/usr/bin/env Rscript
#source("http://bioconductor.org/biocLite.R")
#biocLite("topGO")
library("topGO")
#biocLite("Rgraphviz")
library("Rgraphviz")


geneID2GO <- readMappings(file = "/Users/tejashree/Documents/Projects/cnv/scripts/output_files/oyster_cnv/topgo/annotations_topgo.txt") 
genesOfInterest <- read.table("/Users/tejashree/Documents/Projects/cnv/scripts/output_files/oyster_cnv/topgo/interestinggenes.txt",header=FALSE)
geneUniverse <- names(geneID2GO)
genesOfInterest <- as.character(genesOfInterest$V1)
head(genesOfInterest)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myBPGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
myMFGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
myCCGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
#Look at the object
myBPGOdata
myMFGOdata
myCCGOdata
#Access significant genes
#BP
sg_bp <- sigGenes(myBPGOdata)
str(sg_bp)
numSigGenes(myBPGOdata) 
#MF
sg_mf <- sigGenes(myMFGOdata)
str(sg_mf)
numSigGenes(myMFGOdata) 
#CC
sg_cc <- sigGenes(myCCGOdata)
str(sg_cc)
numSigGenes(myCCGOdata)

############ Performing enrichment tests ############ 
#Fishers exact test
#ontology = BP
resultFisherBP <- runTest(myBPGOdata, algorithm="classic", statistic="fisher") 
resultFisherBP
#Selecting 'algorithm=classic' means that the GO hierarchy isn't taken into account, 
#so each GO term is tested independently.
resultFisherBP #The p-values have not been corrected for multiple testing.
#ontology = MF
resultFisherMF <- runTest(myMFGOdata, algorithm="classic", statistic="fisher") 
resultFisherMF #The p-values have not been corrected for multiple testing.
#ontology = CC
resultFisherCC <- runTest(myCCGOdata, algorithm="classic", statistic="fisher") 
resultFisherCC

#The GenTable function returns a data frame containing the top topNodes GO terms identified by the classic algorithm,
#BP
FisherResBP <- GenTable(myBPGOdata, classicFisher = resultFisherBP, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 25)
FisherResBP
myterms = c("GO:0007165", "GO:0044281", "GO:0050789")
mygenes <- genesInTerm(myBPGOdata, myterms)
for (i in 1:length(myterms))
{
  myterm <- myterms[i]
  mygenesforterm <- mygenes[myterm][[1]]
  mygenesforterm <- paste(mygenesforterm, collapse=',')
  print(paste("Term",myterm,"genes:",mygenesforterm))
}
write.csv(FisherResBP, file = "FisherResBP.csv", quote = F, row.names = F)
#MF
FisherResMF <- GenTable(myMFGOdata, classicFisher = resultFisherMF, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 22)
write.csv(FisherResMF, file = "FisherResMF.csv", quote = F, row.names = F)
#CC
FisherResCC <- GenTable(myCCGOdata, classicFisher = resultFisherCC, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 22)
write.csv(FisherResCC, file = "FisherResCC.csv", quote = F, row.names = F)

# 
# resultweight01BP <- runTest(myBPGOdata, algorithm="weight01", statistic="fisher") 
# resultweight01BP
# resultweight01MF <- runTest(myMFGOdata, algorithm="weight01", statistic="fisher") 
# resultweight01MF
# resultweight01CC <- runTest(myCCGOdata, algorithm="weight01", statistic="fisher") 
# resultweight01CC
# 
# #BP
# Fisher_wt_BP <- GenTable(myBPGOdata, topgoFisher = resultweight01BP, orderBy = "topgoFisher", ranksOf = "topgoFisher", topNodes = 22)
# score(Fisher_wt_BP)
# write.csv(Fisher_wt_BP, file = "Fisher_wt_BP.csv")
# #MF
# Fisher_wt_MF <- GenTable(myMFGOdata, topgoFisher = resultweight01MF, orderBy = "topgoFisher", ranksOf = "topgoFisher", topNodes = 22)
# write.csv(Fisher_wt_MF, file = "Fisher_wt_MF.csv")
# #CC
# Fisher_wt_CC <- GenTable(myCCGOdata, topgoFisher = resultweight01CC, orderBy = "topgoFisher", ranksOf = "topgoFisher", topNodes = 22)
# write.csv(Fisher_wt_CC, file = "Fisher_wt_CC.csv")



# resultElimBP <- runTest(myBPGOdata, algorithm="elim", statistic="fisher") 
# resultElimBP
# resultElimMF <- runTest(myMFGOdata, algorithm="elim", statistic="fisher") 
# resultElimMF
# resultElimCC <- runTest(myCCGOdata, algorithm="elim", statistic="fisher") 
# resultElimCC
# 
# result_PC_BP <- runTest(myBPGOdata, algorithm="parentchild", statistic="fisher") 
# result_PC_BP
# result_PC_MF <- runTest(myMFGOdata, algorithm="parentchild", statistic="fisher") 
# result_PC_MF
# result_PC_CC <- runTest(myCCGOdata, algorithm="parentchild", statistic="fisher") 
# result_PC_CC
# 
# # see how many results we get where weight01 gives a P-value <= 0.001:
# mysummary_BP <- summary(attributes(resultweight01BP)$score <= 0.001)
# numsignif_BP <- as.integer(mysummary_BP[[3]]) # how many terms is it true that P <= 0.001
# 
# mysummary_MF <- summary(attributes(resultweight01MF)$score <= 0.001)
# numsignif_MF <- as.integer(mysummary_MF[[3]]) # how many terms is it true that P <= 0.001
# 
# mysummary_CC <- summary(attributes(resultweight01CC)$score <= 0.001)
# numsignif_CC <- as.integer(mysummary_CC[[3]]) # how many terms is it true that P <= 0.001
# 
# 
# # print out the top 'numsignif' results:
# allRes_BP <- GenTable(myBPGOdata, classicFisher = resultFisherBP, elimFisher = resultElimBP, 
#                       topgoFisher = resultweight01BP, parentchildFisher = result_PC_BP, 
#                       orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif_BP)
# allRes_BP
# 
# allRes_MF <- GenTable(myMFGOdata, classicFisher = resultFisherMF, elimFisher = resultElimMF, 
#                       topgoFisher = resultweight01MF, parentchildFisher = result_PC_MF, 
#                       orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif_MF)
# allRes_MF
# 
# allRes_CC <- GenTable(myCCGOdata, classicFisher = resultFisherCC, elimFisher = resultElimCC, 
#                       topgoFisher = resultweight01CC, parentchildFisher = result_PC_CC, 
#                       orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif_CC)
# allRes_CC
