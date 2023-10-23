setwd("/home/javier/Documentos/Repositorios/Doctorado-Platelmintos/Snail_Genome/Functional_Annotation")

library(topGO)

Files <- gsub("_GO.tab", "" ,list.files(path=".", pattern="_GO.tab$"))
go_classes <- c("CC","BP","MF")
pvalue_good <- .05

readMappings("Full_Genome.tab") -> anotacionGO  # Tabla de anotación $1 GenID; $2 TermGO/nada 
readLines("Full_Genome.id") -> totalGenes # Archivo con las IDs de todos los genes (1/linea)

for (spe in Files) {
  Interes_genes <- paste (spe,"_GO.tab", sep = '')
  print(spe)
  readLines(Interes_genes) -> genes.interes # Archivo con las IDs de los genes de interes (1/linea)

  geneList <- factor(as.integer(totalGenes %in% genes.interes))
  names(geneList) <- totalGenes
  
  for (go_c in go_classes) {
    sampleGOdata <- new("topGOdata", ontology=go_c, allGenes=geneList, nodeSize = 10, annot = annFUN.gene2GO, gene2GO = anotacionGO)
    resultFisher <- runTest(sampleGOdata, statistic = "fisher", algorithm = "Weight01")
    
    sig_go <- sum(score(resultFisher) <= pvalue_good)
    if (sig_go > 0) {
      GOtermsdeInteres.data.frame <- GenTable(sampleGOdata, P_Value = resultFisher, orderBy = "Weight01Fisher", ranksOf = "Weight01Fisher", topNodes = sig_go)
    } else {
      GOtermsdeInteres.data.frame = data.frame()
    }
    
    print_tab <- print(paste(spe,"_",go_c,".tab", sep = ''))
    print_tab
    write.table(GOtermsdeInteres.data.frame,print_tab, sep="\t")
  }
}

