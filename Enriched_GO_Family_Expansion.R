setwd("C:/Users/El3ti/Documents/Trabajo_Ciencia/Doctorado-Platelmintos/Snail_Genome/Notes_Miniorthofinder/CAFE_Expansion/GO_TERM")

library(topGO)

Files <- gsub(".txt", "" ,list.files(path=".", pattern=".txt$"))
go_classes <- c("CC","BP","MF")
pvalue_good <- .05

readMappings("Background_PHO_All.in") -> anotacionGO  # Tabla de anotación $1 GenID; $2 TermGO/nada 
readLines("All_PHO.id") -> totalGenes # Archivo con las IDs de todos los genes (1/linea)

for (spe in Files) {
  Interes_genes <- paste (spe,".txt", sep = '')
  print(spe)
  readLines(Interes_genes) -> genes.interes # Archivo con las IDs de los genes de interes (1/linea)

  geneList <- factor(as.integer(totalGenes %in% genes.interes))
  names(geneList) <- totalGenes
  
  for (go_c in go_classes) {
    sampleGOdata <- new("topGOdata", ontology=go_c, allGenes=geneList, nodeSize = 5, annot = annFUN.gene2GO, gene2GO = anotacionGO)
    resultFisher <- runTest(sampleGOdata, statistic = "fisher", algorithm = "classic")
    
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

