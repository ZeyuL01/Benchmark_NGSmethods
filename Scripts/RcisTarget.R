#BiocManager::install("RcisTarget")

library(RcisTarget)
library(data.table)

motif_database <- "path/to/motif/database"

# Load the motif ranking
motifRankings <- importRankings(paste0(motif_database,"hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"))
# Load the annotation to human transcription factors
data(motifAnnotations_hgnc)

work_dir<-"path/to/input/directory"
work_files<-list.files(work_dir)

output_names<-sapply(strsplit(work_files,".",fixed=TRUE),function(x){return(x[[1]])})
output_dir<-"path/to/output/directory"

for(i in 1:length(work_files)){
	print(i)
	geneset<-fread(paste0(work_dir,work_files[i]),header=FALSE)[[1]]
	geneList<-list(TestSet=geneset)

	motifEnrichmentTable_wGenes <- cisTarget(geneList, 
         motifRankings,
         motifAnnot=motifAnnotations,nes=0)

    anotatedTfs <- lapply(split(motifEnrichmentTable_wGenes$TF_highConf,
                            motifEnrichmentTable_wGenes$geneSet),
                      function(x) {
                        genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
                        genesSplit <- unique(unlist(strsplit(genes, "; ")))
                        return(genesSplit)
                        })

    result_list <- anotatedTfs$TestSet
    result_df<-data.frame(rankings=result_list)
    fwrite(result_df,paste0(output_dir,output_names[i],".csv"))
}



