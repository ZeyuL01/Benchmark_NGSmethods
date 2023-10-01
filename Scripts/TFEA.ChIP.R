require(TFEA.ChIP)
require(GSA)
require(data.table)

work_dir<-"path/to/input/directory"
work_files<-list.files(work_dir)

output_dir<-"path/to/output/directory"

tf_symbol<-sapply(strsplit(work_files,"_",fixed=TRUE),function(x){return(x[[1]])})

for(i in 1:length(work_files)){
	genes<-fread(paste0(work_dir,work_files[i]),header=FALSE)[[1]]
	gene_ids<-GeneID2entrez(genes)
	gene_cmt<-contingency_matrix(gene_ids)
	tb_cmt<-getCMstats(gene_cmt)
	TF_ranking <- rankTFs(tb_cmt,rankMethod = "gsea")
	fwrite(tb_cmt,paste0(output_dir,"cm_pval/",tf_symbol[i],"_",i,".txt"))
	fwrite(TF_ranking,paste0(output_dir,"tfrank/",tf_symbol[i],"_",i,".txt"))
}

