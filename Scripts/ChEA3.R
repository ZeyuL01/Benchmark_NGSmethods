library(httr)
library(jsonlite)
library(data.table)

####Gene set should be named as TF_id.txt

work_dir<-"path/to/input/directory"
work_files<-list.files(work_dir)

export_dir<-"path/to/output/directory"

file_names<-sapply(strsplit(work_files,".",fixed=TRUE),function(x){return(x[[1]])})
tf_symbol<-sapply(strsplit(file_names,"_",fixed=TRUE),function(x){return(x[[1]])})
tf_id<-sapply(strsplit(file_names,"_",fixed=TRUE),function(x){return(x[[2]])})

for(i in 1:length(work_files)){
	gene_list<-fread(paste0(work_dir,work_files[i]),header=FALSE)[[1]]
	url = "https://maayanlab.cloud/chea3/api/enrich/"
	encode = "json"
	payload = list(query_name = "myQuery", gene_set = gene_list)
	response = POST(url = url, body = payload, encode = encode)
	json = content(response, "text")
	results = fromJSON(json)
	tb_names<-names(results)
	for(j in tb_names){
		file_names<-paste0(export_dir,tf_symbol[i],".",j,".",tf_id[i],".txt")
		fwrite(results[[j]],file_names)
	}
}	
