import requests
import os
import json
import re
import time

def retrieve_result(wabi_id_list,name_list,export_dir):
	export_dir=export_dir
	for i in range(len(wabi_id_list)):
		text=requests.get('https://ddbj.nig.ac.jp/wabi/chipatlas/'+wabi_id_list[i]+'?info=result&format=tsv')
		open(export_dir+name_list[i]+".tsv", 'wb').write(text.content)

work_dir="path/to/input/directory"
work_files=os.listdir(work_dir)
export_dir="path/to/output/directory"

tf_names=[k.split(".")[0] for k in work_files]

chipatlas_url='https://ddbj.nig.ac.jp/wabi/chipatlas/'

wabi_id_ls=[]
name_list=[]
for i in range(len(work_files)):
	with open(work_dir+work_files[i]) as file:
		lines = [line.rstrip() for line in file]
	gene_id='\n'.join(lines)
	form_data={'format':'text','result':'www','genome':'hg38','antigenClass':'TFs and others','cellClass':'All cell types','threshold':'50','typeA':'gene','bedAFile':gene_id,'typeB':'refseq','bedBFile':'empty','permTime':1,'title':'test_api','descriptionA':'test_gene','descriptionB':'refseq','distanceUp':'5000','distanceDown':'5000','qsubOptions':'-N test'}
	server=requests.post(chipatlas_url,data=form_data,allow_redirects=True)
	output=server.text
	wabi_id=re.search(r'wabi_chipatlas_.*[0-9]{6}',output).group()
	wabi_id_ls.append(wabi_id)
	name_list.append(work_files[i])
	if (i>0)&(i%40==0 or i==(len(work_files)-1)):
		time.sleep(1800)
		retrieve_result(wabi_id_ls,name_list,export_dir)
		wabi_id_ls=[]
		name_list=[]





