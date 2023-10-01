import os
import re
from bs4 import BeautifulSoup

work_dir="path/to/input/directory"
work_files=os.listdir(work_dir)

export_dir="path/to/output/directory"

TF_names=[x.split("_")[0] for x in work_files]
rank_lists=[]
for j in range(len(work_files)):
	print(j)
	val="NA"
	with open(work_dir+work_files[j]+"/report.html","r") as file:
		data=file.read()
	html_file=BeautifulSoup(data,'html.parser')
	strhtml=html_file.prettify()
	test=re.findall("Description:.*\n",strhtml)
	test=[x.replace("\n","") for x in test]
	TF_rank=[x.replace(u'\xa0',u' ').split(' ')[1] for x in test]
	TF_rank2=[x.replace(u'\xa0',u' ').split(' ')[-1] for x in test]
	for i in range(len(TF_rank)):
		if(TF_rank[i]=="ChIP-seq"):
			TF_rank[i]=TF_rank2[i]
	for i in range(len(TF_rank)):
		if TF_rank[i]==TF_names[j]:
			val=i+1
	rank_lists.append(str(val))

with open(export_dir+"icistarget_rank_200.txt","w") as f:
	for item in rank_lists:
		f.write("%s\n" % item)
	print("Done")




