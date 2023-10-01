import requests
import os
import json
import re

work_dir="path/to/input/directory"
work_files=os.listdir(work_dir)

tf_names=[k.split("_")[0] for k in work_files]
file_names=[k.split(".")[0] for k in work_files]

export_dir="path/to/output/directory"

cscan_url='http://159.149.160.88/pscan/'

for i in range(len(work_files)):
	with open(work_dir+work_files[i]) as file:
		lines = [line.rstrip() for line in file]
	gene_id='\n'.join(lines)
	form_data={'id_text':gene_id,'custom_matrices':'','organism':'Homo_sapiens','region':'-450 +50','matrix':'Jaspar_2020_NR','vc':'0'}
	server=requests.post(cscan_url,data=form_data)
	search_str=re.search(r'\./Jaspar_2020_NR/.*\.id\.res',server.text)
	down_link='http://159.149.160.88/pscan'+search_str.group()[1:]
	down=requests.get(down_link,allow_redirects=True)
	open(export_dir+file_names[i]+".txt", 'wb').write(down.content)


