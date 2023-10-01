import json
import requests
import os

work_dir='path/to/input/directory'
work_files=os.listdir(work_dir)

export_dir='path/to/output/directory'

tf_names=[k.split("_")[0] for k in work_files]

tf_libraries=['ARCHS4_TFs_Coexp','ChEA_2022','ENCODE_TF_ChIP-seq_2015','TRRUST_Transcription_Factors_2019']

ENRICHR_URL_addlist = 'https://maayanlab.cloud/Enrichr/addList'
ENRICHR_URL_export = 'https://maayanlab.cloud/Enrichr/export'

for i in range(len(work_files)):
    with open(work_dir+work_files[i]) as file:
        lines = [line.rstrip() for line in file]
    genes_str = '\n'.join(lines)
    description = 'test gene list'
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }
    response = requests.post(ENRICHR_URL_addlist, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')
    request_data = json.loads(response.text)
    user_list_id = request_data['userListId']
    for j in tf_libraries:
        ENRICHR_URL_export = 'https://maayanlab.cloud/Enrichr/export'
        query_string = '?userListId=%s&filename=%s&backgroundType=%s'
        filename = tf_names[i]+'.'+j
        gene_set_library = j
        url = ENRICHR_URL_export + query_string % (user_list_id, filename, gene_set_library)
        response = requests.get(url, stream=True)
        with open(export_dir+filename +"."+str(i)+'.txt', 'wb') as f:
            for chunk in response.iter_content(chunk_size=1024): 
                if chunk:
                    f.write(chunk)

