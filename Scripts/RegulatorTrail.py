import json
import os
import time
import pandas as pd
# Load method definitions
from graviton import *
# Load assignment of samples into groups

path="path/to/input/directory"
files=os.listdir(path)
path_test=path+files[1]
export_path='path/to/output/directory'

for i in range(len(files))[1:]:
	file_path=path+files[i]
	key=getSession()
	time.sleep(3)
	job_id=uploadFile(key,file_path)['id']
	regulatorORA(key,job_id,'ora','0')
	result=runJob(key)['regulator_ora']['id']
	downloadResult(key,result,export_path+"temp200.json")
	temp_tf=pd.read_json(export_path+"temp200.json")
	temp_tf.to_csv(export_path+"knockTF_200/"+files[i])
	os.remove(export_path+"temp200.json")
	time.sleep(1)
	print(files[i]+" finished")

