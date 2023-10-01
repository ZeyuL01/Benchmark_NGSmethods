# Review_NGS_MTR
All files can be accessed from box folder

Link to box folder: https://smu.box.com/s/hlofza4oe08ob4bb4i02nu3wl2zy39zj

#### Description of files in the folder:

knockTF.zip - knockTF TF perturbation experiment derived genes, including the original data sourced from the knockTF database and benchmarking used top 200/600/1000 genes.








### Command and scripts
Command and scripts run on Ubuntu 22.04, SMU M3 high-performance computing (HPC) cluster.

#### Command
$FILE_PATH = "/path/to/destination/geneset.txt"
$OUTPUT_PATH = "/path/to/output/destination/geneset.txt"

BART: bart2 geneset -i $FILE_PATH -s hg38 --outdir $OUTPUT_PATH

