
#!/bin/bash

#quantify only cdna i assume
path_to_created_index=/media/rna/VAQUITA/snRNAseq_potato/00_META
path_to_data=/media/rna/VAQUITA/snRNAseq_potato/00_RAW_DATA

cd ${path_to_data}


 kb count -x 10xv3\
 -o /media/rna/VAQUITA/snRNAseq_potato/02_KB_count/\
 -i ${path_to_created_index}/index.idx\
 -g ${path_to_created_index}/t2g.txt\
 --batch-barcodes samples.txt\
 --verbose


 exit