for i in `ls *.sorted.bam`;

do

/data/programs/stringtie-1.3.5.Linux_x86_64/stringtie $i -G /data/home/psingh/poppy_RNAseq/genome_assembly/annotation_qiushi_v3.2/poppy_v3.2/poppy_V3.2_fix_inversion_top_level.gff3 -o $i".gtf" -p 2 --rf & 

done
