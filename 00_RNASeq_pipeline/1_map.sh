for i in `cat input_files`;

do

/data/programs/hisat2-2.1.0/hisat2 -p 2 --dta --new-summary --summary-file $i"_hisatsumm" --met-file $i"_hisatmetrics" -x /data/home/psingh/poppy_RNAseq/genome_assembly/hisat2_index/PGA_fixinv -1 $i"_1P.fq.gz" -2 $i"_2P.fq.gz" -S $i"_PGA.sam" &

done
