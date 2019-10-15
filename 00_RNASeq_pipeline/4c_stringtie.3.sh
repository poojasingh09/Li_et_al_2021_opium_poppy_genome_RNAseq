for i in `ls *.bam | cut -f -9 -d"."`;

do

#mkdir $i

/data/programs/stringtie-1.3.5.Linux_x86_64/stringtie $i".bam" -M 0.05 -e -B -G poppy_V3.2_fix_inversion_top_level.mrna.filt.MIA.final.109.strandchecked.gff3 -A ./$i/gene_abundances.tsv -o ./$i/transcripts.gtf -p 2 --rf & 

done
