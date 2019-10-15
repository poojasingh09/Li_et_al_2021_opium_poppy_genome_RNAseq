for i in `cat bamlist`;

do

samtools view -h -F 256 $i | samtools sort > $i".filt.bam" &

done
