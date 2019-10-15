for i in `ls *sam | cut -f -6 -d"."`;

do

#echo $i
samtools view -Sb $i".sam" | samtools sort - >  $i".sorted.bam" &


done
