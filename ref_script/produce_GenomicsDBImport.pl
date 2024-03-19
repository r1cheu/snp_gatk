use strict;
open OUT,">","GenomicsDBImport_newHybrid_mt.qsub" or die "Can't open OUT:$!\n";
print OUT "#!/bin/bash\n";
print OUT "#PBS -N GenomicsDBImport_newHybrid_mt\n";
print OUT "#PBS -l nodes=1:ppn=1\n";
print OUT "#PBS -j oe\n";
print OUT "#Run your executable\n";
print OUT "date\n";
#my $var;
#if($ARGV[0]<=9){
#  $var="0".$ARGV[0];
#}else{
#  $var=$ARGV[0];
#}
print OUT "/data6/home/zlgu/miniconda3/bin/gatk GenomicsDBImport --genomicsdb-workspace-path /data5/home/zlgu/data2/male_sterile_gene/newHybrd_fastq/newhybrid_mt_db --batch-size 100 --reader-threads 5 -L chr${var} ".'$(for i in $(ls /data6/home/zlgu/gvcf_pools/gvcfs/*g.vcf);do echo "-V $i";done)'."\n";
print OUT "date\n";
