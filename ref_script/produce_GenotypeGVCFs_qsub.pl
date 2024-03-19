use strict;
open OUT,">","GenotypeGVCFs_chr$ARGV[0].qsub" or die "Can't open OUT:$!\n";
print OUT "#!/bin/bash\n";
print OUT "#PBS -N GenotypeGVCFs_chr$ARGV[0]\n";
print OUT "#PBS -l nodes=1:ppn=1\n";
print OUT "#PBS -j oe\n";
print OUT "#Run your executable\n";
print OUT "date\n";
my $var;
if($ARGV[0]<=9){
  $var="0".$ARGV[0];
}else{
  $var=$ARGV[0];
}
my $path = $ARGV[1];
my $name = $ARGV[2];
print OUT "/data5/home/zlgu/miniconda3/bin/gatk GenotypeGVCFs -R /data5/home/zlgu/important_basic_data/IRGSP_v1/IRGSP-1.0_genome.fasta -O ${path}${name}_chr${var}.vcf -V gendb://${path}${name}_chr${var}_db -L chr${var} --only-output-calls-starting-in-intervals true 2>${path}chr${var}_nohup.out\n";
print OUT "date\n";
