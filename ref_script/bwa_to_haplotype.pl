use strict;
use File::Basename;
##$ARGV[0]=input fasta file path and output path, end up with /
##$ARGV[1]=accession ID
##$ARGV[2]=pair end fastq *_1.fastq
##$ARGV[3]=pair end fastq *_2.fastq
##$ARGV[4]=reference
##bwa alignment
   my $filename;
   print "/data5/tool/jre1.8.0_144/bin/java -jar /data5/home/zlgu/tool/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 -phred33 $ARGV[2] $ARGV[3] -baseout $ARGV[0]$ARGV[1].fastq ILLUMINACLIP:/data5/home/zlgu/tool/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10:2:true MAXINFO:50:0.6 MINLEN:51\n";  
   system "/data5/tool/jre1.8.0_144/bin/java -jar /data5/home/zlgu/tool/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 -phred33 $ARGV[2] $ARGV[3] -baseout $ARGV[0]$ARGV[1].fastq ILLUMINACLIP:/data5/home/zlgu/tool/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10:2:true MAXINFO:50:0.6 MINLEN:51";
   
   system "head -35000000 $ARGV[0]${ARGV[1]}_1P.fastq >$ARGV[0]${ARGV[1]}_5x_1P.fastq";
   system "head -35000000 $ARGV[0]${ARGV[1]}_2P.fastq >$ARGV[0]${ARGV[1]}_5x_2P.fastq";

   system "echo ...Alignment with Bwa mem...";
   system "/data5/tool/bwa-0.7.1/bwa mem $ARGV[4] $ARGV[0]${ARGV[1]}_5x_1P.fastq $ARGV[0]${ARGV[1]}_5x_2P.fastq -t 4 -M -R \"\@RG\tID:${ARGV[1]}ID\tPL:Illumina\tSM:$ARGV[1]\" >$ARGV[0]$ARGV[1].sam";

   system "echo ...sam transformed to bam with samtools view...";
   system "/data6/tool/samtools-1.9/bin/samtools view -F 0x100 -Sb $ARGV[0]$ARGV[1].sam >$ARGV[0]$ARGV[1].bam";

   system "echo ...bam sorted with gatk SortSam...";
   system "/data5/home/zlgu/miniconda3/bin/gatk SortSam --INPUT $ARGV[0]$ARGV[1].bam --OUTPUT $ARGV[0]${ARGV[1]}_sort.bam --SORT_ORDER \"coordinate\"";

   system "echo ...MarkDuplicates...";
   system "/data5/home/zlgu/miniconda3/bin/gatk MarkDuplicates --INPUT $ARGV[0]${ARGV[1]}_sort.bam --OUTPUT $ARGV[0]${ARGV[1]}_sort_MDup.bam --METRICS_FILE $ARGV[0]${ARGV[1]}_sort_MDup_metrix.txt --VALIDATION_STRINGENCY SILENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 --ASSUME_SORT_ORDER \"coordinate\" --CREATE_INDEX true";

   system "echo ...HaplotypeCaller...";
   system "/data5/home/zlgu/miniconda3/bin/gatk HaplotypeCaller -R $ARGV[4] -I $ARGV[0]${ARGV[1]}_sort_MDup.bam -O $ARGV[0]${ARGV[1]}.g.vcf -ERC GVCF";
   $filename = $ARGV[0].${ARGV[1]}.".g.vcf";

   system "echo ...remove middle file...";  
   system "rm $ARGV[0]${ARGV[1]}_sort_MDup.bam $ARGV[0]${ARGV[1]}_sort_MDup_metrix.txt $ARGV[0]${ARGV[1]}_sort_MDup.bai $ARGV[0]${ARGV[1]}_sort.bam $ARGV[0]${ARGV[1]}.bam $ARGV[0]${ARGV[1]}.sam $ARGV[0]${ARGV[1]}_1P.fastq $ARGV[0]${ARGV[1]}_2P.fastq $ARGV[0]${ARGV[1]}_1U.fastq $ARGV[0]${ARGV[1]}_2U.fastq $ARGV[0]${ARGV[1]}_5x_1P.fastq $ARGV[0]${ARGV[1]}_5x_2P.fastq"
     if -s $filename;
#   $filename = $ARGV[0].${ARGV[1]}."_sort_MDup.bam";
#   system "rm $ARGV[0]${ARGV[1]}_sort.bam $ARGV[0]${ARGV[1]}.bam $ARGV[0]${ARGV[1]}.sam"
#      if -s $filename; 
