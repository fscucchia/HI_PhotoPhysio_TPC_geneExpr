#!/bin/bash
#SBATCH --job-name=Mcap-2pass
#SBATCH --nodes=1 --ntasks-per-node=48
#SBATCH --mem=200G  # Requested Memory
#SBATCH -t 36:00:00
#SBATCH -o slurm-Mcap-2pass.out  # %j = job ID
#SBATCH -e slurm-Mcap-2pass.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Mcap

#load modules
echo "Loading programs" $(date)
module load uri/main
module load STAR/2.7.11b-GCC-12.3.0

#### define parameters for genome index

#genome='/work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.assembly.fasta'
#star_genome_dir='/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/GenomeIndex_Mcap'
#gtf='/work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.genes_fixed.gff3'
#genomeSAindexNbases=13
#limitSjdbInsertNsj=1200000 # default 1000000

### genome index generation

#STAR --runMode genomeGenerate \
#     --runThreadN 20 \
#     --genomeFastaFiles \"$genome\" \
#     --genomeDir \"$star_genome_dir\" \
#	 --sjdbGTFfile \"$gtf\" \
#     --outFileNamePrefix \"$star_genome_dir/star\" \
#     --genomeSAindexNbases $genomeSAindexNbases \
#	 --limitSjdbInsertNsj $limitSjdbInsertNsj


# ### alignment

#  echo "Starting read alignment." $(date)

#  #loop through all files to align them to genome
#  for i in /work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/output/Cleaned_reads/trim.Mcap*_R1_001.fastq.gz; do

#  # Define the corresponding R2 file by replacing _R1_ with _R2_
#      r2_file="${i/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    
#  # Define output prefix
#      output_prefix="/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Mcap/$(basename "${i%_R1_001.fastq.gz}")_"
    
#   # Run STAR alignment
#      STAR --runMode alignReads \
#          --genomeDir /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/GenomeIndex_Mcap \
#          --runThreadN 10 \
#          --readFilesCommand zcat \
#          --readFilesIn "$i" "$r2_file" \
#          --outSAMtype BAM SortedByCoordinate \
#          --outSAMunmapped Within \
#          --outSAMattributes Standard \
#          --genomeSAindexNbases 13 \
#          --outFileNamePrefix "$output_prefix" \
#          --quantMode GeneCounts
#  done

#  echo "Alignment of Trimmed Seq data complete." $(date)



### alignment - --twopassMode ON

  echo "Starting read alignment." $(date)

 #loop through all files to align them to genome
 for i in /work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/output/Cleaned_reads/trim.Mcap*_R1_001.fastq.gz; do

 # Define the corresponding R2 file by replacing _R1_ with _R2_
     r2_file="${i/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    
 # Define output prefix
     output_prefix="/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Mcap/$(basename "${i%_R1_001.fastq.gz}")_"
    
  # Run STAR alignment
     STAR --runMode alignReads \
         --genomeDir /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/GenomeIndex_Mcap \
         --runThreadN 10 \
         --twopassMode Basic \
         --readFilesCommand zcat \
         --readFilesIn "$i" "$r2_file" \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --outSAMattributes Standard \
         --genomeSAindexNbases 13 \
         --outFileNamePrefix "$output_prefix" \
         --quantMode GeneCounts
 done

 echo "Alignment of Trimmed Seq data complete." $(date)
