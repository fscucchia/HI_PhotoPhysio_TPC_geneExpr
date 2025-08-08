#!/bin/bash
#SBATCH --job-name=concatenate
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=200G  # Requested Memory
#SBATCH -t 12:00:00
#SBATCH -o slurm-acr-concatenate.out  # %j = job ID
#SBATCH -e slurm-acr-concatenate.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Symbiont_Genomes_concat/Cladocopium_Durusdinium

### Genomes

# C_goreaui_SCF055='/work/pi_hputnam_uri_edu/Symbiont_Genomes/Cladocopium_goreaui_SCF055/Cladocopium_goreaui/Cladocopium_goreaui.genome.fa'
# Cladocopium_C15='/work/pi_hputnam_uri_edu/Symbiont_Genomes/Cladocopium_sp_C15/SymbC15_plutea_v2.1.fna'
# Cladocopium_C92='/work/pi_hputnam_uri_edu/Symbiont_Genomes/Cladocopium_sp_C92/Cladocopium_sp_C92/Cladocopium_sp_C92.genome.fa'
# C_goreaui_C1='/work/pi_hputnam_uri_edu/Symbiont_Genomes/C_goreaui_cladeC1/SymbC1.Genome.Scaffolds.fasta'
# Symbiodinium_CladeC='/work/pi_hputnam_uri_edu/Symbiont_Genomes/Symbiodinium_CladeC/symC_scaffold_40.fasta'
# Durusdinium_sp='/work/pi_hputnam_uri_edu/Symbiont_Genomes/Durusdinium_sp/102_symbd_genome_scaffold.fa' 
# Durusdinium_trenchii_CCMP2556='/work/pi_hputnam_uri_edu/Symbiont_Genomes/Durusdinium_trenchii_CCMP2556_SCF082/Dtrenchii_CCMP2556_ASSEMBLY_fasta' 
# Durusdinium_trenchii_SCF082='/work/pi_hputnam_uri_edu/Symbiont_Genomes/Durusdinium_trenchii_CCMP2556_SCF082/Dtrenchii_SCF082_ASSEMBLY_fasta' 

#### concatenate Cladocopium sp genomes
# cd /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Symbiont_Genomes_concat/Cladocopium

# cat "$C_goreaui_SCF055" "$Cladocopium_C15" "$Cladocopium_C92" "$C_goreaui_C1" "$Symbiodinium_CladeC" > Cladocopium_concat.fasta

# #### concatenate Durusdinium sp genomes
# cd /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Symbiont_Genomes_concat/Durusdinium

# cat "$Durusdinium_sp" "$Durusdinium_trenchii_CCMP2556" "$Durusdinium_trenchii_SCF082" > Durusdinium_concat.fasta

#### concatenate Cladocopium + Durusdinium sp genomes
# cd /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Symbiont_Genomes_concat/Cladocopium_Durusdinium

# cat "$C_goreaui_SCF055" "$Cladocopium_C15" "$Cladocopium_C92" "$C_goreaui_C1" "$Symbiodinium_CladeC" "$Durusdinium_sp" "$Durusdinium_trenchii_CCMP2556" "$Durusdinium_trenchii_SCF082" > Cladocopium_Durusdinium_concat.fasta


### Gff

C_goreaui_SCF055='/work/pi_hputnam_uri_edu/Symbiont_Genomes/Cladocopium_goreaui_SCF055/Cladocopium_goreaui/Cladocopium_goreaui.gff'
Cladocopium_C15='/work/pi_hputnam_uri_edu/Symbiont_Genomes/Cladocopium_sp_C15/SymbC15_plutea_v2.1.fna.evm.final.gff3'
Cladocopium_C92='/work/pi_hputnam_uri_edu/Symbiont_Genomes/Cladocopium_sp_C92/Cladocopium_sp_C92/Cladocopium_sp_C92.gff'
C_goreaui_C1='/work/pi_hputnam_uri_edu/Symbiont_Genomes/C_goreaui_cladeC1/SymbC1.Gene_Models.GFF3'
Symbiodinium_CladeC='/work/pi_hputnam_uri_edu/Symbiont_Genomes/Symbiodinium_CladeC/40_symb.gff'
Durusdinium_sp='/work/pi_hputnam_uri_edu/Symbiont_Genomes/Durusdinium_sp/102_symbd.gff' 
Durusdinium_trenchii_CCMP2556='/work/pi_hputnam_uri_edu/Symbiont_Genomes/Durusdinium_trenchii_CCMP2556_SCF082/Dtrenchii_CCMP2556_ANNOT_gff' 
Durusdinium_trenchii_SCF082='/work/pi_hputnam_uri_edu/Symbiont_Genomes/Durusdinium_trenchii_CCMP2556_SCF082/Dtrenchii_SCF082_ANNOT_gff' 

#### concatenate Cladocopium sp gff files
cd /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Symbiont_Genomes_concat/Cladocopium_Durusdinium

cat "$C_goreaui_SCF055" "$Cladocopium_C15" "$Cladocopium_C92" "$C_goreaui_C1" "$Symbiodinium_CladeC" "$Durusdinium_sp" "$Durusdinium_trenchii_CCMP2556" "$Durusdinium_trenchii_SCF082" > CladocDurusd_concat.gff
