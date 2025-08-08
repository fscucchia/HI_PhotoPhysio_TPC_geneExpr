#!/bin/bash
#SBATCH --job-name=mapping_qc
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=200G  # Requested Memory
#SBATCH -t 10:00:00
#SBATCH -o slurm-fastqc_raw.out  # %j = job ID
#SBATCH -e slurm-fastqc_raw.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications

#load modules 
module load uri/main
module load fastqc/0.12.1
module load MultiQC/1.12-foss-2021b

#generate multiqc report Mcap
#multiqc /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Mcap --filename multiqc_report_mapping_Mcap.html -o /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Mcap/multiqc

#generate multiqc report Pacu
#multiqc /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Pacu --filename multiqc_report_mapping_Pacu.html -o /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Pacu/multiqc

#generate multiqc report Pcom
#multiqc /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Pcom --filename multiqc_report_mapping_Pcom.html -o /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Pcom/multiqc

#generate multiqc report Pcom-Cladocopium
#multiqc /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Cladocopium_Pcom --filename multiqc_report_mapping_PcomCladoGTF.html -o /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Cladocopium_Pcom/multiqc

#generate multiqc report Pacu-Cladocopium
#multiqc /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Cladocopium_Pacu --filename multiqc_report_mapping_PacuCladoGTF.html -o /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Cladocopium_Pacu/multiqc

#generate multiqc report Mcap-Clado-Durusd
multiqc /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Cladoc_Durusd_Mcap --filename multiqc_report_mapping_McapCladoDurusdGTF.html -o /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Cladoc_Durusd_Mcap/multiqc

#generate multiqc report Mcap-align-2pass
#multiqc /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Mcap --filename multiqc_report_mapping_Mcap2pass.html -o /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Mcap/multiqc

#generate multiqc report Mcap-align-2pass
#multiqc /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Pcom --filename multiqc_report_mapping_Pcom2pass.html -o /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Pcom/multiqc

#generate multiqc report Mcap-align-2pass
#multiqc /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Pacu --filename multiqc_report_mapping_Pacu2pass.html -o /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Pacu/multiqc

echo "QC of cleaned seq data complete." $(date)
