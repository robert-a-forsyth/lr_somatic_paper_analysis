#!/usr/bin/env bash
#SBATCH -J lrsomatic_run
#SBATCH --nodes=1
#SBATCH --cores=1
#SBATCH --time=120:00:00
#SBATCH --mem=4G
#SBATCH --cluster=clustery_the_high_performance_cluster
#SBATCH --partition=my_favorite_part_of_the_cluster
#SBATCH --account=cpu_account
#SBATCH --error=%x.e1%a
#SBATCH --output=%x.log1%a


MASTER_DIR=CCS15_LRSomatic_Testing
CUSTOM_VEP=data/annotation/cosmic_tier_rank.sorted.vcf.gz
CUSTOM_VEP_TBI=data/annotation/cosmic_tier_rank.sorted.vcf.gz.tbi

nextflow -log $MASTER_DIR/gr38_CCS15/data/lr_somatic/out/validation.log run IntGenomicsLab/lr_somatic \
 -r 1e737a3 \
 -c /staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/custom.config \
 --fasta /staging/leuven/stg_00096/references/GRCh38.alt-masked-V2-noALT/fasta/Homo_sapiens_assembly38_masked_noALT.fasta \
 -profile vsc_kul_uhasselt,wice \
 --outdir /$MASTER_DIR/gr38_CCS15/data/lr_somatic/out \
 -w $MASTER_DIR/gr38_CCS15/data/work \
 --input /staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/CC15.csv \
 --vep_cache_version 113 \
 --vep_args "--everything --filter_common --per_gene --total_length --offline --format vcf \
 --custom file=$CUSTOM_VEP,short_name=COSMIC_CENSUS,format=vcf,type=overlap,overlap_cutoff=80,fields=MUTATION_SIGNIFICANCE_TIER%RANK%DESCRIPTION" \
 --vep_custom $CUSTOM_VEP \
 --vep_custom_tbi $CUSTOM_VEP_TBI \
 --genome GRCh38 \
 --ascat_pdf_plots "TRUE" \
 --skip_wakhan \
 -resume
# -dump-hashes