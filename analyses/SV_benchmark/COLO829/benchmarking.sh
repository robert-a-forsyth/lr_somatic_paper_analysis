cd /staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/COLO829
conda activate minda
RAW_COLO_ONT=lr_somatic_paper_analysis/analyses/SV_benchmark/data/COLO829-ONT/severus_somatic.vcf.gz
RAW_COLO_ONT_TO=lr_somatic_paper_analysis/analyses/SV_benchmark/data/COLO829-ONT-TO/severus_somatic.vcf.gz
RAW_COLO_PB=lr_somatic_paper_analysis/analyses/SV_benchmark/data/COLO829-PB/severus_somatic.vcf.gz
RAW_COLO_PB_TO=lr_somatic_paper_analysis/analyses/SV_benchmark/data/COLO829-PB-TO/severus_somatic.vcf.gz

TRUTH_SET=lr_somatic_paper_analysis/analyses/SV_benchmark/COLO829/truthset_somaticSVs_COLO829_hg38lifted_chr_alt_vaf.vcf.gz

### ADD VAF TO RAW FILES
VAF_SCRIPT_PATH=lr_somatic_paper_analysis/analyses/SV_benchmark/add_vaf.py

VAF_SEVERUS_COMP=lr_somatic_paper_analysis/analyses/SV_benchmark/COLO829/vaf_vcfs/severus_pb.vcf.gz
VAF_COLO_ONT=lr_somatic_paper_analysis/analyses/SV_benchmark/COLO829/vaf_vcfs/ONT_vaf.vcf.gz
VAF_COLO_ONT_TO=lr_somatic_paper_analysis/analyses/SV_benchmark/COLO829/vaf_vcfs/ONT_TO_vaf.vcf.gz
VAF_COLO_PB=lr_somatic_paper_analysis/analyses/SV_benchmark/COLO829/vaf_vcfs/PB_vaf.vcf.gz
VAF_COLO_PB_TO=lr_somatic_paper_analysis/analyses/SV_benchmark/COLO829/vaf_vcfs/PB_TO_vaf.vcf.gz

python3 $VAF_SCRIPT_PATH $RAW_COLO_ONT $VAF_COLO_ONT
python3 $VAF_SCRIPT_PATH $RAW_COLO_ONT_TO $VAF_COLO_ONT_TO
python3 $VAF_SCRIPT_PATH $RAW_COLO_PB $VAF_COLO_PB
python3 $VAF_SCRIPT_PATH $RAW_COLO_PB_TO $VAF_COLO_PB_TO

### RUN BENCHMARK
$MINDA_PATH truthset --base $TRUTH_SET --vcfs $VAF_SEVERUS_COMP $VAF_COLO_ONT $VAF_COLO_ONT_TO $VAF_COLO_PB $VAF_COLO_PB_TO --min_size 50 --tolerance 500 --out_dir minda_out --vaf 0.1
