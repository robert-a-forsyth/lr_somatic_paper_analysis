
COLO829_PB_PAIR_SEVERUS="/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/data/lr_somatic/out/COLO829-PB/variants/severus/somatic_SVs/severus_somatic.vcf.gz"

bcftools filter -i 'FILTER="PASS"' $COLO829_PB_PAIR_SEVERUS \
| bcftools filter -i 'VAF>=0.1' | bcftools view -s COLO829-PB_tumor -o COLO829-PB_severus.filtered.vcf