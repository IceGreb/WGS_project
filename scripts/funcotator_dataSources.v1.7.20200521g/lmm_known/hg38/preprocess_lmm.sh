set -e

VCF=LMM_Path_LP_VUS5-variants-6-12-18.sorted_liftover_b38.vcf
OUTPUT=LMM_Path_LP_VUS5-variants-6-12-18.sorted_liftover_b38.corrected.vcf
date
echo "First egrep .... "
egrep "^#" ${VCF} | sed -r "s/^##contig=<ID=/##contig=<ID=chr/g" > ${OUTPUT}
date
echo "Second egrep (+ sed) .... "
egrep -v "^#" ${VCF} | sed -r "s/^/chr/g" >>${OUTPUT}
date