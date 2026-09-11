awk -F'\t' '
    $3=="gene" && ($1=="chr05" || $1=="05") && $4 <= 52800000 && $5 >= 51800000 {
        split($9, arr, ";");
        for (i in arr) {
            if (arr[i] ~ /^ID=/) {
                sub(/^ID=/, "", arr[i]);
                print arr[i]
            }
        }
    }
' DM8.1_gene.gff3 > /public/home/zhaoshuo/work1/01_work/06_candidate/GABA/03_candidate/C05_51.8_52.8_gene.txt

awk '
    NR==FNR {
        if (NR > 1) expr[$1] = 1
        next
    }
    $1 in expr {
        print $1
    }
' zz88_TPM_filtered.txt /public/home/zhaoshuo/work1/01_work/06_candidate/GABA/03_candidate/C05_51.8_52.8_gene.txt > /public/home/zhaoshuo/work1/01_work/06_candidate/GABA/03_candidate/C05_51.8_52.8_candidate_gene.txt

awk -F'\t' '
    NR==FNR {
        genes[$1] = 1
        next
    }
    $3=="gene" {
        split($9, arr, ";")
        id = ""
        for (i in arr) {
            if (arr[i] ~ /^ID=/) {
                id = arr[i]
                sub(/^ID=/, "", id)
            }
        }
        if (id in genes) {
            chrom = $1
            start = $4 - 1
            end = $5
            strand = $7
            print chrom "\t" start "\t" end "\t" id "\t.\t" strand
        }
    }
' /public/home/zhaoshuo/work1/01_work/06_candidate/GABA/03_candidate/C05_51.8_52.8_candidate_gene.txt DM8.1_gene.gff3 > /public/home/zhaoshuo/work1/01_work/06_candidate/GABA/03_candidate/C05_51.8_52.8_candidate_gene.bed

awk '
    NR==FNR {
        genes[$1] = 1
        next
    }
    NR==1 {
        print
        next
    }
    $1 in genes {
        print
    }
' /public/home/zhaoshuo/work1/01_work/06_candidate/GABA/03_candidate/C05_51.8_52.8_candidate_gene.txt zz88_TPM_filtered.txt > /public/home/zhaoshuo/work1/01_work/06_candidate/GABA/03_candidate/C05_51.8_52.8_candidate_expression.txt

awk -F'\t' 'BEGIN {OFS="\t"} {
    start = $2 - 2000; 
    if (start < 0) start = 0; 
    end = $3 + 2000; 
    print $1, start, end
}' /public/home/zhaoshuo/work1/01_work/06_candidate/GABA/03_candidate/C05_51.8_52.8_candidate_gene.bed > /public/home/zhaoshuo/work1/01_work/06_candidate/GABA/03_candidate/C05_51.8_52.8_candidate_genes_slop2k.bed

bcftools view -R /public/home/zhaoshuo/work1/01_work/06_candidate/GABA/03_candidate/C05_51.8_52.8_candidate_genes_slop2k.bed \
              -Oz -o /public/home/zhaoshuo/work1/01_work/06_candidate/GABA/03_candidate/C05_51.8_52.8_candidate_SV.vcf.gz \
              /public/home/zhaoshuo/work1/01_work/05_GWAS/01_SV_GWAS/02_zz/01_PPG1.4_VCF/zz_PPG1.4_SV_sorted_modified.vcf.gz
bcftools index -t /public/home/zhaoshuo/work1/01_work/06_candidate/GABA/03_candidate/C05_51.8_52.8_candidate_SV.vcf.gz

vcftools --gzvcf C05_51.8_52.8_candidate_SV.vcf.gz --maf 0.05 --recode --recode-INFO-all --stdout | bgzip -c > C05_51.8_52.8_candidate_SV_maf0.05.vcf.gz
bcftools index -t C05_51.8_52.8_candidate_SV_maf0.05.vcf.gz

bcftools view C05_51.8_52.8_candidate_SV_maf0.05.vcf.gz | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$3="sv" ++i; print}' | bgzip -c > C05_51.8_52.8_candidate_SV_maf0.05.vcf.gz.tmp && mv C05_51.8_52.8_candidate_SV_maf0.05.vcf.gz.tmp C05_51.8_52.8_candidate_SV_maf0.05.vcf.gz && bcftools index -t C05_51.8_52.8_candidate_SV_maf0.05.vcf.gz

# 1. 将 SV VCF 转换为临时 BED 文件（处理 0-based 坐标）
bcftools query -f '%CHROM\t%POS\t%POS\t%ID\n' C05_51.8_52.8_candidate_SV_maf0.05.vcf.gz | awk 'BEGIN {OFS="\t"} {print $1, $2-1, $2, $4}' > /tmp/sv_temp.bed

# 2. 使用 bedtools window 寻找上下游 2000bp 内的对应关系，并输出 Gene 与 SV ID
bedtools window -a C05_51.8_52.8_candidate_gene.bed -b /tmp/sv_temp.bed -w 2000 | awk 'BEGIN {OFS="\t"; print "Gene_ID", "SV_ID"} {print $4, $10}' > C05_51.8_52.8_gene_sv_mapping.txt

# 3. 清理临时文件
rm /tmp/sv_temp.bed