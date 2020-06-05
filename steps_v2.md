#######################################
## Commands for population genetics ###
#######################################
# 0- After mapping to host, need to sort bam files before stats
module load samtools/1.9

samtools sort ${SAMPLE}_map2host.bam > ${SAMPLE}_sorted.bam

# 1- BAM stats
# Sen utiliza Qualimap, que da mucha información sobre el mapeo, pero es difícil de visualizar todas las muestras en conjunto y parece tener ciertos problemas al numerar duplicados. Se puede usar Qualimap, o extraer los datos y luego visualizarlos en R. Añado las dos versiones.

#1.A- Qualimap - lo copio tal cual lo hace Sen
rule bam_qc:
        input:
                bam_all="{WorkDir}/03_MapToRef/{sample}_all.bam",
                bam_map2host="{WorkDir}/03_MapToRef/{sample}_map2host.bam",
                bam_metaG="{WorkDir}/03_MapToRef/{sample}_metaG.bam"
        output:
                bam_metaG="{WorkDir}/04_BamQC/{sample}_metaG.bam"
        threads: 4
        params:
                prefix="{WorkDir}/04_BamQC/{sample}"
        shell:
                """
                module load intel/perflibs/2019_update5 gcc/8.2.0 java/1.8.0-openjdk R/3.6.1 qualimap/2.2.1
                qualimap bamqc -bam {input.bam_all} -outdir {params.prefix}_all -outformat html --java-mem-size=16G
                qualimap bamqc -bam {input.bam_map2host} -outdir {params.prefix}_map2host -outformat html --java-mem-size=16G
                qualimap bamqc -bam {input.bam_metaG} -outdir {params.prefix}_metaG -outformat html --java-mem-size=16G
                mv {input.bam_metaG} {output.bam_metaG}
                """

#1.A- Stats Sofi - Son distintas posibilidades
module load samtools/1.9
module load java/1.8.0 picard-tools/2.9.1
module load anaconda3/4.4.0 gatk/4.1.7.0
module load bedtools/2.28.0

  # 1.1- Mapping stats
  # yo creo que lo más importante son los siguientes parámetros:
    # - Depth of coverage: el número de veces una base del genoma es cubierto. Esto técnicamente es llamado per-base depth. Overall depth sería la media de profundidad del alineamiento en todo el genoma.
    # - Breadth of coverage: el porcentage de bases del genoma de referencia que han sido cubiertos con una profundidad predeterminanda.
    # - Duplicates: la cantidad de reads duplicados después del alineamineto.
    # - Insert size: la longitud de la secuencia entre los adaptadores.
    # - Percentage mapped: el porcentaje de reads que ha sido alineado con el genoma de referencia.

  # 1.1- Alignment summary
  samtools flagstat ${SAMPLE}_map2host.bam > ${SAMPLE}.flag
  # Da un brief summary del mapeo: el porcentaje de mapped reads y la cantidad, cantidad de duplicados, la cantidad/porcentage de mate mapped y singletons.

  # 1.2- Insert size --> la media de insert size y la desviación estándar. Se pueden plotear histogramas etc en R
  picard CollectInsertSizeMetrics -I ${SAMPLE}_map2host.bam -O ${SAMPLE}.insert H=insert_size_histogram.pdf M=0.5

  # 1.3- Duplicates --> da duplicate ratio, y distintas estadísticas. Se pueden plotear histogramas etc en R
  picard MarkDuplicates -I ${SAMPLE}_map2host.bam -O ${SAMPLE}_marked_sup.bam -M ${SAMPLE}_mark.txt

  # 1.4- Depth
  gatk DepthofCoverage #todavía no consigo hacerlo funcionar. Me da error

  # 1.5- Breadth/Depth con samtools --> con el flag -a le pides que enumere los bases cubiertos y no cubiertos. De aquí se pueden sacar estos datos
  samtools depth -a ${SAMPLE}_map2host.bam > ${SAMPLE}.depth
  less ${SAMPLE}.depth | awk '{c++;s+=$3}END{print s/c}' # overall depth
  less ${SAMPLE}.depth | awk -F "[(|%]" 'NR == 3 {print $2}' # overall breadth

  # 1.6- Breadth --> también como samtools depth
  bedtools genomecov -ibam .bam -d g ${REF}.fai -max 20 > ${SAMPLE}.dcov 2> >(tee "$logfile") # con -d registra las bases cubiertas y no cubiertas
  bedtools genomecov -ibam .bam -bg g ${REF}.fai -max 20 > ${SAMPLE}.dcov 2> >(tee "$logfile") # solo registra bases cubiertas

# 2- Calling variants and filtering
module load bcftools/1.9

find *_sorted.bam > samples_list.txt # tiene que ser un file con los nombres de las muestras, pero el nombre que tienen en el SAM header
chromosomes.txt # tienen que aparecer los nombres de los cromosomas como aparecen en NCBI, pero creo que es más recomendable hacer el variant cal para cada cromosoma aparte, para después estudiar las estadísticas. Yo hago separado.
bcftools mpileup -C 50 -q 30 -Q 20 -Ou -f $REF --regions-file chromosomes.txt --samples-file samples_list.txt | bcftools call -Ou -m -v | bcftools filter -s LowQual -e '%QUAL<30 || DP<(AVG(DP)*3)'-Oz -o call_all_${CHR}.vcf.gz
gatk IndexFeatureFile -I call_all_${CHR}.vcf.gz

# 3- Variant stats
module load perl/5.30.2 vcftools/0.1.8

# allele freq
vcftools --gzvcf $SUBSET_VCF --min-alleles 2 --max-alleles 2 --freq2 --out $OUT

# depth per ind
vcftools --gzvcf $SUBSET_VCF --min-alleles 2 --max-alleles 2 --depth --out $OUT

# depth per site
vcftools --gzvcf $SUBSET_VCF --min-alleles 2 --max-alleles 2 --site-mean-depth --out $OUT

# het per ind
vcftools --gzvcf $SUBSET_VCF --min-alleles 2 --max-alleles 2 --het --out $OUT

# site quality
vcftools --gzvcf $SUBSET_VCF --min-alleles 2 --max-alleles 2 --site-quality --out $OUT

# missing ind
vcftools --gzvcf $SUBSET_VCF --min-alleles 2 --max-alleles 2 --missing-indv --out $OUT

# missing site
vcftools --gzvcf $SUBSET_VCF --min-alleles 2 --max-alleles 2 --missing-site --out $OUT

# site quality
vcftools --gzvcf $SUBSET_VCF --min-alleles 2 --max-alleles 2 --site-quality --out $OUT


# 4- Improve filtering steps if needed

# 5- Merge all VCF files
module load bcftools/1.9

find *${CHR}.vcf.gz > chrs_list.txt
bcftools concat -f chrs_list.txt -Oz -o all_variants.vcf.gz

# 6- Select SNPs
bcftools view -m2 -M2 -v snps -Oz -o chickens.vcf.gz raw_SNPs_only.vcf.gz

# 6.1- Separate subpopulations
# subir dos text files con los nombres de los individuos de cada población
bcftools view -S Cobbs.txt -Oz -o ${RAW}_SNPs_cobb.vcf.gz ${RAW}_SNPs_only.vcf.gz
bcftools view -S Rosses.txt -Oz -o ${RAW}_SNPs_ross.vcf.gz ${RAW}_SNPs_only.vcf.gz

# 7- Filtering SNPs

# 8- Annotate SNPs
module load perl/5.30.2 snpeff/4.3r
snpeff -c /home/people/sofbas/CustomSoftware/snpEff/snpEff.config -v GRCg6a all_raw_SNPs.vcf.gz > all_raw_SNPs_ann_snpEff.vcf
bgzip all_raw_SNPs_ann_snpEff.vcf
gatk IndexFeatureFile -I all_raw_SNPs_ann_snpEff.vcf.gz
