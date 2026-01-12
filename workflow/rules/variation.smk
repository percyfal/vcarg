rule bcftools_view_biallelic_variant_sites:
    """Generate biallelic variants only output file"""
    output:
        vcf="<results>/biallelic-{callmode}/{callset}.vcf.gz",
        csi="<results>/biallelic-{callmode}/{callset}.vcf.gz.csi",
    input:
        vcf="<results>/gatk-gather-vcfs-{callmode}/{callset}.vcf.gz",
        tbi="<results>/gatk-gather-vcfs-{callmode}/{callset}.vcf.gz.tbi",
    conda:
        "../envs/variation.yaml"
    benchmark:
        "<benchmarks>/bcftools_view_variant_sites/<results>/biallelic-{callmode}/{callset}.vcf.gz.benchmark.txt"
    log:
        "<logs>/bcftools_view_variant_sites/<results>/biallelic-{callmode}/{callset}.vcf.gz.log",
    threads: 1
    shell:
        """
        bcftools filter -i "QUAL>0" {input.vcf} | bcftools view -v snps -m 2 -M 2 - -o {output.vcf} -O z 2>{log}
        bcftools index {output.vcf}
        """
