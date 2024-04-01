// STAR FeatureCounts pipline

params.readsFile = ""
//params.sp = ""
params.genome = ""
params.gtf = ""
params.prefix = ""

//if( params.sp == "human" ) {
//    genome = "/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/"
//    gtf = "/SGRNJ01/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.92.chr.gtf"
//}
//else {
//    genome = "/SGRNJ/Public/Database/genome/Mouse/ensembl_92/"
//    gtf = "/SGRNJ01/Public/Database/genome/Mouse/ensembl_92/Mus_musculus.GRCm38.92.chr.gtf"
//}

prfx = params.prefix
genome = params.genome
gtf = params.gtf

starOutdir = file("$launchDir/out/star_out/$prfx/")
featureCountsOutdir = file("$launchDir/out/featureCounts_out/$prfx/")
starOutdir.mkdirs()
featureCountsOutdir.mkdirs()


process runSTAR {
    beforeScript 'source ~/.profile'
    cpus 4
    input:
        path reads
        path genome
    output:
        val "$launchDir/out/star_out/$prfx/${prfx}_Aligned_Sorted.out.bam"

    """
    STAR --genomeDir $genome \
        --readFilesIn $reads \
        --outFileNamePrefix "$launchDir/out/star_out/$prfx/${prfx}_" \
        --outSAMtype BAM Unsorted \
        --runThreadN $task.cpus \
        --outFilterMultimapNmax 1 \
        --outFilterMatchNmin 0 \
        --readFilesCommand zcat

    samtools sort $launchDir/out/star_out/${prfx}/${prfx}_Aligned.out.bam -o $launchDir/out/star_out/$prfx/${prfx}_Aligned_Sorted.out.bam --threads $task.cpus
    samtools index $launchDir/out/star_out/${prfx}/${prfx}_Aligned_Sorted.out.bam
    """
}

process runfeatureCounts {
    beforeScript 'source ~/.profile'
    cpus 4
    input:
        val bamfile
    output:
        val "$launchDir/out/featureCounts_out/$prfx/${prfx}_Aligned_Sorted.out.bam.featureCounts.bam"
    """
    featureCounts \
        -s 1 \
        -a $gtf \
        -o $launchDir/out/featureCounts_out/${prfx}/${prfx} \
        -R BAM \
        -T $task.cpus \
        -t gene \
        $bamfile
    """
}

process add_tag {
    beforeScript 'source ~/.profile'

    input:
        val bamfile
    
    """
    python /SGRNJ06/randd/USER/liuzihao/mydata/shells/add_gene_tag.py --inbam $bamfile --outbam $launchDir/out/featureCounts_out/${prfx}/${prfx}_featureCounts_tag_final.bam --gtf $gtf
    """
}

workflow {
    def reads_ch = Channel.fromPath(params.readsFile)
    runSTAR(reads_ch, genome) | runfeatureCounts | add_tag
}

