// fetchcosmic pipline

params.fa = ""
params.ge = ""
params.p = ""

reads = params.fa
genome = params.ge
prfx = params.p

myDir = file("$launchDir/out/")
myDir.mkdirs()


process runSTAR2bed {
    beforeScript 'source ~/.profile'
    cpus 4
    input:
        path reads
        path genome

    output:
        val "$launchDir/out/${prfx}.bed"

    """
    STAR --genomeDir $genome --readFilesIn $reads --outFileNamePrefix "$launchDir/out/star_out/out_" --outSAMtype BAM SortedByCoordinate --runThreadN $task.cpus
    bedtools bamtobed -i $launchDir/out/star_out/out_Aligned.sortedByCoord.out.bam -bed12 > $launchDir/out/${prfx}.bed
    """
}

process bed2gtf {
    beforeScript 'source activate cellrank && source ~/.profile'
    input:
        val bedfile
    """
    bedToGenePred $bedfile $launchDir/out/${prfx}.pred
    genePredToGtf file $launchDir/out/${prfx}.pred $launchDir/out/${prfx}.gtf 
    
    sed -i 's/${prfx}.pred/{$prfx}_genome/g'  $launchDir/out/${prfx}.gtf 
    awk -F '\t' -v OFS='\t' '{if (\$3=="transcript") \$3="gene"; print \$0}' ${prfx}.gtf  > ${prfx}_2.gtf

    """
}

workflow {
    def reads_ch = Channel.fromPath(params.fa)
    def ge_ch = Channel.fromPath(params.ge)
    runSTAR2bed(reads_ch, ge_ch) | bed2gtf
}
