// fetchcosmic pipline

params.readsFile = ""
params.sp = "human" // 目前只能用于人的数据！

if( sp = "human" ) {
    genome = ""
}
else {
    genome = ""
}

myDir = file("$launchDir/out/star_out/")
myDir.mkdirs()


process runSTAR {
    beforeScript 'source ~/.profile'
    cpus 4
    input:
        path reads
        path genome
    output:
        val "$launchDir/out/out.bed"

    """
    STAR --genomeDir $genome --readFilesIn $reads --outFileNamePrefix "$launchDir/out/star_out/out_" --outSAMtype BAM SortedByCoordinate --runThreadN $task.cpus
    bedtools bamtobed -i $launchDir/out/star_out/out_Aligned.sortedByCoord.out.bam  > $launchDir/out/out.bed
    """
}

process fetchCosmic {
    beforeScript 'source ~/.profile'
    input:
        val bedfile
    """
    #!/SGRNJ/Public/Software/conda_env/cenv_lzh/bin/python
    from tk.fetch_cosmic import CosmicFtechAndAnnoate
    r = CosmicFtechAndAnnoate("$bedfile", "$launchDir/out/70")
    r.run()
    """
}

workflow {
    def reads_ch = Channel.fromPath(params.readsFile)
    runSTAR(reads_ch, genome) | fetchCosmic
}
