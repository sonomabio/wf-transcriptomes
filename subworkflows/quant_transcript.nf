// import build_minimap_index_transcriptome process from differential_expression.nf

include {build_minimap_index_transcriptome; map_transcriptome; count_transcripts; mergeCounts; mergeTPM} from './differential_expression.nf'


workflow quant_transcript {
    take: 
        ref_transcriptome
        full_len_reads
        sample_sheet
    main:
        sample_sheet = Channel.fromPath(sample_sheet)
        t_index = build_minimap_index_transcriptome(ref_transcriptome)
        mapped = map_transcriptome(full_len_reads.combine(t_index)
                .map{meta, fastq, reference, transcriptome -> tuple(meta, fastq, reference) })
        count_transcripts(mapped.bam.combine(t_index.map{ mmi, reference -> reference}))
        merged = mergeCounts(count_transcripts.out.counts.collect())
        merged_TPM = mergeTPM(count_transcripts.out.counts.collect())
    emit:
        merge_counts = merged
        merged_TPM = merged_TPM
}