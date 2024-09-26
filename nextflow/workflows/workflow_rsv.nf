#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {interleave_bbmap} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {qc_fastqc} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {qc_fastqscreen} from "${params.nfpath}/nextflow/modules/module_rsv.nf"

include {filter_srahumanscrubber} from "${params.nfpath}/nextflow/modules/module_rsv.nf"

include {downsample_bbmap} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {trim_cutadapt} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {premap_bwamem} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {premap_verif} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {map_bwamem} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {postmap_picard} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {capdepth_jvarkit} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {postmap_abra2} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {postmap_ampliconclip} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {postmap_fixbam} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {callvar_bcftools} from "${params.nfpath}/nextflow/modules/module_rsv.nf"

include {gentsv_vcf} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {genlist_countsummary} from "${params.nfpath}/nextflow/modules/module_rsv.nf"

include {genbed_bedtools} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {genfasta_consensus} from "${params.nfpath}/nextflow/modules/module_rsv.nf"

include {gentable_variant} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {merge_variant} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {seek_contaminant} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {seek_coinf} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {merge_contaminant} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {merge_coverage} from "${params.nfpath}/nextflow/modules/module_rsv.nf"

include {merge_consensus} from "${params.nfpath}/nextflow/modules/module_rsv.nf"

include {merge_covsummary} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {gentsv_nextclade} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {merge_fastqstat} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {merge_bamstat} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {merge_countsummary} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {merge_coinfsummary} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {merge_poscsummary} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {merge_veriftable} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {gentsv_summary} from "${params.nfpath}/nextflow/modules/module_rsv.nf"

include {merge_summary} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {gentsv_validation} from "${params.nfpath}/nextflow/modules/module_rsv.nf"

include {merge_fastfinder} from "${params.nfpath}/nextflow/modules/module_rsv.nf"

include {qc_multiqc} from "${params.nfpath}/nextflow/modules/module_rsv.nf"

include {merge_bbmerge} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {assign_kraken2} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {assign_tax2lin} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {import_q2table} from "${params.nfpath}/nextflow/modules/module_rsv.nf"

include {qc_covplot} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {merge_covplot} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {gentsv_blast} from "${params.nfpath}/nextflow/modules/module_rsv.nf"
include {merge_blast} from "${params.nfpath}/nextflow/modules/module_rsv.nf"

workflow denovo {
    take:
        ch_fastq

    main:
    
        ch_fastq.map{ it[0] }.set{ ch_sampleId }

        if ( params.readtype == 'paired-end' ) {
            interleave_bbmap(ch_fastq)
            interleave_bbmap.out.set{ ch_fastq }
        }

        qc_fastqc(ch_fastq)
        qc_fastqc.out.fastqc.set{ ch_multiqc }
        qc_fastqc.out.fastqstat.map{ it[1] }.set{ ch_fastqstat }

        qc_fastqscreen(ch_fastq)
        ch_multiqc.mix(qc_fastqscreen.out).set{ ch_multiqc }

        if ( params.downsample ) {
            downsample_bbmap(ch_fastq)
            downsample_bbmap.out.set{ ch_fastq }
        }

        if ( params.dehost ) {
            filter_srahumanscrubber(ch_fastq)
            filter_srahumanscrubber.out.set{ ch_fastq }
        }

        if ( params.trim ) {
            trim_cutadapt(ch_fastq)
            trim_cutadapt.out.fastq.set{ ch_fastq }
            ch_multiqc.mix(trim_cutadapt.out.qc).set{ ch_multiqc }
        }

        if ( params.readtype == 'paired-end' ) {
            merge_bbmerge(ch_fastq)
        }

        if ( params.readtype == 'paired-end' ) {
            merge_bbmerge.out.set{ ch_fastqDeNovo }
        } else {
            ch_fastq.set{ ch_fastqDeNovo }
        }

        if ( params.assign == 'kraken2' ) {
            assign_kraken2(ch_fastqDeNovo)
            assign_tax2lin(assign_kraken2.out.k2report)
            import_q2table(assign_tax2lin.out.collect())
        } 

        if ( params.refmode == 'list' ) {
            ch_sampleId.combine(Channel.fromPath(params.premap_ref).collect().toSortedList()).set{ ch_premapref }
            ch_sampleId.combine(Channel.fromPath(params.premap_poscref).collect().toSortedList()).set{ ch_premapposcref }
            ch_sampleId.combine(Channel.fromPath(params.map_ref)).set{ ch_mapref }
        }

        ch_mapref.map{ it[1..-1] }.collect().flatten().unique().map{ [ (it =~  /(.+)\/([^\.\/]+)\.(.+)/)[0][2] , it ] }.set{ ch_summary }

        if ( params.premap == 'bwamem' ) {
            premap_bwamem(ch_fastq.join(ch_premapref).join(ch_premapposcref))
            premap_bwamem.out.coveragetable.set{ ch_premap_coveragetable }
            premap_bwamem.out.poscsummary.set{ ch_premap_poscsummary }
        } 

        if ( params.map == 'bwamem' ) {
            map_bwamem(ch_fastq.join(ch_premap_coveragetable).combine(ch_mapref, by: 0))
            map_bwamem.out.bam.set{ ch_bamTuple }
            map_bwamem.out.fastqstat.map{ [it[0], it[2]] }.set{ ch_bamstat }
        }

        if ( params.capdepth ) {
            capdepth_jvarkit(ch_bamTuple)
            capdepth_jvarkit.out.set{ ch_bamTuple }
        }

        if ( params.realign ) {
            postmap_abra2(ch_bamTuple)
            postmap_abra2.out.set{ ch_bamTuple }
        }

        if ( params.clipbam == 'ampliconclip' ) {
            postmap_ampliconclip(ch_bamTuple.combine(Channel.from(params.postmap_ampliconclip.primer), by: 0))
            postmap_ampliconclip.out.set{ ch_bamTuple }
        }

        if ( params.fixbam ) {
            postmap_fixbam(ch_bamTuple)
            postmap_fixbam.out.bam.set{ ch_bamTuple }
            ch_bamstat.mix(postmap_fixbam.out.bamstat.map{ [it[0], it[2]] }).set{ ch_bamstat }
        }

        if ( params.calling == 'bcftools' ) {
            callvar_bcftools(ch_bamTuple.combine(Channel.from(params.callvar_bcftools.toexclude), by: 0))
            callvar_bcftools.out.set{ ch_calling }
        } 

        gentsv_vcf(ch_calling.gentsv)

        genlist_countsummary(ch_calling.gentsv)

        genbed_bedtools(ch_bamTuple.join(ch_calling.rawvcf, by:[0,1]))
        ch_bamstat.mix(genbed_bedtools.out.bamstat).set{ ch_bamstat }

        genfasta_consensus(ch_calling.callcons.join(genbed_bedtools.out.callcons, by:[0,1]).combine(Channel.from(params.empty)))

        gentsv_blast(genfasta_consensus.out.consensus)

        qc_covplot(genbed_bedtools.out.callcons)

        merge_consensus(genfasta_consensus.out.consensus.map{ [it[0], it[2]] }.groupTuple())
        ch_summary.join(merge_consensus.out).set{ ch_summary }

        if ( params.premap_verifref.size() > 0 ) {
            ch_sampleId.combine(Channel.fromPath(params.premap_verifref).collect().toSortedList()).set{ ch_premapverifref }
            premap_verif(ch_fastq.join(ch_premapverifref))
            premap_verif.out.veriftable.set{ ch_premap_veriftable }
            merge_veriftable(ch_premap_veriftable.collect())
            ch_summary.combine(merge_veriftable.out).set{ ch_summary }
        }

        if ( params.refmode == 'list' ) {
            gentable_variant(ch_calling.gentsv.combine(Channel.from(params.gentable_variant.gff), by: 0))
            merge_variant(gentable_variant.out.table.groupTuple().combine(Channel.from(params.merge_variant.level)))

            seek_contaminant((ch_calling.contaminated.join(genbed_bedtools.out.contaminated, by:[0,1])).combine(ch_calling.contaminant.groupTuple(), by:[0]))
            seek_coinf((ch_calling.contaminated.join(genbed_bedtools.out.contaminated, by:[0,1])).combine(Channel.from(params.seek_coinf.vcf), by:[0]))

            merge_contaminant(seek_contaminant.out.table.groupTuple())
            ch_summary.join(merge_contaminant.out.merged).set{ ch_summary }

            merge_coinfsummary(seek_coinf.out.raw.groupTuple())
            if ( params.seek_coinf.vcf.size() > 0 ) {
                ch_summary.join(merge_coinfsummary.out).set{ ch_summary }
            }
        }

        merge_coverage(ch_premap_coveragetable.map{ it[1] }.collect())

        merge_covsummary(genbed_bedtools.out.covsummary.groupTuple())
        ch_summary.join(merge_covsummary.out).set{ ch_summary }

        merge_fastqstat(ch_fastqstat.collect())
        ch_summary.combine(merge_fastqstat.out).set{ ch_summary }

        merge_bamstat(ch_bamstat.groupTuple())
        ch_summary.join(merge_bamstat.out).set{ ch_summary }

        merge_countsummary(genlist_countsummary.out.groupTuple())
        ch_summary.join(merge_countsummary.out).set{ ch_summary }

        merge_poscsummary(ch_premap_poscsummary.collect())
        ch_summary.combine(merge_poscsummary.out).set{ ch_summary }

        merge_blast(gentsv_blast.out.blast.map{ [it[0], it[2]] }.groupTuple())
        ch_summary.join(merge_blast.out).set{ ch_summary }

        merge_covplot(qc_covplot.out.groupTuple())

        gentsv_summary(ch_summary.map{ [it[0], it[1..-1]] })

        merge_summary(gentsv_summary.out.map{ it[1] }.collect())

        qc_multiqc(ch_multiqc.map{ it[1] }.collect())
        
        gentsv_nextclade(merge_consensus.out.groupTuple().join(Channel.from(params.gentsv_nextclade.dataset)))

        if ( params.gentsv_validation["mode"] == 'rsv' ) {
            gentsv_validation(gentsv_summary.out.join(gentsv_nextclade.out).join(merge_variant.out).map{ [it[0], it[1..-1]] }.join(Channel.from(params.gentsv_validation.matrix)))
            merge_fastfinder(gentsv_validation.out.validation.collect())
        }

}