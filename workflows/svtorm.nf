#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                        IMPORT PLUGINS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                                                              } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                        IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DELLY                                                                         } from '../modules/local/delly/main'
include { MANTA                                                                         } from '../modules/local/manta/main'
include { SVABA                                                                         } from '../modules/local/svaba/main'
include { CROSSV                                                                        } from '../modules/local/crossv/main'
include { DRAWSV                                                                        } from '../modules/local/drawsv/main'
include { GRIDSS                                                                        } from '../modules/local/gridss/main'
include { SVSOMF                                                                        } from '../modules/local/svsomf/main'
include { MULTIQC                                                                       } from '../modules/nf-core/multiqc/main'
include { RECALL_SV                                                                     } from '../modules/local/recallsv/main'
include { BAM_PAIRED                                                                    } from '../modules/local/bampaired/main'
include { IANNOTATESV                                                                   } from '../modules/local/iannotatesv/main'
include { SURVIVOR_MERGE                                                                } from '../modules/local/survivor/merge/main'
include { SURVIVOR_STATS                                                                } from '../modules/local/survivor/stats/main'
include { SURVIVOR_FILTER                                                               } from '../modules/local/survivor/filter/main'
include { SERACARE_CHECKUP                                                              } from '../modules/local/seracare/checkup/main'
include { GATK4_BEDTOINTERVALLIST                                                       } from '../modules/local/gatk4/bedtointervallist/main'
include { PICARD_COLLECTMULTIPLEMETRICS                                                 } from '../modules/nf-core/picard/collectmultiplemetrics/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                     IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMultiqc                                                          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                                        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                                        } from '../subworkflows/local/utils_nfcore_svtorm_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    CREATE CUSTOM CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_fai                                          = Channel.fromPath(params.fai).map                              { it -> [[id:it.Name], it] }.collect()
ch_dict                                         = Channel.fromPath(params.dict).map                             { it -> [[id:it.Name], it] }.collect()
ch_fasta                                        = Channel.fromPath(params.fasta).map                            { it -> [[id:it.Name], it] }.collect()
ch_known_sites                                  = Channel.fromPath(params.known_sites).map                      { it -> [[id:it.Name], it] }.collect()
ch_targets_bed                                  = Channel.fromPath(params.intervals_bed_gunzip).map             { it -> [[id:it.Name], it] }.collect()
ch_known_sites_tbi                              = Channel.fromPath(params.known_sites_tbi).map                  { it -> [[id:it.Name], it] }.collect()
ch_targets_bed_tbi                              = Channel.fromPath(params.intervals_bed_gunzip_index).map       { it -> [[id:it.Name], it] }.collect()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                      RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SVTORM {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

   //
    // MODULE: Determine BAM pairedness for fastq conversion
    //
    BAM_PAIRED(ch_samplesheet)
    ch_versions = ch_versions.mix(BAM_PAIRED.out.versions)
    ch_bam_with_pairedness = BAM_PAIRED.out.bam_bai_issingleend.map { meta, bam, bai, single_end ->
        meta["single_end"] = single_end.text.trim().toBoolean()
        [meta, bam, bai]
    }

    //
    // Extract bed files (prior) from samplesheet for CROSSV
    //
    ch_bed_files = ch_samplesheet
        .map { meta, files ->
            def bed_path = meta.containsKey('prior') ? meta.prior : null
            if (bed_path) {
                try {
                    def bed_file = file(bed_path.toString())
                    if (!bed_file.exists()) {
                        log.warn "WARNING: Bed file does not exist: ${bed_path}"
                    }
                    [meta.patient, bed_file]
                } catch (Exception e) {
                    log.warn "WARNING: Could not create file object from bed path: ${bed_path} - ${e.message}"
                    [meta.patient, null]
                }
            } else {
                log.warn "WARNING: No 'prior' field found in metadata for patient ${meta.patient}, sample ${meta.id}"
                [meta.patient, null]
            }
        }
        .filter { patient, bed_file -> bed_file != null }
        .unique()
        .set { ch_bed_files_checked }

    //
    // MODULE: Run Picard Collect Multiple Metrics
    //
    PICARD_COLLECTMULTIPLEMETRICS(ch_bam_with_pairedness, ch_fasta, ch_fai)
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)
    ch_multiqc_files  = ch_multiqc_files.mix(PICARD_COLLECTMULTIPLEMETRICS.out.metrics.collect{it[1]}.ifEmpty([]))

    //
    // Pair tumour samples with normal samples by patient meta key if only tumour use backup normal
    //
    def control_normal_bam = file(params.normal_uncollapsed_bam)
    def control_normal_bai = file(params.normal_uncollapsed_bai)

    ch_bam_with_pairedness
        .map { meta, bam, bai -> [meta.patient, meta, bam, bai] }
        .groupTuple(by: 0)
        .map { patient, meta_list, bam_list, bai_list -> 
            def samples = []
            meta_list.eachWithIndex { meta, i ->
                samples << [meta, bam_list[i], bai_list[i]]
            }

            def tumour = samples.find { it[0].tumour == true  }
            def normal = samples.find { it[0].tumour == false }

            if (tumour && normal) {
                return [tumour, normal]
            } else if (tumour && tumour[0].matched == false) {
                return [tumour, null]
            } else {
                return null
            }
        }
        .filter { it != null }
        .map { tumour, normal ->
            def meta_t = tumour[0]
            def bam_t  = tumour[1]
            def bai_t  = tumour[2]

            if (normal) {
                def bam_n = normal[1]
                def bai_n = normal[2]
                return [meta_t, bam_t, bai_t, bam_n, bai_n]
            } else {
                return [meta_t, bam_t, bai_t, control_normal_bam, control_normal_bai]
            }
        }
        .set { ch_bam_pairs }

    //
    // MODULE: Run Delly Call
    //
    DELLY(ch_bam_pairs, ch_fasta, ch_fai, params.exclude_bed)
    ch_versions = ch_versions.mix(DELLY.out.versions)
    ch_delly_vcf = DELLY.out.vcf
    ch_delly_vcf = ch_delly_vcf.map { meta, vcf -> tuple(meta.patient, meta, vcf) }

    //
    // MODULE: Run Gridds (Extract overlapping fragments & calling)
    //
    GRIDSS(ch_bam_pairs, ch_fai, ch_fasta, params.intervals, params.blocklist_bed, params.bwa, params.kraken2db)
    ch_versions = ch_versions.mix(GRIDSS.out.versions)
    ch_gridss_vcf = GRIDSS.out.vcf
    ch_gridss_vcf = ch_gridss_vcf.map { meta, vcf -> tuple(meta.patient, meta, vcf) }

    //
    // MODULE: Run Manta in Somatic Mode
    //
    MANTA(ch_bam_pairs, ch_targets_bed, ch_targets_bed_tbi, ch_fasta, ch_fai, [])
    ch_versions = ch_versions.mix(MANTA.out.versions)
    ch_manta_vcf = MANTA.out.vcf
    ch_manta_vcf = ch_manta_vcf.map { meta, vcf -> tuple(meta.patient, meta, vcf) }

    //
    // MODULE: Run SvABA Note: version 1.2.0
    //
    SVABA(ch_bam_pairs, ch_fasta, ch_fai, params.known_sites_tbi, params.known_sites, params.intervals, params.bwa)
    ch_versions = ch_versions.mix(SVABA.out.versions)
    ch_svaba_vcf = SVABA.out.vcf
    ch_svaba_vcf = ch_svaba_vcf.map { meta, vcf -> tuple(meta.patient, meta, vcf) }

    //
    // Combine the vcf by patient
    //
    ch_vcf_merged = ch_delly_vcf
        .join(ch_gridss_vcf)
        .join(ch_manta_vcf)
        .join(ch_svaba_vcf)
        .map { patient, meta_delly, delly_vcf, meta_gridss, gridss_vcf, meta_manta, manta_vcf, meta_svaba, svaba_vcf ->
            tuple(
                meta_delly, 
                meta_delly,  delly_vcf,
                meta_gridss, gridss_vcf,
                meta_manta,  manta_vcf,
                meta_svaba,  svaba_vcf
            )
        }
    
    //
    // MODULE: Run Survivor to merge Unfiltered VCFs
    //
    SURVIVOR_MERGE(ch_vcf_merged, params.chromosomes, 1000, 2, 1, 1, 0, 250)
    ch_versions = ch_versions.mix(SURVIVOR_MERGE.out.versions)
    ch_merged_bed = SURVIVOR_MERGE.out.bed
    ch_merged_vcf = SURVIVOR_MERGE.out.vcf

    //
    // MODULE: Run GATK4 to Convert BED into Interval List
    //
    GATK4_BEDTOINTERVALLIST(ch_merged_bed, ch_dict)
    ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions)
    ch_merged_int_list = GATK4_BEDTOINTERVALLIST.out.interval_list
    ch_merged_int_list = ch_merged_int_list.map { meta, interval_list -> tuple(meta.patient, meta, interval_list) }

    //
    // Join interval lists with BAM pairs based on patient
    //
    ch_bam_pairs_by_patient = ch_bam_pairs.map {meta, bam_t, bai_t, bam_n, bai_n -> tuple(meta.patient, meta, bam_t, bai_t, bam_n, bai_n) }
    ch_recall_input = ch_bam_pairs_by_patient
        .join(ch_merged_int_list)
        .map { patient, meta_b, bam_t, bai_t, bam_n, bai_n, meta_i, interval_list ->
            tuple(
                meta_b, 
                meta_b, bam_t, bai_t, 
                meta_b, bam_n, bai_n, 
                meta_i, interval_list
            )
        }

    //
    // MODULE: Run Gridds in ReCall mode
    //
    RECALL_SV(ch_recall_input, ch_fasta, ch_fai, ch_known_sites, ch_known_sites_tbi, params.blocklist_bed, params.bwa, params.kraken2db, params.pon_directory, params.refflat, params.intervals)
    ch_versions = ch_versions.mix(RECALL_SV.out.versions)
    ch_prefilter_vcf = RECALL_SV.out.vcf

    //
    // MODULE: Run SV Somatic Filter
    //
    SVSOMF(ch_prefilter_vcf, ch_fasta, ch_fai, params.pon_directory)
    ch_versions = ch_versions.mix(SVSOMF.out.versions)
    ch_recall_vcf = SVSOMF.out.vcf
    ch_recall_vcf = ch_recall_vcf.map { meta, vcf -> tuple(meta.patient, meta, vcf) }

    //
    // MODULE: Run Survivor to filter Unfiltered VCFs
    //
    ch_survivor_filter_input = ch_delly_vcf
        .join(ch_gridss_vcf)
        .join(ch_manta_vcf)
        .join(ch_recall_vcf)
        .join(ch_svaba_vcf)
        .map { patient, meta_delly, delly_vcf, meta_gridss, gridss_vcf, meta_manta, manta_vcf, meta_recall, recall_vcf, meta_svaba, svaba_vcf ->
            tuple(
                meta_delly, 
                meta_delly , delly_vcf,
                meta_gridss, gridss_vcf,
                meta_manta , manta_vcf,
                meta_recall, recall_vcf,
                meta_svaba , svaba_vcf
            )
        }
    SURVIVOR_FILTER(ch_survivor_filter_input, 1000, 3, 1, 1, 0, 250, params.allowlist_bed)
    ch_versions = ch_versions.mix(SURVIVOR_FILTER.out.versions)
    ch_filtered_all = SURVIVOR_FILTER.out.filtered_all
    ch_filtered_vcf = SURVIVOR_FILTER.out.filtered_vcf
    ch_filtered_tsv = SURVIVOR_FILTER.out.filtered_tsv

    //
    // MODULE: Run Survivor Stats
    //
    SURVIVOR_STATS(ch_filtered_vcf, -1, -1, -1)
    ch_versions = ch_versions.mix(SURVIVOR_STATS.out.versions)
    ch_multiqc_files  = ch_multiqc_files.mix(SURVIVOR_STATS.out.stats.collect{it[1]}.ifEmpty([]))

    //
    // MODULE: Run CrosSV
    //
    ch_bam_pairs_for_crossv = ch_bam_pairs
        .map { meta, bam_t, bai_t, bam_n, bai_n -> 
            tuple(meta.patient, meta, bam_t, bai_t, bam_n, bai_n) 
        }
    ch_filtered_tsv_for_crossv = ch_filtered_tsv
        .map { meta, tsv -> tuple(meta.patient, meta, tsv) }
    ch_crossv_input = ch_bam_pairs_for_crossv
        .join(ch_delly_vcf)
        .join(ch_gridss_vcf)
        .join(ch_manta_vcf)
        .join(ch_recall_vcf)
        .join(ch_svaba_vcf)
        .join(ch_filtered_tsv_for_crossv)
        .join(ch_bed_files_checked)
        .map { patient, meta_bam, bam_t, bai_t, bam_n, bai_n, 
               meta_delly, delly_vcf, meta_gridss, gridss_vcf, 
               meta_manta, manta_vcf, meta_recall, recall_vcf, 
               meta_svaba, svaba_vcf, meta_tsv, tsv, bed_file -> 
               [meta_bam, bam_t, bai_t, bam_n, bai_n, 
                delly_vcf, gridss_vcf, manta_vcf, recall_vcf, svaba_vcf, tsv, bed_file]
        }
    CROSSV(ch_crossv_input)
    ch_versions = ch_versions.mix(CROSSV.out.versions)
    ch_crossv_tsv = CROSSV.out.tsv
    ch_crossv_annote = CROSSV.out.annote
    ch_crossv_bam_bai = CROSSV.out.bam.join(CROSSV.out.bai)

    //
    // MODULE: Run iAnnotateSV 
    //
    ch_crossv_tsv_mapped = ch_crossv_tsv
        .map { meta, tsv -> tuple(meta.patient, meta, tsv) }
    ch_crossv_annote_mapped = ch_crossv_annote
        .map { meta, annote -> tuple(meta.patient, meta, annote) }
    ch_annotate_input_tsv = ch_crossv_tsv_mapped
        .ifEmpty { ch_filtered_tsv_for_crossv }
    ch_annotate_input_annote = ch_crossv_annote_mapped
        .ifEmpty { ch_filtered_all.map { meta, vcf, vcf_tbi, tsv, annote_input -> tuple(meta.patient, meta, annote_input) } }
    ch_annotate_input = ch_filtered_all
        .map { meta, vcf, vcf_tbi, tsv, annote_input -> tuple(meta.patient, meta, vcf, vcf_tbi, annote_input) }
        .join(ch_annotate_input_tsv)
        .join(ch_annotate_input_annote)
        .map { patient, meta_f, vcf, vcf_tbi, annote_input_old, meta_t, tsv, meta_a, annote_input ->
            tuple(meta_f, vcf, vcf_tbi, tsv, annote_input)
        }
    IANNOTATESV(ch_annotate_input, params.genome)
    ch_versions = ch_versions.mix(IANNOTATESV.out.versions)
    ch_annotated_tsv = IANNOTATESV.out.tsv
    ch_annotated_ann = IANNOTATESV.out.ann
    ch_annotated_tsv_with_svs = ch_annotated_tsv
        .filter { meta, tsv ->
            tsv.readLines().size() > 1
        }
        .map { meta, tsv -> tuple(meta.patient, meta, tsv) }

    //
    // MODULE: Run DrawSV
    //
    ch_bam_pairs_by_patient = ch_bam_pairs.map {meta, bam_t, bai_t, bam_n, bai_n -> tuple(meta.patient, meta, bam_t, bai_t, bam_n, bai_n) }
    ch_drawsv_input = ch_bam_pairs_by_patient
        .join(ch_annotated_tsv_with_svs)
        .map { patient, meta_b, bam_t, bai_t, bam_n, bai_n, meta_t, tsv ->
            tuple(
                meta_b, 
                meta_b, bam_t, bai_t, 
                meta_b, bam_n, bai_n, 
                meta_t, tsv
            )
        }
    DRAWSV(ch_drawsv_input, params.annotations, params.cytobands, params.drawsv_chr, params.protein_domains)
    ch_versions = ch_versions.mix(DRAWSV.out.versions)
    ch_drawsv_pdf = DRAWSV.out.pdf

     //
    // MODULE: Run SeraCare Check-Up
    //
    ch_seracare_sample = ch_annotated_tsv
        .filter { meta, file -> 
            meta.id.contains("SeraCare") || meta.id.contains("SRCR")
        }
    SERACARE_CHECKUP(ch_seracare_sample)
    ch_versions = ch_versions.mix(SERACARE_CHECKUP.out.versions)

    //
    // MODULE: MultiQC
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'svtorm_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()
    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList()
    versions       = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                            THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
