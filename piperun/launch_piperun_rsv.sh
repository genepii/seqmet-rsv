#!/bin/bash

PIPERUN_DIR="${PWD}"
NF_PATH="${PIPERUN_DIR%%/piperun}"

for ARGS in "$@"
do
PIPERUN="${ARGS%/}"
echo "${PIPERUN} - $(date) - Starting"
echo "${PIPERUN} - $(date) - Starting" >> "${NF_PATH}/piperun/launch_piperun.log"
cd "${NF_PATH}/piperun/${PIPERUN}"
RUN_ID="${ARGS%/}"
DATA_TYPE="${PIPERUN##*_}"
    for JSON in *.json
    do
    rm -rf work .nextflow* report.html*
    LOG=""
    "${NF_PATH}/nextflow/nextflow" -C "${NF_PATH}/nextflow/nextflow_rsv.config" run "${NF_PATH}/nextflow/main_rsv.nf" -params-file "${JSON}" -with-trace --prefix "${RUN_ID}" || LOG="error"
    #chmod -fR 777 "${NF_PATH}/piperun/${PIPERUN}"
    if [[ "${LOG}" != "" ]]; then echo "${RUN_ID} - $(date) - Failed" >> "${NF_PATH}/piperun/launch_piperun.log"; continue; fi
    head -n 1 trace.txt > "${NF_PATH}/result/${PIPERUN}/${DATA_TYPE}_trace.tsv"
    for f in trace*; do tail -n+2 "${f}" >> "${NF_PATH}/result/${PIPERUN}/${DATA_TYPE}_trace.tsv"; done
    cp -r "${JSON}" "${NF_PATH}/result/${PIPERUN}/"
        if [[ "${PIPERUN%%_*}" == "000000" ]]
        then
        echo "${PIPERUN} - $(date) - Completed" >> "${NF_PATH}/piperun/launch_piperun.log"
        for f in ${NF_PATH}/result/${PIPERUN}/interleaved/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/merged/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/trimmed/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/downsampled/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/filtered/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/metashotgun/kraken2/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/varcall/bowtie2/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/varcall/bwamem/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/varcall/minimap2/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/varcall/picard/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/varcall/abra2/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        #find "${NF_PATH}/result/${PIPERUN}/" -type f -exec chmod 755 {} + -o -type d -exec chmod 777 {} +
        #mkdir -p "/srv/net/ngs-stockage.chu-lyon.fr/NGS_IAI/redmine/${DATA_TYPE}/"
        #cp -prn "${NF_PATH}/result/${PIPERUN}/" "/srv/net/ngs-stockage.chu-lyon.fr/NGS_IAI/redmine/${DATA_TYPE}/" && echo "${PIPERUN} - $(date) - Copied" >> "${NF_PATH}/piperun/launch_piperun.log"
        elif [[ ${DATA_TYPE} == "rsv"  ]]
        then
        echo "${PIPERUN} - $(date) - Completed" >> "${NF_PATH}/piperun/launch_piperun.log"
        for f in ${NF_PATH}/result/${PIPERUN}/interleaved/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/merged/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/trimmed/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/downsampled/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/filtered/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/metashotgun/kraken2/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/varcall/bowtie2/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/varcall/bwamem/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/varcall/minimap2/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/varcall/picard/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/varcall/abra2/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        #find "${NF_PATH}/result/${PIPERUN}/" -type f -exec chmod 755 {} + -o -type d -exec chmod 777 {} +
        mkdir -p "/srv/net/ngs-stockage.chu-lyon.fr/NGS_Virologie/seqmet/${DATA_TYPE}/"
        cp -prn "${NF_PATH}/result/${PIPERUN}/" "/srv/net/ngs-stockage.chu-lyon.fr/NGS_Virologie/seqmet/${DATA_TYPE}/" && echo "${PIPERUN} - $(date) - Copied" >> "${NF_PATH}/piperun/launch_piperun.log"
        else
        echo "${PIPERUN} - $(date) - Completed" >> "${NF_PATH}/piperun/launch_piperun.log"
        for f in ${NF_PATH}/result/${PIPERUN}/interleaved/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/merged/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/trimmed/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/downsampled/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/filtered/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/metashotgun/kraken2/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/varcall/bowtie2/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/varcall/bwamem/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/varcall/minimap2/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/varcall/picard/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/varcall/abra2/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        #find "${NF_PATH}/result/${PIPERUN}/" -type f -exec chmod 755 {} + -o -type d -exec chmod 777 {} +
        mkdir -p "/srv/net/ngs-stockage.chu-lyon.fr/NGS_IAI/redmine/${DATA_TYPE}/"
        cp -prn "${NF_PATH}/result/${PIPERUN}/" "/srv/net/ngs-stockage.chu-lyon.fr/NGS_IAI/redmine/${DATA_TYPE}/" && echo "${PIPERUN} - $(date) - Copied" >> "${NF_PATH}/piperun/launch_piperun.log"
        fi
    done
done
