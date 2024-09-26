#!/bin/bash

WORK_PATH="${PWD}"

declare -a TARGETS=$TARGETS
if [[ "${TARGETS[0]}" == "" ]]; then TARGETS=( WG ); fi
if [[ "${TARGETS[0]}" != "WG" ]]; then COVCOL="_${TARGETS[0]}"; else COVCOL=""; fi
if [[ "$RUNREF" == "ALL" ]]; then runref=$(echo "$@" | sed "s/\([^ ]\+\)\/\(\ \|$\)/\1 /g" | sed 's/ *$//g'); else runref=$(echo "$1" | sed "s/\([^ ]\+\)\/\(\ \|$\)/\1/g"); fi

rm -rf "${WORK_PATH}/comp_summary.tsv"
touch "${WORK_PATH}/comp_summary.tsv"

for ARGS in "$@"
do
PIPERUN="${ARGS%/}"
echo "${PIPERUN} - $(date) - Starting"
echo "${PIPERUN} - $(date) - Starting" >> "${WORK_PATH}/consensuscomp.log"
cd "${ARGS}"
    for TSV in *_summary.tsv
    do
    LOG=""
    SAMPLEID_COL=$(sed -n $'1s/\\\t/\\\n/gp' "${TSV}" | grep -nx "sample_id" | cut -d: -f1) || LOG="error"
    REFERENCEID_COL=$(sed -n $'1s/\\\t/\\\n/gp' "${TSV}" | grep -nx "reference_id" | cut -d: -f1) || LOG="error"
    PERCCOVERAGE_COL=$(sed -n $'1s/\\\t/\\\n/gp' "${TSV}" | grep -nx "consensus_perccoverage${COVCOL}" | cut -d: -f1) || LOG="error"
    SUP_COL=$(sed -n $'1s/\\\t/\\\n/gp' "${TSV}" | grep -nx "vcf_dpcount" | cut -d: -f1) || LOG="error"
    #GLIMSID="$( printf $SAMPLEID | cut -d$'-' -f2 | cut -d$'_' -f1  | sed -e 's/\([0-9]*\)\(.*\)/\1/g' )"
    singularity exec /srv/scratch/iai/seqmet/singularity/varcall.sif awk -v piperun=${PIPERUN} -v sampleid_col="${SAMPLEID_COL}" -v referenceid_col="${REFERENCEID_COL}" -v perccoverage_col="${PERCCOVERAGE_COL}" -v sup_col="${SUP_COL}" 'BEGIN { FS="\t"; OFS="\t"} FNR > 1 { split($sampleid_col,splitA,"-"); split(splitA[2],splitB,"_"); glimsid=splitB[1]; glimsid=gensub(/([0-9]*)(.*)/, "\\1", "g", glimsid); if ($perccoverage_col>=90 && $perccoverage_col!="NA" && $referenceid_col!="NA" && $sup_col<999999) print piperun,glimsid,$sampleid_col,$referenceid_col,$perccoverage_col }' "${TSV}" >> "${WORK_PATH}/comp_summary.tsv"
    done
done

cd "${WORK_PATH}"
tac "${WORK_PATH}/comp_summary.tsv" | awk -F"\t" '!_[$3]++' > "${WORK_PATH}/comp_summary_uniq.tsv"
cat "${WORK_PATH}/comp_summary_uniq.tsv" | awk -v runref="$runref" 'BEGIN { FS="\t"; OFS="\t"; split(runref, runlist, " "); for (i in runlist) runvalues[runlist[i]] = "" } $1 in runvalues { print $2 }' > comp_runref.txt
cat "${WORK_PATH}/comp_summary_uniq.tsv" | grep -f comp_runref.txt | sort -t$'\t' -k2,2 | awk 'BEGIN { FS="\t"; OFS="\t"} { count[$2]++ } END { for ( i in count ) { if (count[i]>=2 && i!="") print i } }' > "${WORK_PATH}/comp_summary_uniq.txt"

while read GLIMSID
do
grep "\-${GLIMSID}" "${WORK_PATH}/comp_summary_uniq.tsv" > "${WORK_PATH}/current_comp.tsv"
    while read LINE
    do
    RUNPATH="$(printf "${LINE}" | cut -d$'\t' -f1)"
    GLIMSID="$(printf "${LINE}" | cut -d$'\t' -f2)"
    SAMPLEID="$(printf "${LINE}" | cut -d$'\t' -f3)"
    REFERENCE="$(printf "${LINE}" | cut -d$'\t' -f4)"
    cd "$RUNPATH"
        for TARGET in ${TARGETS[*]}
        do
        cat "varcall/cons/${REFERENCE}/${SAMPLEID}.fna" | grep "${TARGET}$" -A1 >> "${WORK_PATH}/${REFERENCE}_${GLIMSID}_${TARGET}.fna"
        done
    done<"${WORK_PATH}/current_comp.tsv"
done<"${WORK_PATH}/comp_summary_uniq.txt"

cd "${WORK_PATH}"

for FNA in *.fna
do
    if [[ "$(cat "${FNA}" | grep ">" | wc -l)" -le 1 ]]
    then
    rm -rf "${FNA}"
    else
    singularity exec /srv/scratch/iai/seqmet/singularity/varcall.sif mafft --thread 16 --threadtb 8 --threadit 0 --adjustdirection --anysymbol --op 5.0 --auto "${FNA}" | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | sed -e "s/\t\(.\{100\}\)/\t$(printf 'N%0.s' {0..100})/g" | sed -e "s/\(.\{100\}\)$/$(printf 'N%0.s' {0..100})/g" | tr '\t' '\n' > "${FNA%%.fna}_mafft.fasta"
    singularity exec /srv/scratch/iai/seqmet/singularity/varcall.sif python /srv/scratch/iai/seqmet/script/fasta_describeregions.py -f "${FNA%%.fna}_mafft.fasta" -o "${FNA%%.fna}_mafftb.tsv" -x "${FNA%%.fna}_mafftb.fasta"
    singularity exec /srv/scratch/iai/seqmet/singularity/matdist.sif ruby /srv/scratch/iai/seqmet/script/regionalrealignment.rb "${FNA%%.fna}_mafftb.tsv" "${FNA%%.fna}_mafftb.fasta" > "${FNA%%.fna}_mafftc.fasta"
    singularity exec /srv/scratch/iai/seqmet/singularity/matdist.sif Rscript /srv/scratch/iai/seqmet/script/gen_matdist.R --input "${FNA%%.fna}_mafftb.fasta" --matrix "${FNA%%.fna}_mafft.pdf" --table "${FNA%%.fna}_mafft.tsv"
    fi
done
