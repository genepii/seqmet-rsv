runpath="/srv/nfs/ngs-stockage/NGS_Virologie/seqmet/ncov/"
outputpath="/srv/nfs/ngs-stockage/NGS_Virologie/redmine/7408/"

#touch ${outputpath}BILAN.tsv
touch ${outputpath}all_consensus.fasta
mkdir -p ${outputpath}/summary/

for run in "$@";
do
    echo $run
    cat ${runpath}${run}/varcall/MN908947_consensus.fna >> ${outputpath}all_consensus.fasta
    cp -r ${runpath}${run}/*_summary.tsv ${outputpath}summary/${run}.tsv

done
