include {retrieveMetadata} from "${params.nfpath}/nextflow/modules/util.nf"

process interleave_bbmap {

    // Convert paired-end fastq to an interleaved fastq

    label 'varcall'
    storeDir params.result
    echo false
    tag { sampleId }

    when:
    params.interleave_bbmap.todo == 1

    input:
    tuple(val(sampleId), path(R1), path(R2))
    
    output:
    tuple val(sampleId), path("interleaved/${sampleId}.fastq.gz")
    
    script:
    memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p interleaved/

    if [[ ! -s ${R1} ]]; then
        touch "interleaved/${sampleId}.fastq.gz"
        
        exit 0
    fi

    reformat.sh -Xmx${memory}G threads=$task.cpus in1="${R1}" in2="${R2}" out="interleaved/${sampleId}.fastq.gz" qin=auto qout=33
    
    rm -rf work
    """
}

process qc_fastqc {
    label 'varcall'
    storeDir params.result
    echo false
    tag { sampleId }

    when :
    params.qc_fastqc.todo == 1

    input :

    tuple(val(sampleId), path(fastq))

    output :

    tuple val(sampleId), path("fastqc/${sampleId}{_fastqc.zip,_fastqc.html,_R1_fastqc.zip,_R1_fastqc.html,_R2_fastqc.zip,_R2_fastqc.html}"), emit: fastqc
    tuple val(sampleId), path("fastqc/${sampleId}_readcount.tsv"), emit: fastqstat

    script :
    memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p fastqc work

    if [[ ! -s "${sampleId}.fastq.gz" ]]; then
        touch "fastqc/${sampleId}_fastqc.zip" "fastqc/${sampleId}_fastqc.html"
        echo -e "RAW\\tfastq_readcount\\t0\\t${sampleId}" > "fastqc/${sampleId}_readcount.tsv"
        
        exit 0
    fi

    if [[ ${params.readtype} == "paired-end" ]]; then
        reformat.sh -Xmx${memory}G threads=$task.cpus in="${sampleId}.fastq.gz" out1="work/${sampleId}_R1.fastq.gz" out2="work/${sampleId}_R2.fastq.gz" int=t qin=auto qout=33
        fastqc "work/${sampleId}_R1.fastq.gz" "work/${sampleId}_R2.fastq.gz" --outdir fastqc -t $task.cpus
        unzip -d work "fastqc/${sampleId}_R1_fastqc.zip"
        unzip -d work "fastqc/${sampleId}_R2_fastqc.zip"
    else
        fastqc "${sampleId}.fastq.gz" --outdir fastqc -t $task.cpus
        unzip -d work "fastqc/${sampleId}_fastqc.zip"
    fi

    cat \$(find work/ -type f -name fastqc_data.txt) >> "work/${sampleId}.txt"
    echo -e "RAW\\tfastq_readcount\\t\$(grep "Total Sequences" "work/${sampleId}.txt" | cut -f 2 | paste -sd+ | bc)\\t${sampleId}" > "fastqc/${sampleId}_readcount.tsv"

    rm -rf work
    """
}

process qc_fastqscreen {
    label 'varcall'
    storeDir params.result
    echo false
    tag { sampleId }

    when :
    params.qc_fastqscreen.todo == 1 && ['-TPOS', 'PCR', '-NT', 'NEG', 'VIDE'].any{ sampleId.toUpperCase().contains(it) }

    input:
    tuple(val(sampleId), path(fastq))

    output :
    tuple val(sampleId), path("fastqscreen/${sampleId}_screen.txt")

    script :
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8

    mkdir -p fastqscreen work

    if [[ ! -s "${sampleId}.fastq.gz" ]]; then
        touch "fastqscreen/${sampleId}.txt"
        
        exit 0
    fi

    fastq_screen --threads $task.cpus --conf ${params.qc_fastqscreen["conf"]} "${sampleId}.fastq.gz" --outdir "fastqscreen/"

    rm -rf work
    """
}

process merge_bbmerge {

    // Merge reads with an overlap, also output non-overlapping reads as an interleaved fastq

    label 'highmemory'
    storeDir params.result
    echo false
    tag { sampleId }

    when:
    params.merge_bbmerge.todo == 1

    input:
    tuple(val(sampleId), path(fastq))
    
    output:
    tuple val(sampleId), path("merged/${sampleId}{.fastq.gz,_RM.fastq.gz}")
    
    script:
    memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p merged

    if [[ ! -s "${sampleId}.fastq.gz" ]]; then
        touch "merged/${sampleId}.fastq.gz" "merged/${sampleId}_RM.fastq.gz"
        
        exit 0
    fi

    BBMERGE_ADAPTER_F="\$(bbmerge.sh in="${sampleId}.fastq.gz" outa=stdout | grep -A 1 'Read1' | grep -v '>' | grep -v '^--\$' | perl -pe 'chomp if eof' | tr '\\n' ',')"
    BBMERGE_ADAPTER_R="\$(bbmerge.sh in="${sampleId}.fastq.gz" outa=stdout | grep -A 1 'Read2' | grep -v '>' | grep -v '^--\$' | perl -pe 'chomp if eof' | tr '\\n' ',')"

    bbmerge.sh -Xmx${memory}G in="${sampleId}.fastq.gz" adapter1="\${BBMERGE_ADAPTER_F}" adapter2="\${BBMERGE_ADAPTER_R}" out="merged/${sampleId}_RM.fastq.gz" outu="merged/${sampleId}.fastq.gz" threads=$task.cpus int=t qin=auto qout=33 ecct
    
    rm -rf work
    """
}

process trim_cutadapt {

    // Trim adapters from each read, drop any read pair if at least one read is shorter than a specified length, multiple adapters can be provided

    label 'varcall'
    storeDir params.result
    echo false
    tag { sampleId }

    when:
    params.trim_cutadapt.todo == 1

    input:
    tuple(val(sampleId), path(fastq))
    
    output:
    tuple val(sampleId), path("trimmed/${sampleId}.fastq.gz"), emit: fastq
    tuple val(sampleId), path("trimmed/${sampleId}.txt"), emit : qc
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p work trimmed

    if [[ ! -s "${sampleId}.fastq.gz" ]]; then
        touch "trimmed/${sampleId}.fastq.gz" "trimmed/${sampleId}.txt"
        
        exit 0
    fi
    
    if [[ "${params.trim_cutadapt["adapter_f"]}" != [] || "${params.trim_cutadapt["adapter_r"]}" != [] ]]; then
        adapter_f=\$(for item in \$(echo ${params.trim_cutadapt["adapter_f"]} | tr -d '[\\[\\],\\"]'); do echo -n "-a \${item} "; done)
        adapter_r=\$(for item in \$(echo ${params.trim_cutadapt["adapter_r"]} | tr -d '[\\[\\],\\"]'); do echo -n "-A \${item} "; done)
    elif [[ "${params.readtype}" == "paired-end" ]]; then
        BBMERGE_ADAPTER_F="\$(bbmerge.sh in="${sampleId}.fastq.gz" outa=stdout | grep -A 1 'Read1' | grep -v '>' | grep -v '^--\$' | perl -pe 'chomp if eof' | tr '\\n' ',' | sed 's/\\(.\\{63\\}\\).*/\\1/')"
        BBMERGE_ADAPTER_R="\$(bbmerge.sh in="${sampleId}.fastq.gz" outa=stdout | grep -A 1 'Read2' | grep -v '>' | grep -v '^--\$' | perl -pe 'chomp if eof' | tr '\\n' ',' | sed 's/\\(.\\{63\\}\\).*/\\1/')"
        if [[ \${BBMERGE_ADAPTER_F} == "N" ]]; then BBMERGE_ADAPTER_F=""; fi
        if [[ \${BBMERGE_ADAPTER_R} == "N" ]]; then BBMERGE_ADAPTER_R=""; fi
        adapter_f=\$(for item in \$(echo \${BBMERGE_ADAPTER_F} | tr -d '[\\[\\],\\"]'); do echo -n "-a \${item} "; done)
        adapter_r=\$(for item in \$(echo \${BBMERGE_ADAPTER_R} | tr -d '[\\[\\],\\"]'); do echo -n "-A \${item} "; done)
    elif [[ "${params.readtype}" == "single-end" ]]; then
        BBMERGE_ADAPTER_F="\$(bbmerge.sh in="${sampleId}.fastq.gz" outa=stdout | grep -v '>' | grep -v '^--\$' | perl -pe 'chomp if eof' | tr '\\n' ',' | sed 's/\\(.\\{63\\}\\).*/\\1/')"
        TOEXCLUDE=("N" "N,N")
        if [[ \$(echo \${TOEXCLUDE[@]} | fgrep -w \$BBMERGE_ADAPTER_F) ]]; then BBMERGE_ADAPTER_F=""; fi
        adapter_f=\$(for item in \$(echo \${BBMERGE_ADAPTER_F} | tr -d '[\\[\\],\\"]'); do echo -n "-a \${item} "; done)
        adapter_r=""
    else
        adapter_f=""
        adapter_r=""
    fi

    front_f=\$(for item in \$(echo ${params.trim_cutadapt["front_f"]} | tr -d '[\\[\\],\\"]'); do echo -n "-g \${item} "; done)
    front_r=\$(for item in \$(echo ${params.trim_cutadapt["front_r"]} | tr -d '[\\[\\],\\"]'); do echo -n "-G \${item} "; done)
    anywhere_f=\$(for item in \$(echo ${params.trim_cutadapt["anywhere_f"]} | tr -d '[\\[\\],\\"]'); do echo -n "-b \${item} "; done)
    anywhere_r=\$(for item in \$(echo ${params.trim_cutadapt["anywhere_r"]} | tr -d '[\\[\\],\\"]'); do echo -n "-B \${item} "; done)

    if [[ ${params.readtype} == "paired-end" ]]; then
        ARGS="--interleaved --pair-filter=any \${adapter_f}\${adapter_r}\${front_f}\${front_r}\${anywhere_f}\${anywhere_r}"
    else
        ARGS="\${adapter_f}\${front_f}\${anywhere_f}"
    fi

    echo "ARGS : \$ARGS"
    echo "BBMERGE_ADAPTER_F : \$BBMERGE_ADAPTER_F"
    echo "BBMERGE_ADAPTER_R : \$BBMERGE_ADAPTER_R"

    cutadapt \${ARGS} --cores $task.cpus --minimum-length ${params.trim_cutadapt["minimum_length"]} --error-rate ${params.trim_cutadapt["error_rate"]} --nextseq-trim ${params.trim_cutadapt["quality"]} --trim-n --overlap ${params.trim_cutadapt["overlap"]} -o "trimmed/${sampleId}.fastq" "${sampleId}.fastq.gz" > "trimmed/${sampleId}.txt"
    pigz -p $task.cpus -6 "trimmed/${sampleId}.fastq"

    rm -rf work
    """
}

process assign_kraken2 {

    // Taxonomic assignation of reads by kraken2, output filtered reads on a taxon
    // TODO test read extraction included in krakentools, update combine_kreports to keep minimizer info, try RAM disk to reduce database loading time
    // TODO Fix bug of combine_kreports crashing when only one line in report

    executor 'local'
    label 'denovo'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.assign_kraken2.todo == 1

    input:
    tuple(val(sampleId), path(fastq))
    
    output:
    tuple val(sampleId), path("metashotgun/kraken2/${sampleId}{.fastq.gz,_RM.fastq.gz}"), emit: k2fastq
    tuple val(sampleId), path("metashotgun/kraken2/${sampleId}.tsv"), emit: k2report
    
    script:
    memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p work/kraken2 metashotgun/kraken2/

    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8

    #cutoff set to 400. Previously was 100 (up to 06-03-2023)
    if [[ ! -s "${sampleId}.fastq.gz" || \$(stat -L -c%s "${sampleId}.fastq.gz") -le 1000 ]]; then
        touch "metashotgun/kraken2/${sampleId}.fastq.gz" "metashotgun/kraken2/${sampleId}.tsv"
        
        exit 0
    fi
    
    if [[ -s "${params.assign_kraken2["taxid_list"]}" ]]; then
        cat "${params.assign_kraken2["taxid_list"]}" > taxid_list.txt
    else
        for f in \$(grep \$"${params.assign_kraken2["taxa"]}" ${params.assign_bracken["ncbitax2lin"]} | cut -f1); do echo "kraken:taxid|\${f}\$" >> taxid_list.txt; done;
    fi

    if [[ ${params.readtype} == "paired-end" ]]; then
        reformat.sh -Xmx${memory}G threads=$task.cpus in="${sampleId}.fastq.gz" out1="work/${sampleId}_R1.fastq.gz" out2="work/${sampleId}_R2.fastq.gz" int=t qin=auto qout=33
        kraken2 --db ${params.assign_kraken2["db"]} --confidence ${params.assign_kraken2["confidence"]} --minimum-hit-groups ${params.assign_kraken2["minimum-hit-groups"]} --paired --classified-out "work/kraken2/${sampleId}#.fastq" --output "work/kraken2/${sampleId}.out" --report-minimizer-data --report "work/kraken2/${sampleId}.tsv" "work/${sampleId}_R1.fastq.gz" "work/${sampleId}_R2.fastq.gz" --threads $task.cpus
        cat "work/kraken2/${sampleId}_1.fastq" | sed 's/\$/\$/g' | grep -F -f taxid_list.txt -A 3 | sed 's/\\\$//g' | grep -v "^--\$" | pigz -6 -p $task.cpus -c > "work/kraken2/${sampleId}_R1.fastq.gz"
        cat "work/kraken2/${sampleId}_2.fastq" | sed 's/\$/\$/g' | grep -F -f taxid_list.txt -A 3 | sed 's/\\\$//g' | grep -v "^--\$" | pigz -6 -p $task.cpus -c > "work/kraken2/${sampleId}_R2.fastq.gz"
        reformat.sh -Xmx${memory}G threads=$task.cpus in1="work/kraken2/${sampleId}_R1.fastq.gz" in2="work/kraken2/${sampleId}_R2.fastq.gz" out="metashotgun/kraken2/${sampleId}.fastq.gz" qin=auto qout=33
    else
        kraken2 --db ${params.assign_kraken2["db"]} --confidence ${params.assign_kraken2["confidence"]} --minimum-hit-groups ${params.assign_kraken2["minimum-hit-groups"]} --classified-out "work/kraken2/${sampleId}.fastq" --output "work/kraken2/${sampleId}.out" --report-minimizer-data --report "work/kraken2/${sampleId}.tsv" "${sampleId}.fastq.gz" --threads $task.cpus
        cat "work/kraken2/${sampleId}.fastq" | sed 's/\$/\$/g' | grep -F -f taxid_list.txt -A 3 | sed 's/\\\$//g' | grep -v "^--\$" | pigz -6 -p $task.cpus -c > "metashotgun/kraken2/${sampleId}.fastq.gz"
    fi

    if [[ -s "${sampleId}_RM.fastq.gz" ]]; then
        kraken2 --db ${params.assign_kraken2["db"]} --confidence ${params.assign_kraken2["confidence"]} --minimum-hit-groups ${params.assign_kraken2["minimum-hit-groups"]} --classified-out "work/kraken2/${sampleId}_RM.fastq" --output "work/kraken2/${sampleId}_RM.out" --report-minimizer-data --report "work/kraken2/${sampleId}_RM.tsv" "${sampleId}_RM.fastq.gz" --threads $task.cpus
        cat "work/kraken2/${sampleId}_RM.fastq" | sed 's/\$/\$/g' | grep -F -f taxid_list.txt -A 3 | sed 's/\\\$//g' | grep -v "^--\$" | pigz -6 -p $task.cpus -c > "metashotgun/kraken2/${sampleId}_RM.fastq.gz"
    fi

    touch "work/kraken2/${sampleId}_RM.tsv"
    python3 ${params.nfpath}/script/rsv/combine_kreports.py --only-combined --no-headers --report-file "work/kraken2/${sampleId}.tsv" "work/kraken2/${sampleId}_RM.tsv" -o "metashotgun/kraken2/${sampleId}.tsv"

    rm -rf work
    """
}

process assign_tax2lin {

    // Reestimate the abudance of assignments done by kraken2 on a given taxonomy level

    executor 'local'
    label 'denovo'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.assign_tax2lin.todo == 1

    input:
    tuple(val(sampleId), path(assignment))
    
    output:
    tuple val(sampleId), path("metashotgun/table/${sampleId}.tsv"), emit: table
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8
    mkdir -p work metashotgun/table/

    if [[ ! -s ${sampleId}.tsv ]]; then
        echo -e "Eukaryota; Chordata; Mammalia; Primates; Hominidae; Homo; Homo sapiens\\t1" >> "metashotgun/table/${sampleId}.tsv"
        
        exit 0
    fi

    if [[ ${params.assign} == "bracken" ]]; then
        awk -F'\\t' 'FNR > 1 { print \$2"\\t"\$6}' ${sampleId}.tsv > work/${sampleId}.tsv
        LC_ALL=C sort -k 1b,1 work/${sampleId}.tsv > work/${sampleId}_sorted.tsv
        join -a1 -o 2.2,1.2 -t\$'\\t' work/${sampleId}_sorted.tsv ${params.assign_tax2lin["ncbitax2lin"]} > metashotgun/table/${sampleId}.tsv
    elif [[ ${params.assign} == "kraken2" ]]; then
        awk -F'\\t' '{ print \$5"\\t"\$3}' ${sampleId}.tsv > work/${sampleId}.tsv
        LC_ALL=C sort -k 1b,1 work/${sampleId}.tsv > work/${sampleId}_sorted.tsv
        join -a1 -o 2.2,1.2 -t\$'\\t' work/${sampleId}_sorted.tsv ${params.assign_tax2lin["ncbitax2lin"]} > metashotgun/table/${sampleId}.tsv
    fi

    rm -rf work
    """
}

process import_q2table {

    // Import table in qiime2 format

    executor 'local'
    label 'qiime2'
    storeDir (params.result)
    echo false

    when:
    params.import_q2table.todo == 1

    input:
    file 'metashotgun/table/*'
    
    output:
    path("${runprefix}_meta.tsv"), emit: tsv
    path("metashotgun/table.qza"), emit: table
    path("metashotgun/faketax.qza"), emit: faketax
    
    script:
    runprefix = (params.prefix =~ /([^\_]+)_(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export MPLCONFIGDIR=${params.tmp}
    export NUMBA_CACHE_DIR=${params.tmp}
    mkdir -p work metashotgun

    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8
    python2.7 ${params.nfpath}/script/rsv/merge_tables.py metashotgun/table/*.tsv > work/table.tsv

    cat work/table.tsv > ${runprefix}_meta.tsv

    awk 'BEGIN { FS="\\t"; OFS="\\t" ; bp=0 } NR>1 { bp+=1 ; \$1="T"bp } { print \$0 }' ${runprefix}_meta.tsv > work/faketable.tsv
    biom convert -i work/faketable.tsv -o work/table.biom --table-type="OTU table" --to-hdf5
    qiime tools import --input-path work/table.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path metashotgun/table.qza

    awk 'BEGIN { FS="\\t"; OFS="\\t" ; bp=0 } NR>1 { bp+=1 } NR>1 { print "T"bp, \$1 }' ${runprefix}_meta.tsv > work/faketax.tsv
    qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path work/faketax.tsv --output-path metashotgun/faketax.qza

    rm -rf work
    """
}

process filter_srahumanscrubber {

    // Convert paired-end fastq to an interleaved fastq

    label 'srahumanscrubber'
    storeDir params.result
    echo false
    tag { sampleId }

    when:
    params.filter_srahumanscrubber.todo == 1

    input:
    tuple(val(sampleId), path(fastq))
    
    output:
    tuple val(sampleId), path("dehosted/${sampleId}.fastq.gz")
    
    script:
    memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p dehosted/

    if [[ ! -s "${sampleId}.fastq.gz" ]]; then
        touch "dehosted/${sampleId}.fastq.gz"
        
        exit 0
    fi

    if [[ ${params.readtype} == "paired-end" ]]; then
        unpigz -p $task.cpus -c "${sampleId}.fastq.gz" | /opt/scrubber/scripts/scrub.sh -p $task.cpus -x -d ${params.filter_srahumanscrubber["db"]} -o - -i - | repair.sh -Xmx${memory}G threads=$task.cpus in=stdin out=stdout | pigz -6 -p $task.cpus -c > "dehosted/${sampleId}.fastq.gz"
    else
        unpigz -p $task.cpus -c "${sampleId}.fastq.gz" | /opt/scrubber/scripts/scrub.sh -p $task.cpus -x -d ${params.filter_srahumanscrubber["db"]} -o - -i - | pigz -6 -p $task.cpus -c > "dehosted/${sampleId}.fastq.gz"
    fi

    #rm -rf work
    """
}

process downsample_bbmap {

    // Downsample fastq pairs randomly, reproductible with a given seed

    label 'varcall'
    storeDir params.result
    echo false
    tag { sampleId }

    when:
    params.downsample_bbmap.todo == 1

    input:
    tuple(val(sampleId), path(fastq))
    
    output:
    tuple val(sampleId), path("downsampled/${sampleId}.fastq.gz")
    
    script:
    //sampleId = retrieveMetadata(params.metadata, sampleId, "externalid")
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p work downsampled

    if [[ ! -s "${sampleId}.fastq.gz" ]]; then
        touch "downsampled/${sampleId}.fastq.gz"
        
        exit 0
    fi

    if [[ ${params.readtype} == "paired-end" ]]; then
        ARGS="int=t"
    else
        ARGS="int=f"
    fi

    reformat.sh in="${sampleId}.fastq.gz" out="downsampled/${sampleId}.fastq.gz" \${ARGS} samplereadstarget=${params.downsample_bbmap["samplereadstarget"]} sampleseed=100 qin=${params.downsample_bbmap["encoding"]} qout=33

    rm -rf work
    """
}

process premap_bwamem {

    // Map reads on a set of reference, permit to choose an appropriate reference for each sample, and verify sequencing success of positive controls
    // TODO Reformat depth tsv to keep positions as column while permitting merging

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.premap_bwamem.todo == 1

    input:
    tuple(val(sampleId), path(fastq), path("ref/*"), path("posc/*"))
    
    output:
    tuple val(sampleId), path("varcall/premap/depth/${sampleId}.tsv"), emit: depthtable
    tuple val(sampleId), path("varcall/premap/coverage/${sampleId}.tsv"), emit: coveragetable
    path("varcall/premap/posc/${sampleId}.tsv"), emit: poscsummary
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p varcall/premap/depth varcall/premap/coverage varcall/premap/posc work/premap/coverage

    if [[ ${params.readtype} == "paired-end" ]]; then
        ARGS="-p"
    else
        ARGS=""
    fi

    for ref in ref/*; do bn=\$(basename "\$ref" | cut -d. -f1); grep -v "^>" "\${ref}" | tr -d '\\015' | tr -d '\\n' | \
    awk -v ref=\${bn%%.*} 'BEGIN { print ">"ref } { print }' >> work/refpm.fna ; done;
    cat work/refpm.fna >> ref/refpm.fna
    for ref in posc/*; do bn=\$(basename "\$ref" | cut -d. -f1); grep -v "^>" "\${ref}" | tr -d '\\015' | tr -d '\\n' | \
    awk -v ref=\${bn%%.*} 'BEGIN { print ">"ref } { print }' >> ref/refpm.fna ; done;
    for ref in posc/*; do bn=\$(basename "\$ref" | cut -d. -f1); echo \$bn >> posc.txt; done;

    if [[ ! -s "${sampleId}.fastq.gz" || \$(grep -v '>' "work/refpm.fna" | wc -c) -le 100 ]]; then
        touch "varcall/premap/coverage/${sampleId}.cov" "varcall/premap/depth/${sampleId}.tsv" "varcall/premap/coverage/${sampleId}.tsv" "varcall/premap/posc/${sampleId}.tsv"
        
        exit 0
    fi
    
    bwa index ref/refpm.fna
    bwa mem \${ARGS} -A ${params.premap_bwamem["A"]} -B ${params.premap_bwamem["B"]} -O ${params.premap_bwamem["O"]} -E ${params.premap_bwamem["E"]} -L ${params.premap_bwamem["L"]} -U ${params.premap_bwamem["U"]} -T ${params.premap_bwamem["T"]} -t $task.cpus ref/refpm.fna "${sampleId}.fastq.gz" | \
    samtools view --threads $task.cpus -b -F 256 - | \
    samtools sort --threads $task.cpus -o "ref/${sampleId}_sorted.bam"
    samtools index "ref/${sampleId}_sorted.bam"

    bedtools genomecov -ibam "ref/${sampleId}_sorted.bam" -d > "${sampleId}.tsv"
    posc=\$(cat posc.txt | perl -pe 'chomp if eof' | tr '\\n' '\\|'); egrep "\$posc" "${sampleId}.tsv" > "work/${sampleId}_posc_depth.tsv"; egrep -v "\$posc" "${sampleId}.tsv" > "work/${sampleId}_ref_depth.tsv"
    for ref in \$(cut -f 1 "work/${sampleId}_ref_depth.tsv" | sort -u); do cov=\$(bc <<< "scale=4; (\$(awk -v ref=\$ref 'BEGIN {bp=0} \$1==ref && \$3>=${params.premap_bwamem["min_depth"]} {bp+=1} {print bp}' "work/${sampleId}_ref_depth.tsv" | tail -1)/\$(awk -v ref=\$ref 'BEGIN {bp=0} \$1==ref {bp+=1} {print bp}' "work/${sampleId}_ref_depth.tsv" | tail -1))"); dep=\$(bc <<< "scale=4; (\$(awk -v ref=\$ref 'BEGIN {bp=0} \$1==ref && \$3>=${params.premap_bwamem["min_depth"]} {bp+=\$3} {print bp}' "work/${sampleId}_ref_depth.tsv" | tail -1)/\$(awk -v ref=\$ref 'BEGIN {bp=0} \$1==ref && \$3>=${params.premap_bwamem["min_depth"]} {bp+=1} {print bp}' "work/${sampleId}_ref_depth.tsv" | tail -1))"); echo -e "\${ref}\\t\${cov}\\t\${dep}" >> "work/premap/coverage/${sampleId}.tsv"; done;
    sort -nr -t\$'\\t' -k 3 "work/premap/coverage/${sampleId}.tsv" | cut -d\$'\\t' -f 1-2 > "varcall/premap/coverage/${sampleId}.tsv"
    for ref in \$(cut -f 1 "work/${sampleId}_posc_depth.tsv" | sort -u); do cov=\$(bc <<< "scale=4; (\$(awk -v ref=\$ref 'BEGIN {bp=0} \$1==ref && \$3>=${params.premap_bwamem["min_depth"]} {bp+=1} {print bp}' "work/${sampleId}_posc_depth.tsv" | tail -1)/\$(awk -v ref=\$ref 'BEGIN {bp=0} \$1==ref {bp+=1} {print bp}' "work/${sampleId}_posc_depth.tsv" | tail -1))"); echo -e "\${ref}\\t\${cov}" >> "varcall/premap/posc/${sampleId}.tsv"; done;
    awk -F'[\\t]' '{print \$1"_"\$2"\\t"\$3}' "${sampleId}.tsv" > "varcall/premap/depth/${sampleId}.tsv"
    #rm -rf ref work
    """
}

process premap_verif {

    // 

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.premap_verif.todo == 1

    input:
    tuple(val(sampleId), path(fastq), path("ref/*"))
    
    output:
    tuple val(sampleId), path("varcall/premap/verif/${sampleId}.tsv"), emit: veriftable
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p varcall/premap/verif work/ ref/

    if [[ ! -s "${sampleId}.fastq.gz" ]]; then
        touch "varcall/premap/verif/${sampleId}.tsv"
        
        exit 0
    fi

    if [[ ${params.readtype} == "paired-end" ]]; then
        BOWTIE2_INPUT="--no-mixed --no-discordant --interleaved ${sampleId}.fastq.gz"
    else
        BOWTIE2_INPUT="-U ${sampleId}.fastq.gz"
    fi

    for ref in ref/*; do cat \$ref >> ref/refpm.fna ; done;
    bowtie2-build -f ref/refpm.fna ref/premap_index
    bowtie2 -q -x ref/premap_index -p $task.cpus ${params.premap_verif["args"]} --score-min ${params.premap_verif["score-min"]} --dpad ${params.premap_verif["dpad"]} --ma ${params.premap_verif["ma"]} --mp ${params.premap_verif["mp"]} --rdg ${params.premap_verif["rdg"]} --rfg ${params.premap_verif["rfg"]} --gbar ${params.premap_verif["gbar"]} \${BOWTIE2_INPUT} | \
    samtools view --threads $task.cpus -b -F 256 - | \
    samtools sort --threads $task.cpus -o "ref/${sampleId}_sorted.bam"
    samtools index "ref/${sampleId}_sorted.bam"
    bedtools genomecov -ibam "ref/${sampleId}_sorted.bam" -d > "work/${sampleId}.tsv"
    for ref in \$(cut -f 1 "work/${sampleId}.tsv" | sort -u); do cov=\$(bc <<< "scale=4; (\$(awk -v ref=\$ref 'BEGIN {bp=0} \$1==ref && \$3>=${params.premap_verif["min_depth"]} {bp+=1} {print bp}' "work/${sampleId}.tsv" | tail -1)/\$(awk -v ref=\$ref 'BEGIN {bp=0} \$1==ref {bp+=1} {print bp}' "work/${sampleId}.tsv" | tail -1))"); echo -e "\${ref}\\t\${cov}" >> "varcall/premap/verif/${sampleId}.tsv"; done;
    rm -rf work
    """
}

process map_bwamem {

    // Map reads on each reference passing a given threshold, multiple reference can be used if max_match > 1

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.map_bwamem.todo == 1

    input:
    tuple(val(sampleId), path(fastq), path("coverage/${sampleId}.tsv"), path(ref))
    
    output:
    tuple val(refName), val(sampleId), path("varcall/bwamem/${ref.simpleName}/${sampleId}.bam"), path("varcall/bwamem/${ref.simpleName}/${sampleId}.bam.bai"), path("varcall/ref/${ref.simpleName}.fna"), emit: bam
    tuple val(refName), val(sampleId), path("varcall/bwamem/${ref.simpleName}/${sampleId}_tomapreadcount.tsv"), emit: fastqstat
    
    script:
    refName = (ref =~  /([^\.]+)\.(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p work/fastqc varcall/ref/ varcall/bwamem/${ref.simpleName}/
    
    bn=\$(basename "${ref}" | cut -d. -f1)
    awk -F '[\\t]' '{print \$1"\\t"\$NF}' "coverage/${sampleId}.tsv" > "coverage/premap.tsv"

    if [[ ! -s "${sampleId}.fastq.gz" || ! \$bn =~ ^(\$(cat "coverage/premap.tsv" | sort -t\$'\\t' -k2 -nr | head -n ${params.map_bwamem["max_match"]} | awk -F'[\\t]' '{ORS="|"; if (\$2>${params.map_bwamem["min_cov"]}) print \$1 }'))\$ ]]; then
        echo -e "${ref.simpleName}\\tfastq_tomapreadcount\\t0\\t${sampleId}" > "varcall/bwamem/${ref.simpleName}/${sampleId}_tomapreadcount.tsv"
        touch "varcall/bwamem/${ref.simpleName}/${sampleId}.bam" "varcall/bwamem/${ref.simpleName}/${sampleId}.bam.bai"
        cp "${ref}" "varcall/ref/${ref.simpleName}.fna"

        exit 0
    fi

    if [[ ${params.readtype} == "paired-end" ]]; then
        ARGS="-p"
        READ_FLAG="12"
    else
        ARGS=""
        READ_FLAG="4"
    fi

    bwa index ${ref}
    bwa mem \${ARGS} -A ${params.map_bwamem["A"]} -B ${params.map_bwamem["B"]} -O ${params.map_bwamem["O"]} -E ${params.map_bwamem["E"]} -L ${params.map_bwamem["L"]} -U ${params.map_bwamem["U"]} -T ${params.map_bwamem["T"]} -t $task.cpus ${ref} "${sampleId}.fastq.gz" | \
    samtools view --threads $task.cpus -b -F \${READ_FLAG} - | \
    samtools sort --threads $task.cpus -o "varcall/bwamem/${ref.simpleName}/${sampleId}.bam"
    samtools index "varcall/bwamem/${ref.simpleName}/${sampleId}.bam"

    fastqc "${sampleId}.fastq.gz" --outdir work/fastqc -t $task.cpus
    unzip -d work "work/fastqc/${sampleId}_fastqc.zip"
    cat \$(find work/ -type f -name fastqc_data.txt) >> "work/${sampleId}.txt"
    echo -e "${ref.simpleName}\\tfastq_tomapreadcount\\t\$(grep "Total Sequences" "work/${sampleId}.txt" | cut -f 2 | paste -sd+ | bc)\\t${sampleId}" > "varcall/bwamem/${ref.simpleName}/${sampleId}_tomapreadcount.tsv"

    cp "${ref}" "varcall/ref/${ref.simpleName}.fna"
    rm -rf work
    """
}

process postmap_picard {

    // Filter any duplicate reads as determined by Picard
    // TODO Remove picard version in path

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.postmap_picard.todo == 1

    input:
    tuple(val(refName), val(sampleId), path(bam), path(bai), path(ref))
    
    output:
    tuple val(refName), val(sampleId), path("varcall/picard/${ref.simpleName}/${sampleId}.bam"), path("varcall/picard/${ref.simpleName}/${sampleId}.bam.bai"), path("varcall/ref/${ref.simpleName}.fna")
    
    script:
    memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p work/${ref.simpleName}/ varcall/ref/ varcall/picard/${ref.simpleName}/

    if [[ ! -s ${bam} ]]; then
        touch "varcall/picard/${ref.simpleName}/${sampleId}.bam" "varcall/picard/${ref.simpleName}/${sampleId}.bam.bai"
        cp "${ref}" "varcall/ref/${ref.simpleName}.fna"

        exit 0
    fi

    java -XX:ParallelGCThreads=$task.cpus -Xmx${memory}G -jar /opt/conda/envs/varcall/bin/picard.jar MarkDuplicates I=${bam} O=varcall/picard/${ref.simpleName}/${sampleId}.bam M=marked_dup_metrics.txt REMOVE_DUPLICATES=true TMP_DIR=${params.tmp}
    samtools index "varcall/picard/${ref.simpleName}/${sampleId}.bam"

    cp "${ref}" "varcall/ref/${ref.simpleName}.fna"
    rm -rf work
    """
}

process capdepth_jvarkit {

    // Downsample a bam file to limit maxdepth by position
    // TODO Verify if downsampling is not biased as it is not random

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.capdepth_jvarkit.todo == 1

    input:
    tuple(val(refName), val(sampleId), path(bam), path(bai), path(ref))
    
    output:
    tuple val(refName), val(sampleId), path("varcall/jvarkit/${ref.simpleName}/${sampleId}.bam"), path("varcall/jvarkit/${ref.simpleName}/${sampleId}.bam.bai"), path("varcall/ref/${ref.simpleName}.fna")
    
    script:
    memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p work/${ref.simpleName}/ varcall/ref/ varcall/jvarkit/${ref.simpleName}/

    if [[ ! -s ${bam} ]]; then
        touch "varcall/jvarkit/${ref.simpleName}/${sampleId}.bam" "varcall/jvarkit/${ref.simpleName}/${sampleId}.bam.bai"
        cp "${ref}" "varcall/ref/${ref.simpleName}.fna"

        exit 0
    fi

    java -XX:ParallelGCThreads=$task.cpus -Xmx${memory}G -jar /opt/conda/envs/varcall/bin/picard.jar CreateSequenceDictionary R=${ref} TMP_DIR=${params.tmp}
    samtools faidx ${ref}

    java -XX:ParallelGCThreads=$task.cpus -Xmx${memory}G -jar /opt/conda/envs/varcall/bin/sortsamrefname.jar --reference ${ref} --samoutputformat SAM ${bam} | \
    java -XX:ParallelGCThreads=$task.cpus -Xmx${memory}G -jar /opt/conda/envs/varcall/bin/biostar154220.jar --depth ${params.capdepth_jvarkit["max_depth"]} --samoutputformat BAM | \
    samtools sort --threads $task.cpus -o "varcall/jvarkit/${ref.simpleName}/${sampleId}.bam"
    samtools index "varcall/jvarkit/${ref.simpleName}/${sampleId}.bam"

    cp "${ref}" "varcall/ref/${ref.simpleName}.fna"
    rm -rf work
    """
}

process postmap_abra2 {

    // Realign soft-clipped reads, trying to align all reads of a given region in the same way, permit to improve detection of indels 
    // TODO Adapt for ONT, remove abra2 version in path

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.postmap_abra2.todo == 1

    input:
    tuple(val(refName), val(sampleId), path(bam), path(bai), path(ref))
    
    output:
    tuple val(refName), val(sampleId), path("varcall/abra2/${ref.simpleName}/${sampleId}.bam"), path("varcall/abra2/${ref.simpleName}/${sampleId}.bam.bai"), path("varcall/ref/${ref.simpleName}.fna")
    
    script:
    memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p work varcall/ref/ varcall/abra2/${ref.simpleName}/

    if [[ ! -s ${bam} ]]; then
        touch "varcall/abra2/${ref.simpleName}/${sampleId}.bam" "varcall/abra2/${ref.simpleName}/${sampleId}.bam.bai"
        cp "${ref}" "varcall/ref/${ref.simpleName}.fna"

        exit 0
    fi

    samtools view ${bam} | awk 'BEGIN {FS="\\t"; OFS="\\t"} \$6 ~ /^[0-9][M]([0-9]+)[DI]/ || \$6 ~ /([0-9]+)[DI][0-9][M]\$/ {print \$1}' > "work/toexclude.txt"
    samtools view -h ${bam} | grep -v -F -f "work/toexclude.txt" | samtools view -b | samtools sort > "work/${sampleId}.bam"
    samtools index "work/${sampleId}.bam"

    java -XX:ParallelGCThreads=$task.cpus -Xmx${memory}G -jar /opt/conda/envs/varcall/bin/abra2.jar ${params.postmap_abra2["args"]} --amq ${params.postmap_abra2["amq"]} --ca ${params.postmap_abra2["ca"]} --mac ${params.postmap_abra2["mac"]} --mad ${params.postmap_abra2["mad"]} --mapq ${params.postmap_abra2["mapq"]} --maxn ${params.postmap_abra2["maxn"]} --mbq ${params.postmap_abra2["mbq"]} --mcl ${params.postmap_abra2["mcl"]} --mcr ${params.postmap_abra2["mcr"]} --mer ${params.postmap_abra2["mer"]} --mmr ${params.postmap_abra2["mmr"]} --mnf ${params.postmap_abra2["mnf"]} --mrn ${params.postmap_abra2["mrn"]} --mrr ${params.postmap_abra2["mrr"]} --msr ${params.postmap_abra2["msr"]} --rcf ${params.postmap_abra2["rcf"]} --in "work/${sampleId}.bam" --out "varcall/abra2/${ref.simpleName}/${sampleId}.bam" --ref ${ref} --threads $task.cpus --tmpdir ${params.tmp}
    samtools index "varcall/abra2/${ref.simpleName}/${sampleId}.bam"

    cp "${ref}" "varcall/ref/${ref.simpleName}.fna"
    rm -rf work
    """
}

process postmap_ampliconclip {

    // Hard-clip all reads in regions specified in a given bed file
    // TODO

    label 'samtools'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.postmap_ampliconclip.todo == 1

    input:
    tuple(val(refName), val(sampleId), path(bam), path(bai), path(ref), path(bed))
    
    output:
    tuple val(refName), val(sampleId), path("varcall/ampliconclip/${ref.simpleName}/${sampleId}.bam"), path("varcall/ampliconclip/${ref.simpleName}/${sampleId}.bam.bai"), path("varcall/ref/${ref.simpleName}.fna")
    
    script:
    memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p work/${ref.simpleName}/ varcall/ref/ varcall/ampliconclip/${ref.simpleName}/

    if [[ ! -s ${bam} ]]; then
        touch "varcall/ampliconclip/${ref.simpleName}/${sampleId}.bam" "varcall/ampliconclip/${ref.simpleName}/${sampleId}.bam.bai"
        cp "${ref}" "varcall/ref/${ref.simpleName}.fna"

        exit 0
    fi
    
    samtools ampliconclip -o "work/${ref.simpleName}/${sampleId}.bam" --hard-clip --both-ends --filter-len ${params.postmap_ampliconclip["filter_len"]} -b ${bed} ${bam}
    samtools sort --threads $task.cpus -o "varcall/ampliconclip/${ref.simpleName}/${sampleId}.bam" "work/${ref.simpleName}/${sampleId}.bam"
    samtools index "varcall/ampliconclip/${ref.simpleName}/${sampleId}.bam"

    cp "${ref}" "varcall/ref/${ref.simpleName}.fna"
    rm -rf work
    """
}

process postmap_fixbam {

    // Fix bam on various criteria, remove reads pair with a read beginning or ending with an indel, remove soft-clipped reads
    // TODO verify impact soft-clipping to hard-clipping

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.postmap_fixbam.todo == 1

    input:
    tuple(val(refName), val(sampleId), path(bam), path(bai), path(ref))
    
    output:
    tuple val(refName), val(sampleId), path("varcall/fixbam/${ref.simpleName}/${sampleId}.bam"), path("varcall/fixbam/${ref.simpleName}/${sampleId}.bam.bai"), path("varcall/ref/${ref.simpleName}.fna"), emit: bam
    tuple val(refName), val(sampleId), path("varcall/fixbam/${ref.simpleName}/${sampleId}_softclippedreadcount.tsv"), emit: bamstat
    
    script:
    memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p work/${ref.simpleName}/ varcall/ref/ varcall/fixbam/${ref.simpleName}/

    if [[ ! -s ${bam} ]]; then
        echo -e "${refName}\\tbam_softclippedreadcount\\t0\\t${sampleId}" > "varcall/fixbam/${ref.simpleName}/${sampleId}_softclippedreadcount.tsv"
        touch "varcall/fixbam/${ref.simpleName}/${sampleId}.bam" "varcall/fixbam/${ref.simpleName}/${sampleId}.bam.bai"
        cp "${ref}" "varcall/ref/${ref.simpleName}.fna"

        exit 0
    fi
    
    echo -e "${refName}\\tbam_softclippedreadcount\\t\$(samtools view -F 260 ${bam} | awk 'BEGIN {FS="\\t"; OFS="\\t"} \$6 ~ /S/ {print \$0}' | wc -l)\\t${sampleId}" > "varcall/fixbam/${ref.simpleName}/${sampleId}_softclippedreadcount.tsv"
    #samtools view ${bam} | awk 'BEGIN {FS="\\t"; OFS="\\t"} \$6 ~ /^([0-9]+)[DI]/ || \$6 ~ /[DI]\$/ || \$6 ~ /S/ {print \$1}' > "work/toexclude.txt"
    samtools view ${bam} | awk 'BEGIN {FS="\\t"; OFS="\\t"} \$6 ~ /^([0-9]+)[DI]/ || \$6 ~ /[DI]\$/ || \$6 ~ /^[0-9][M]([0-9]+)[DI]/ || \$6 ~ /([0-9]+)[DI][0-9][M]\$/ {print \$1}' > "work/toexclude.txt"
    samtools view -h ${bam} | grep -v -F -f "work/toexclude.txt" | java -XX:ParallelGCThreads=$task.cpus -Xmx${memory}G -jar /opt/conda/envs/varcall/bin/biostar84452.jar --samoutputformat bam | samtools sort > "varcall/fixbam/${ref.simpleName}/${sampleId}.bam"
    samtools index "varcall/fixbam/${ref.simpleName}/${sampleId}.bam" || exit 1

    cp "${ref}" "varcall/ref/${ref.simpleName}.fna"
    rm -rf work
    """
}

process callvar_freebayes {

    // Call and filter variants based on a given set of quality criteria, done separately for major and minor variants
    // TODO When filtering minor variants, evaluate strand bias appropriately

    label 'freebayes'
    storeDir (params.result)
    echo false
    tag { sampleId }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.callvar_freebayes.todo == 1

    input:
    tuple(val(refName), val(sampleId), path(bam), path(bai), path(ref))
    
    output:
    path "varcall/ref/${ref.simpleName}.fna"
    path "varcall/raw/${ref.simpleName}/${sampleId}.vcf"
    tuple val(refName), val(sampleId), path("varcall/raw/${ref.simpleName}/${sampleId}.vcf"), emit: rawvcf
    tuple val(refName), val(sampleId), path("varcall/ref/${ref.simpleName}.fna"), path("varcall/vcf/${ref.simpleName}/${sampleId}.vcf"), emit: gentsv
    tuple val(refName), val(sampleId), path("varcall/ref/${ref.simpleName}.fna"), path("varcall/vcf/${ref.simpleName}/${sampleId}.vcf"), emit: callcons
    tuple val(refName), val(sampleId), path("varcall/ref/${ref.simpleName}.fna"), path("varcall/vcf/${ref.simpleName}/${sampleId}.vcf"), emit: contaminated
    tuple val(refName), path("varcall/vcf/${ref.simpleName}/${sampleId}.vcf"), emit: contaminant
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}

    mkdir -p work varcall/ref/ varcall/vcf/${ref.simpleName}/ varcall/raw/${ref.simpleName}/ work/${ref.simpleName}/
    perl -i.bkp -pe 's/\\r\$//' ${ref}

    if [[ ! -s ${bam} ]]; then
        touch "varcall/vcf/${ref.simpleName}/${sampleId}.vcf" "varcall/raw/${ref.simpleName}/${sampleId}.vcf"

        cp "${ref}" "varcall/ref/${ref.simpleName}.fna"

        exit 0
    fi

    freebayes --theta ${params.callvar_freebayes["theta"]} --ploidy ${params.callvar_freebayes["ploidy"]} --report-all-haplotype-alleles --pooled-continuous --use-best-n-alleles ${params.callvar_freebayes["use-best-n-alleles"]} --allele-balance-priors-off --haplotype-length ${params.callvar_freebayes["haplotype_length"]} --use-duplicate-reads --genotyping-max-iterations ${params.callvar_freebayes["genotyping-max-iterations"]} --genotyping-max-banddepth ${params.callvar_freebayes["genotyping-max-banddepth"]} --min-mapping-quality ${params.callvar_freebayes["min_mapping_quality"]} --min-base-quality ${params.callvar_freebayes["min_base_quality"]} -F 0.5 -C ${params.callvar_freebayes["min_var_depth"]} --min-coverage ${params.callvar_freebayes["min_depth"]} --f ${ref} -b ${bam} > "work/${ref.simpleName}/${sampleId}_cons.vcf"

    vt decompose -s "work/${ref.simpleName}/${sampleId}_cons.vcf" -o "work/${ref.simpleName}/${sampleId}_decomposed_snps.vcf"
    python2.7 ${params.nfpath}/script/rsv/recalc_fbvcf.py -i "work/${ref.simpleName}/${sampleId}_decomposed_snps.vcf" -o "work/${ref.simpleName}/${sampleId}_recalc_snps.vcf" -n 0.0 -x 1.0 -t snp,del,ins,mnp,complex
    bcftools filter -i 'STRLEN(REF)==STRLEN(ALT) & AF >= 0.5 & ${params.callvar_freebayes["bcftools_filter_snps"]}' "work/${ref.simpleName}/${sampleId}_recalc_snps.vcf" > "work/${ref.simpleName}/${sampleId}_filtered_snps.vcf"
    #bcftools filter -i 'STRLEN(REF)==STRLEN(ALT) & AF >= 0.5 & ((SAF>=${params.callvar_freebayes["saf"]} & SAR>=${params.callvar_freebayes["sar"]}) | AF >= 0.9) & INFO/QA>=${params.callvar_freebayes["qa"]}' "work/${ref.simpleName}/${sampleId}_recalc_snps.vcf" > "work/${ref.simpleName}/${sampleId}_filtered_snps.vcf"
    vt normalize "work/${ref.simpleName}/${sampleId}_filtered_snps.vcf" -r "${ref}" -o "work/${ref.simpleName}/${sampleId}_normalized_snps.vcf"

    bgzip < "work/${ref.simpleName}/${sampleId}_normalized_snps.vcf" > "work/${ref.simpleName}/${sampleId}_snps.vcf.gz"
    bcftools index "work/${ref.simpleName}/${sampleId}_snps.vcf.gz"
    cat ${ref} | bcftools consensus "work/${ref.simpleName}/${sampleId}_snps.vcf.gz" > "work/${ref.simpleName}/${sampleId}_snps.fna"

    vt decompose -s "work/${ref.simpleName}/${sampleId}_cons.vcf" -o "work/${ref.simpleName}/${sampleId}_decomposed_complex.vcf"
    python2.7 ${params.nfpath}/script/rsv/recalc_fbvcf.py -i "work/${ref.simpleName}/${sampleId}_decomposed_complex.vcf" -o "work/${ref.simpleName}/${sampleId}_recalc_complex.vcf" -n 0.0 -x 1.0 -t snp,del,ins,mnp,complex
    bcftools filter -i 'STRLEN(REF)!=STRLEN(ALT) & AF >= 0.5 & ${params.callvar_freebayes["bcftools_filter_complex"]}' "work/${ref.simpleName}/${sampleId}_recalc_complex.vcf" > "work/${ref.simpleName}/${sampleId}_filtered_complex.vcf"
    #bcftools filter -i 'STRLEN(REF)!=STRLEN(ALT) & AF >= 0.5 & ((SAF>=${params.callvar_freebayes["saf"]} & SAR>=${params.callvar_freebayes["sar"]}) | AF >= 0.9) & INFO/QA>=${params.callvar_freebayes["qa"]}' "work/${ref.simpleName}/${sampleId}_recalc_complex.vcf" > "work/${ref.simpleName}/${sampleId}_filtered_complex.vcf"
    vt normalize "work/${ref.simpleName}/${sampleId}_filtered_complex.vcf" -r "${ref}" -o "work/${ref.simpleName}/${sampleId}_normalized_complex.vcf"

    freebayes --theta ${params.callvar_freebayes["theta"]} --ploidy ${params.callvar_freebayes["ploidy"]} --report-all-haplotype-alleles --pooled-continuous --use-best-n-alleles ${params.callvar_freebayes["use-best-n-alleles"]} --allele-balance-priors-off --haplotype-length ${params.callvar_freebayes["haplotype_length"]} --use-duplicate-reads --genotyping-max-iterations ${params.callvar_freebayes["genotyping-max-iterations"]} --genotyping-max-banddepth ${params.callvar_freebayes["genotyping-max-banddepth"]} --min-mapping-quality ${params.callvar_freebayes["min_mapping_quality"]} --min-base-quality ${params.callvar_freebayes["min_base_quality"]} -F 0.01 -C ${params.callvar_freebayes["min_var_depth"]} --min-coverage ${params.callvar_freebayes["min_depth"]} --f "work/${ref.simpleName}/${sampleId}_snps.fna" -b ${bam} > "work/${ref.simpleName}/${sampleId}_var.vcf"

    vt decompose -s "work/${ref.simpleName}/${sampleId}_var.vcf" -o "work/${ref.simpleName}/${sampleId}_decomposed_var.vcf"
    python2.7 ${params.nfpath}/script/rsv/recalc_fbvcf.py -i "work/${ref.simpleName}/${sampleId}_decomposed_var.vcf" -o "work/${ref.simpleName}/${sampleId}_recalc_var.vcf" -n 0.0 -x 0.499 -t snp,del,ins,mnp,complex
    bcftools filter -i 'AF < 0.5 & ${params.callvar_freebayes["bcftools_filter_minor"]}' "work/${ref.simpleName}/${sampleId}_recalc_var.vcf" > "work/${ref.simpleName}/${sampleId}_filtered_var.vcf"
    #bcftools filter -i 'AF < 0.5 & AF >= ${params.callvar_freebayes["min_freq"]} & ( (SAF/(SAF+SAR))>=((SRF/(SRF+SRR))-${params.callvar_freebayes["strandbias_tolerance"]}) & (SAF/(SAF+SAR))<=((SRF/(SRF+SRR))+${params.callvar_freebayes["strandbias_tolerance"]}) & (SAF/(SAF+SAR))>=0.01 & (SAF/(SAF+SAR))<=0.99 ) & SAF>=${params.callvar_freebayes["saf"]} & SAR>=${params.callvar_freebayes["sar"]} & INFO/QA>=${params.callvar_freebayes["qa"]}' "work/${ref.simpleName}/${sampleId}_recalc_var.vcf" > "work/${ref.simpleName}/${sampleId}_filtered_var.vcf"
    vt normalize "work/${ref.simpleName}/${sampleId}_filtered_var.vcf" -r "work/${ref.simpleName}/${sampleId}_snps.fna" -o "work/${ref.simpleName}/${sampleId}_normalized_var.vcf"

    #Prepare final filtered vcf file
    cp "work/${ref.simpleName}/${sampleId}_normalized_snps.vcf" "work/${ref.simpleName}/${sampleId}_unsorted_cons.vcf"
    grep -v '#' "work/${ref.simpleName}/${sampleId}_normalized_complex.vcf" >> "work/${ref.simpleName}/${sampleId}_unsorted_cons.vcf"
    bcftools sort -Ov -o "work/${ref.simpleName}/${sampleId}_sorted_cons.vcf" "work/${ref.simpleName}/${sampleId}_unsorted_cons.vcf"

    cp "work/${ref.simpleName}/${sampleId}_sorted_cons.vcf" "work/${ref.simpleName}/${sampleId}_sorted_concat.vcf"
    grep -v '#' "work/${ref.simpleName}/${sampleId}_normalized_var.vcf" >> "work/${ref.simpleName}/${sampleId}_sorted_concat.vcf"
    bcftools sort -Ov -o "work/${ref.simpleName}/${sampleId}_sorted.vcf" "work/${ref.simpleName}/${sampleId}_sorted_concat.vcf"
    vcfuniq "work/${ref.simpleName}/${sampleId}_sorted.vcf" > "varcall/vcf/${ref.simpleName}/${sampleId}.vcf"

    #Prepare vcf file without any filter post-calling
    vt normalize "work/${ref.simpleName}/${sampleId}_recalc_snps.vcf" -r "${ref}" -o "work/${ref.simpleName}/${sampleId}_normalized_raw_snps.vcf"
    vt normalize "work/${ref.simpleName}/${sampleId}_recalc_complex.vcf" -r "${ref}" -o "work/${ref.simpleName}/${sampleId}_normalized_raw_complex.vcf"
    vt normalize "work/${ref.simpleName}/${sampleId}_recalc_var.vcf" -r "work/${ref.simpleName}/${sampleId}_snps.fna" -o "work/${ref.simpleName}/${sampleId}_normalized_raw_var.vcf"
    cp "work/${ref.simpleName}/${sampleId}_normalized_raw_snps.vcf" "work/${ref.simpleName}/${sampleId}_unsorted_raw_cons.vcf"
    grep -v '#' "work/${ref.simpleName}/${sampleId}_normalized_raw_complex.vcf" >> "work/${ref.simpleName}/${sampleId}_unsorted_raw_cons.vcf"
    bcftools sort -Ov -o "work/${ref.simpleName}/${sampleId}_sorted_raw_cons.vcf" "work/${ref.simpleName}/${sampleId}_unsorted_raw_cons.vcf"

    cp "work/${ref.simpleName}/${sampleId}_sorted_raw_cons.vcf" "work/${ref.simpleName}/${sampleId}_sorted_raw_concat.vcf"
    grep -v '#' "work/${ref.simpleName}/${sampleId}_normalized_raw_var.vcf" >> "work/${ref.simpleName}/${sampleId}_sorted_raw_concat.vcf"
    bcftools sort -Ov -o "work/${ref.simpleName}/${sampleId}_sorted_raw.vcf" "work/${ref.simpleName}/${sampleId}_sorted_raw_concat.vcf"
    vcfuniq "work/${ref.simpleName}/${sampleId}_sorted_raw.vcf" > "varcall/raw/${ref.simpleName}/${sampleId}.vcf"

    cp "${ref}" "varcall/ref/${ref.simpleName}.fna"
    #rm -rf work
    """
}

process callvar_bcftools {

    // Call and filter variants based on a given set of quality criteria, done separately for major and minor variants
    // TODO When filtering minor variants, evaluate strand bias appropriately

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.callvar_bcftools.todo == 1

    input:
    tuple(val(refName), val(sampleId), path(bam), path(bai), path(ref), path(toexclude))
    
    output:
    path "varcall/ref/${ref.simpleName}.fna"
    path "varcall/raw/${ref.simpleName}/${sampleId}.vcf"
    tuple val(refName), val(sampleId), path("varcall/raw/${ref.simpleName}/${sampleId}.vcf"), emit: rawvcf
    tuple val(refName), val(sampleId), path("varcall/ref/${ref.simpleName}.fna"), path("varcall/vcf/${ref.simpleName}/${sampleId}.vcf"), emit: gentsv
    tuple val(refName), val(sampleId), path("varcall/ref/${ref.simpleName}.fna"), path("varcall/vcf/${ref.simpleName}/${sampleId}.vcf"), emit: callcons
    tuple val(refName), val(sampleId), path("varcall/ref/${ref.simpleName}.fna"), path("varcall/vcf/${ref.simpleName}/${sampleId}.vcf"), emit: contaminated
    tuple val(refName), path("varcall/vcf/${ref.simpleName}/${sampleId}.vcf"), emit: contaminant
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}

    mkdir -p work varcall/ref/ varcall/vcf/${ref.simpleName}/ varcall/raw/${ref.simpleName}/ work/${ref.simpleName}/
    perl -i.bkp -pe 's/\\r\$//' ${ref}

    if [[ ! -s ${bam} ]]; then
        touch "varcall/vcf/${ref.simpleName}/${sampleId}.vcf" "varcall/raw/${ref.simpleName}/${sampleId}.vcf"

        cp "${ref}" "varcall/ref/${ref.simpleName}.fna"

        exit 0
    fi

    if [[ ${params.callvar_bcftools["mode"]} == "ont" ]]; then
        ARGS="--no-BAQ"
    else
        ARGS=""
    fi

    bcftools mpileup \${ARGS} -f ${ref} --annotate ${params.callvar_bcftools["annotate"]} -A --max-depth ${params.callvar_bcftools["max-depth"]} --min-MQ ${params.callvar_bcftools["min-MQ"]} --min-BQ ${params.callvar_bcftools["min-BQ"]} --max-BQ ${params.callvar_bcftools["max-BQ"]} --ambig-reads ${params.callvar_bcftools["ambig-reads"]} --ext-prob ${params.callvar_bcftools["ext-prob"]} --gap-frac ${params.callvar_bcftools["gap-frac"]} --indel-bias ${params.callvar_bcftools["indel-bias"]} --indel-size ${params.callvar_bcftools["indel-size"]} --max-idepth ${params.callvar_bcftools["max-idepth"]} --min-ireads ${params.callvar_bcftools["min-ireads"]} --max-read-len ${params.callvar_bcftools["max-read-len"]} --open-prob ${params.callvar_bcftools["open-prob"]} -Ov ${bam} | \
    bcftools filter -i 'SUM(FORMAT/AD[0:1-]) >= 4 & (SUM(FORMAT/AD[0:1-]))/((SUM(FORMAT/AD[0:1-])+(FORMAT/AD[0:0]))) >= 0.01' | \
    bcftools call --threads $task.cpus --keep-alts --multiallelic-caller --ploidy ${params.callvar_bcftools["ploidy"]} --prior ${params.callvar_bcftools["prior"]} -Ob | \
    bcftools annotate --remove INFO/INDEL - | \
    bcftools +fill-tags - -Ob -- -t FORMAT/VAF,FORMAT/VAF1 | \
    bcftools +fill-tags - -Ov -o "work/${ref.simpleName}/${sampleId}.vcf" -- -t 'RAF=1.0-(FORMAT/VAF1)'

    vt decompose -s "work/${ref.simpleName}/${sampleId}.vcf" -o "work/${ref.simpleName}/${sampleId}_decomposed.vcf"
    vt normalize "work/${ref.simpleName}/${sampleId}_decomposed.vcf" -r "${ref}" -o "work/${ref.simpleName}/${sampleId}_normalized.vcf"
    python3 ${params.nfpath}/script/rsv/recalc_bcftoolsvcf.py -v "work/${ref.simpleName}/${sampleId}_normalized.vcf" -o "varcall/raw/${ref.simpleName}/${sampleId}.vcf"
    vcfintersect --intersect-vcf "${toexclude}" --reference "${ref}" "varcall/raw/${ref.simpleName}/${sampleId}.vcf" | grep -v '^#' | cut --delimiter=\$'\t' --fields=1,2,3,4,5,6,7 > "work/${ref.simpleName}/${sampleId}_toexclude.vcf"

    bcftools filter -i 'STRLEN(REF)==STRLEN(ALT) & AF >= 0.5 & ${params.callvar_bcftools["bcftools_filter_snps"]}' "varcall/raw/${ref.simpleName}/${sampleId}.vcf" | grep -vFf "work/${ref.simpleName}/${sampleId}_toexclude.vcf" > "work/${ref.simpleName}/${sampleId}_filtered_snps.vcf"
    bcftools filter -i 'STRLEN(REF)!=STRLEN(ALT) & AF >= 0.5 & ${params.callvar_bcftools["bcftools_filter_complex"]}' "varcall/raw/${ref.simpleName}/${sampleId}.vcf" | grep -vFf "work/${ref.simpleName}/${sampleId}_toexclude.vcf" > "work/${ref.simpleName}/${sampleId}_filtered_complex.vcf"
    bcftools filter -i 'AF < 0.5 & ${params.callvar_bcftools["bcftools_filter_minor"]}' "varcall/raw/${ref.simpleName}/${sampleId}.vcf" > "work/${ref.simpleName}/${sampleId}_filtered_minor.vcf"

    grep '^##' "varcall/raw/${ref.simpleName}/${sampleId}.vcf" > "work/${ref.simpleName}/${sampleId}_sorted_raw_concat.vcf"
    cat "work/${ref.simpleName}/${sampleId}_filtered_snps.vcf" "work/${ref.simpleName}/${sampleId}_filtered_complex.vcf" "work/${ref.simpleName}/${sampleId}_filtered_minor.vcf" | grep '^##bcftools_filterCommand' >> "work/${ref.simpleName}/${sampleId}_sorted_raw_concat.vcf"
    grep '^#CHROM' "work/${ref.simpleName}/${sampleId}_filtered_snps.vcf" >> "work/${ref.simpleName}/${sampleId}_sorted_raw_concat.vcf"
    cat "work/${ref.simpleName}/${sampleId}_filtered_snps.vcf" "work/${ref.simpleName}/${sampleId}_filtered_complex.vcf" "work/${ref.simpleName}/${sampleId}_filtered_minor.vcf" | grep -v '^#' >> "work/${ref.simpleName}/${sampleId}_sorted_raw_concat.vcf"
    bcftools sort -Ov -o "work/${ref.simpleName}/${sampleId}_sorted_raw.vcf" "work/${ref.simpleName}/${sampleId}_sorted_raw_concat.vcf"
    vcfuniq "work/${ref.simpleName}/${sampleId}_sorted_raw.vcf" > "varcall/vcf/${ref.simpleName}/${sampleId}.vcf"

    cp "${ref}" "varcall/ref/${ref.simpleName}.fna"
    #rm -rf work
    """
}

process gentsv_vcf {

    // Convert vcf files to a more user-friendly tsv format

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.gentsv_vcf.todo == 1 && vcf.size() > 0

    input:
    tuple(val(refName), val(sampleId), path(ref), path(vcf))
    
    output:
    tuple val(refName), val(sampleId), path("varcall/ref/${ref.simpleName}.fna"), path("varcall/raw/${refName}/${sampleId}.tsv")
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p varcall/raw/${refName}/ varcall/ref/ work/${refName}/ 

    cp "${ref}" "varcall/ref/${ref.simpleName}.fna"
    vcf2tsv -n "." "${sampleId}.vcf" > "varcall/raw/${refName}/${sampleId}.tsv"
    """
}

process genlist_countsummary {

    // Generate a list of variant counts with multiple thresholds

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.genlist_countsummary.todo == 1

    input:
    tuple(val(refName), val(sampleId), path(ref), path(vcf))
    
    output:
    tuple val(refName), path("varcall/raw/${refName}/${sampleId}_varcount.tsv")
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p varcall/raw/${refName}/ work/${refName}/

    echo -e "${refName}\\tvcf_af05-10count\\t\$(bcftools filter -i 'AF >= 0.05 && AF < 0.10' "${vcf}" | grep -v "^#" | wc -l)\\t${sampleId}" > "varcall/raw/${ref.simpleName}/${sampleId}_varcount.tsv"
    echo -e "${refName}\\tvcf_af10-20count\\t\$(bcftools filter -i 'AF >= 0.10 && AF < 0.20' "${vcf}" | grep -v "^#" | wc -l)\\t${sampleId}" >> "varcall/raw/${ref.simpleName}/${sampleId}_varcount.tsv"
    echo -e "${refName}\\tvcf_af20-50count\\t\$(bcftools filter -i 'AF >= 0.2 && AF < 0.5' "${vcf}" | grep -v "^#" | wc -l)\\t${sampleId}" >> "varcall/raw/${ref.simpleName}/${sampleId}_varcount.tsv"
    echo -e "${refName}\\tvcf_af05-50count\\t\$(bcftools filter -i 'AF >= 0.05 && AF < 0.5' "${vcf}" | grep -v "^#" | wc -l)\\t${sampleId}" >> "varcall/raw/${ref.simpleName}/${sampleId}_varcount.tsv"
    """
}

process genbed_bedtools {

    // Generate coverage files in two different format

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.genbed_bedtools.todo == 1

    input:
    tuple(val(refName), val(sampleId), path(bam), path(bai), path(ref), path(vcf))
    
    output:
    path "varcall/ref/${ref.simpleName}.fna"
    path "varcall/bga/${ref.simpleName}/${sampleId}.bed"
    tuple val(refName), val(sampleId), path("varcall/bga/${ref.simpleName}/${sampleId}.bed"), emit: callcons
    tuple val(refName), val(sampleId), path("varcall/bga/${ref.simpleName}/${sampleId}.bed"), emit: contaminated
    tuple val(refName), path("varcall/bga/${ref.simpleName}/${sampleId}.cov"), emit: covsummary
    tuple val(refName), path("varcall/bga/${ref.simpleName}/${sampleId}_readcount.tsv"), emit: bamstat
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}

    mkdir -p work varcall/ref/ varcall/cons/${ref.simpleName}/ varcall/vcf/${ref.simpleName}/ work/${ref.simpleName}/ varcall/raw/${ref.simpleName}/ varcall/bga/${ref.simpleName}/
    perl -i.bkp -pe 's/\\r\$//' ${ref}

    if [[ ! -s ${bam} ]]; then
        touch "varcall/bga/${ref.simpleName}/${sampleId}.bed"
        echo -e "${refName}\\tbam_readcount\\t0\\t${sampleId}" > "varcall/bga/${ref.simpleName}/${sampleId}_readcount.tsv"

        length=\$(grep -v ">" ${ref} | tr -d '\n' | wc -m)
        #for i in \$(seq 1 \${length}); do echo -e "${ref.simpleName}\t\${i}\t0\t${sampleId}" >> "varcall/bga/${ref.simpleName}/${sampleId}.cov"; done;
        touch "varcall/bga/${ref.simpleName}/${sampleId}.cov"
        for line in \$(grep ">" ${ref} | tr -d '>'); do echo -e "${ref.simpleName}\\tmean_depth_\${line}\\t0\\t${sampleId}" >> "varcall/bga/${ref.simpleName}/${sampleId}.cov"; done;

        cp "${ref}" "varcall/ref/${ref.simpleName}.fna"

        exit 0
    fi

    if [[ "${params.calling}" == "bcftools" ]]; then
        #python3 ${params.nfpath}/script/rsv/gen_vcftogenomecov.py -v ${vcf} -r ${ref} -o "varcall/bga/${ref.simpleName}/${sampleId}.bed" -m bga
        bedtools genomecov -ibam ${bam} -bga > "varcall/bga/${ref.simpleName}/${sampleId}.bed"
    else
        bedtools genomecov -ibam ${bam} -bga > "varcall/bga/${ref.simpleName}/${sampleId}.bed"
    fi

    echo -e "${refName}\\tbam_readcount\\t\$(samtools view -c -F 260 ${bam})\\t${sampleId}" > "varcall/bga/${ref.simpleName}/${sampleId}_readcount.tsv"
    
    for line in \$(samtools coverage ${bam} | tail -n +2 | tr '\\t' ','); do SEG="\${line%%,*}"; echo -e "${ref.simpleName}\\tmean_depth_\${SEG}\\t\$(echo \$line | cut -d ',' -f 7)\\t${sampleId}" >> "varcall/bga/${ref.simpleName}/${sampleId}.cov"; done;
    cp "${ref}" "varcall/ref/${ref.simpleName}.fna"
    #rm -rf work
    """
}

process genfasta_consensus {

    // Generate a consensus sequence, then mask regions with an insufficient depth and following a additional bed file

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.genfasta_consensus.todo == 1

    input:
    tuple(val(refName), val(sampleId), path(ref), path(vcf), path(bga), path(mask))
    
    output:
    tuple val(refName), val(sampleId), path("varcall/cons/${ref.simpleName}/${sampleId}.fna"), emit: consensus

    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}

    mkdir -p work/${ref.simpleName}/ varcall/cons/${ref.simpleName}/
    perl -i.bkp -pe 's/\\r\$//' ${ref}

    if [[ ! -s ${vcf} ]]; then
        echo ">${sampleId}" > "varcall/cons/${ref.simpleName}/${sampleId}.fna"
        seq \$(grep -v ">" ${ref} | tr -d '\n' | wc -m) | sed "c N" | tr -d '\n' | sed -e '\$a\\' >> "varcall/cons/${ref.simpleName}/${sampleId}.fna"

        exit 0
    fi

    bcftools filter -i 'AF >= 0.5 && STRLEN(REF)>STRLEN(ALT)' "${sampleId}.vcf" | grep -v '^#' > "work/${ref.simpleName}/${sampleId}_del.vcf"
    awk 'BEGIN {FS="\\t"; OFS="\\t"} {\$2=\$2--; for (i = 1; i <= length(\$4)-length(\$5); i++) print \$1,\$2+i-1,\$2+i,0 }' "work/${ref.simpleName}/${sampleId}_del.vcf" > "work/${ref.simpleName}/${sampleId}_del.bed"

    awk '\$4 < ${params.genfasta_consensus["min_depth"]}' "${sampleId}.bed" > "work/${ref.simpleName}/${sampleId}_lowcov.bed"
    bedtools subtract -a "work/${ref.simpleName}/${sampleId}_lowcov.bed" -b "work/${ref.simpleName}/${sampleId}_del.bed" > "work/${ref.simpleName}/${sampleId}.bed"

    if [[ -s "${mask}" ]]; then
        cat "${mask}" >> "work/${ref.simpleName}/${sampleId}.bed"
    fi

    awk 'BEGIN {FS="\\t"; OFS="\\t"} { if (\$2-${params.genfasta_consensus["maskwiden"]}<=0) \$2=0; else \$2=\$2-${params.genfasta_consensus["maskwiden"]} } { \$3=\$3+${params.genfasta_consensus["maskwiden"]} } {print \$0}' "work/${ref.simpleName}/${sampleId}.bed" > "work/tmp.bed" && mv "work/tmp.bed" "work/${ref.simpleName}/${sampleId}.bed"

    vcftools --vcf "${sampleId}.vcf" --exclude-bed "work/${ref.simpleName}/${sampleId}.bed" --recode --recode-INFO-all --stdout > "work/${ref.simpleName}/${sampleId}_filtered_cons.vcf"
    bcftools filter -i 'AF >= 0.5' "work/${ref.simpleName}/${sampleId}_filtered_cons.vcf" > "work/${ref.simpleName}/${sampleId}_unsorted_cons.vcf"
    bcftools sort -Ov -o "work/${ref.simpleName}/${sampleId}_sorted_cons.vcf" "work/${ref.simpleName}/${sampleId}_unsorted_cons.vcf"
    bgzip < "work/${ref.simpleName}/${sampleId}_sorted_cons.vcf" > "work/${ref.simpleName}/${sampleId}_cons.vcf.gz"
    bcftools index "work/${ref.simpleName}/${sampleId}_cons.vcf.gz"
    bedtools maskfasta -fi ${ref} -bed "work/${ref.simpleName}/${sampleId}.bed" -fo "work/${ref.simpleName}/${sampleId}_cons_masked.fna" -soft
    cat "work/${ref.simpleName}/${sampleId}_cons_masked.fna" | bcftools consensus "work/${ref.simpleName}/${sampleId}_cons.vcf.gz" > "work/${ref.simpleName}/${sampleId}_cons.fna"
    sed 's/\\(>[^|]*\\)/>${sampleId}/' "work/${ref.simpleName}/${sampleId}_cons.fna" > "work/${ref.simpleName}/${sampleId}_rehead.fna"
    python2.7 ${params.nfpath}/script/rsv/fasta_masklowercase.py -i "work/${ref.simpleName}/${sampleId}_rehead.fna" -o "varcall/cons/${ref.simpleName}/${sampleId}.fna"

    #rm -rf work
    """
}

process gentsv_blast {

    // Verify consensus result by blasting it

    label 'denovo'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.gentsv_blast.todo == 1

    input:
    tuple(val(refName), val(sampleId), path(cons))
    
    output:
    tuple val(refName), val(sampleId), path("varcall/blast/${refName}/${sampleId}.tsv"), emit: blast
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p work varcall/blast/${refName}/
    
    if [[ ! -s "${cons}" || \$(stat -L -c%s "${cons}") -le 100 ]]; then
        touch "varcall/blast/${refName}/${sampleId}.tsv"
        
        exit 0
    fi

    blastn -num_threads $task.cpus -query "${cons}" -task megablast -db ${params.gentsv_blast["blastdb"]} -perc_identity ${params.gentsv_blast["perc_identity"]} -max_target_seqs ${params.gentsv_blast["max_target_seqs"]} -outfmt "${params.gentsv_blast["outfmt"]}" >> "work/${sampleId}.tsv"
    echo -e "QUERY_SEQID\\tQUERY_LENGTH\\tQUERY_ALIGNMENT_START\\tQUERY_ALIGNMENT_END\\tSUBJECT_LENGTH\\tSUBJECT_ALIGNMENT_START\\tSUBJECT_ALIGNMENT_END\\tSUBJECT_TITLE\\tSUBJECT_ACCID\\tSUBJECT_TAXID\\tEVALUE\\tBITSCORE\\tPERC_IDENT\\tNUM_IDENT" >> "varcall/blast/${refName}/${sampleId}.tsv"
    python2.7 ${params.nfpath}/script/rsv/blast_getcustomtitle.py -i "work/${sampleId}.tsv" -o "varcall/blast/${refName}/${sampleId}.tsv" -d "${params.gentsv_blast["accid2title"]}"
    #awk -F'\\t' '{("grep ^" \$9 "\\t ${params.gentsv_blast["accid2title"]} | cut -f2") | getline id ; print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6"\\t"\$7"\\t"id"\\t"\$9"\\t"\$10"\\t"\$11"\\t"\$12"\\t"\$13"\\t"\$14}' "work/${sampleId}.tsv" | tr \$'%' 'P' >> "varcall/blast/${refName}/${sampleId}.tsv"
    tail -n +2 "varcall/blast/${refName}/${sampleId}.tsv" | sort -t\$'\\t' -k${params.gentsv_blast["sorting_field"]} -gr | head -n 1 > varcall/blast/${refName}/${sampleId}.txt

    #rm -rf work
    """
}

process gentable_variant {

    // Generate a table of variants in a consistent and comparable way

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.gentable_variant.todo == 1 && vcf.size() > 0

    input:
    tuple(val(refName), val(sampleId), path(ref), path(vcf), path(gff))
    
    output:
    tuple val(refName), path("varcall/tsv/${refName}/${sampleId}.tsv"), emit: table
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p varcall/tsv/${refName}/ work/${refName}/

    python2.7 ${params.nfpath}/script/rsv/vartsv_table.py -r "${ref}" -g "${gff}" -v "${vcf}" -o "work/${refName}/${sampleId}.tsv" -d ${params.gentable_variant["datatype"]}

    awk -F'[\\t]' '!seen[\$1]++' "work/${refName}/${sampleId}.tsv" > "varcall/tsv/${refName}/${sampleId}.tsv"
    """
}

process merge_variant {

    // Merge variant tables generated for the same reference

    label 'varcall'
    storeDir params.result
    echo false
    tag { refName }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.merge_variant.todo == 1

    input:
    tuple(val(refName), path("*"), val(level))
    
    output:
    tuple val(refName), path("${refName}_table${criteria}.tsv")
    
    script:
    if( level == 0 ) {
    criteria = ''
    }
    else {
    criteria = '_' + level
    }
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8

    if [[ ${level} != 0 ]]; then for it in *.tsv; do sampleId="\${it%%.tsv}" ; cat "\${it}" | while read line ; do echo -n \${line} | cut -f 1  | cut -d '|' -f 1-${level} | tr -d \$'\n' >> "\${sampleId}.txt"; echo -n -e '\\t' >> "\${sampleId}.txt"; echo "\${line}" | cut -f 2 >> "\${sampleId}.txt"; done; done; else for it in *.tsv; do sampleId="\${it%%.tsv}"; cat "\${it}" >> "\${sampleId}.txt"; done; fi;
    python2.7 ${params.nfpath}/script/rsv/merge_tables.py *.txt > ${refName}_table${criteria}.tsv
    """
}

process seek_contaminant {

    // Compare vcf files to count minor variants in a potentially contaminated sample found as consensus in potentially contaminant samples 

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.seek_coinf.todo == 1 && target.size() > 0

    input:
    tuple(val(refName), val(sampleId), path(ref), path(target), path(bga), path("varcall/vcf/${refName}/*"))
    
    output:
    tuple val(refName), path("varcall/conta/${refName}/${sampleId}.tsv"), emit: table
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p varcall/conta/${refName}

    for file in varcall/vcf/${refName}/*.vcf; do
    bn=\$(basename "\${file}")
    if [[ \${bn} != ${target} && -s \${file} ]]; then
        python3 ${params.nfpath}/script/rsv/compare_vcf.py --reference \${file} --af_field ${params.seek_contaminant["af_field"]} --target ${target} --depth ${bga} --min_freq ${params.seek_contaminant["min_freq"]} --mode "raw" --output "varcall/conta/${refName}/${sampleId}.tsv"
    fi
    done

    touch varcall/conta/${refName}/${sampleId}.tsv
    """
}

process seek_coinf {

    // Compare vcf files to count minor variants in a sample found as consensus in any lineage

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.seek_coinf.todo == 1 && target.size() > 0

    input:
    tuple(val(refName), val(sampleId), path(ref), path(target), path(bga), val(vcf))
    
    output:
    tuple val(refName), path("varcall/coinf/${refName}/${sampleId}_01.tsv"), emit: coinf01
    tuple val(refName), path("varcall/coinf/${refName}/${sampleId}_02.tsv"), emit: coinf02
    tuple val(refName), path("varcall/raw/${refName}/${sampleId}_coinf.tsv"), emit: raw
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p work varcall/coinf/${refName} varcall/contadb/${refName} varcall/raw/${refName}
    tbn=\$(basename "${target}")

    for file in ${vcf}*.vcf; do
    bn=\$(basename "\${file}")
    if [[ \${bn} != \${tbn} && -s \${file} ]]; then
        python3 ${params.nfpath}/script/rsv/compare_vcf.py --reference \${file} --af_field ${params.seek_coinf["af_field"]} --target ${target} --depth ${bga} --min_depth ${params.seek_coinf["min_depth"]} --min_freq ${params.seek_coinf["min_freq"]} --mode "maj" --output "varcall/coinf/${refName}/\${tbn%%.*}_01.tsv"
    fi
    done

    echo -e "${refName}\\tvcf_coinf01match\\t\$(cat "varcall/coinf/${refName}/\${tbn%%.*}_01.tsv" | sort -t\$'\\t' -k7,7nr -k2,2n | awk -F "\\t" '{print \$1}' | head -n 1 | tr -d '\\n')\\t${sampleId}" > "varcall/raw/${refName}/${sampleId}_coinf.tsv"
    echo -e "${refName}\\tvcf_coinf01count\\t\$(cat "varcall/coinf/${refName}/\${tbn%%.*}_01.tsv" | sort -t\$'\\t' -k7,7nr -k2,2n | awk -F "\\t" '{print \$5}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_coinf.tsv"
    echo -e "${refName}\\tvcf_coinf01score\\t\$(cat "varcall/coinf/${refName}/\${tbn%%.*}_01.tsv" | sort -t\$'\\t' -k7,7nr -k2,2n | awk -F "\\t" '{print \$7}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_coinf.tsv"
    echo -e "${refName}\\tvcf_coinf01median\\t\$(cat "varcall/coinf/${refName}/\${tbn%%.*}_01.tsv" | sort -t\$'\\t' -k7,7nr -k2,2n | awk -F "\\t" '{print \$8}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_coinf.tsv"
    echo -e "${refName}\\tvcf_coinf01iqr\\t\$(cat "varcall/coinf/${refName}/\${tbn%%.*}_01.tsv" | sort -t\$'\\t' -k7,7nr -k2,2n | awk -F "\\t" '{print \$9}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_coinf.tsv"

    MAJ_MATCH=\$(cat "varcall/raw/${refName}/${sampleId}_coinf.tsv" | grep 'vcf_coinf01match' | awk -F "\\t" '{print \$3}')

    for file in ${vcf}*.vcf; do
    bn=\$(basename "\${file}")
    if [[ \${bn} != \${tbn} && -s \${file} ]]; then
        python3 ${params.nfpath}/script/rsv/compare_vcf.py --reference \${file} --af_field ${params.seek_coinf["af_field"]} --target ${target} --exclusion "${vcf}\${MAJ_MATCH}.vcf" --depth ${bga} --min_depth ${params.seek_coinf["min_depth"]} --min_freq ${params.seek_coinf["min_freq"]} --mode "coinf" --output "varcall/coinf/${refName}/\${tbn%%.*}_02.tsv"
    elif [[ \${MAJ_MATCH} == "NA" && -s \${file} ]]; then
        python3 ${params.nfpath}/script/rsv/compare_vcf.py --reference \${file} --af_field ${params.seek_coinf["af_field"]} --target ${target} --depth ${bga} --min_depth ${params.seek_coinf["min_depth"]} --min_freq ${params.seek_coinf["min_freq"]} --mode "coinf" --output "varcall/coinf/${refName}/\${tbn%%.*}_02.tsv"
    fi
    done

    echo -e "${refName}\\tvcf_coinf02match\\t\$(cat "varcall/coinf/${refName}/\${tbn%%.*}_02.tsv" | sort -t\$'\\t' -k7,7nr -k2,2n | awk -F "\\t" '{print \$1}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_coinf.tsv"
    echo -e "${refName}\\tvcf_coinf02count\\t\$(cat "varcall/coinf/${refName}/\${tbn%%.*}_02.tsv" | sort -t\$'\\t' -k7,7nr -k2,2n | awk -F "\\t" '{print \$5}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_coinf.tsv"
    echo -e "${refName}\\tvcf_coinf02score\\t\$(cat "varcall/coinf/${refName}/\${tbn%%.*}_02.tsv" | sort -t\$'\\t' -k7,7nr -k2,2n | awk -F "\\t" '{print \$7}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_coinf.tsv"
    echo -e "${refName}\\tvcf_coinf02median\\t\$(cat "varcall/coinf/${refName}/\${tbn%%.*}_02.tsv" | sort -t\$'\\t' -k7,7nr -k2,2n | awk -F "\\t" '{print \$8}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_coinf.tsv"
    echo -e "${refName}\\tvcf_coinf02iqr\\t\$(cat "varcall/coinf/${refName}/\${tbn%%.*}_02.tsv" | sort -t\$'\\t' -k7,7nr -k2,2n | awk -F "\\t" '{print \$9}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_coinf.tsv"

    touch  varcall/coinf/${refName}/${sampleId}_01.tsv varcall/coinf/${refName}/${sampleId}_02.tsv varcall/raw/${refName}/${sampleId}_coinf.tsv

    rm -rf work
    """
}

process merge_contaminant {

    // Merge contamination tables generated for the same reference

    label 'varcall'
    storeDir params.result
    echo false
    tag { refName }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.merge_contaminant.todo == 1

    input:
    tuple(val(refName), path("varcall/conta/${refName}/*"))
    
    output:
    tuple val(refName), path("varcall/${refName}_conta.tsv"), emit: merged
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8

    for f in varcall/conta/${refName}/*; do mv \$f "\${f%%.tsv}.tsv"; done;
    if [[ \$( (ls -1 "varcall/conta/${refName}/"* 2>/dev/null || true) | wc -l) -gt 1 ]]; then
        python2.7 ${params.nfpath}/script/rsv/merge_tables.py varcall/conta/${refName}/*.tsv > varcall/${refName}_conta.tsv
    else
        touch varcall/${refName}_conta.tsv
    fi
    """
}

process merge_coverage {

    // Merge premapping coverage files generated for each sample

    label 'varcall'
    storeDir params.result
    echo false
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.merge_coverage.todo == 1

    input:
    file "varcall/premap/coverage/*"
    
    output:
    file "${runprefix}_premap_coverage.tsv"
    
    script:
    runprefix = (params.prefix =~ /([^\_]+)_(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8

    python2.7 ${params.nfpath}/script/rsv/merge_tables.py varcall/premap/coverage/*.tsv > ${runprefix}_premap_coverage.tsv
    """
}

process merge_consensus {

    // Merge consensus files generated for the same reference

    label 'varcall'
    storeDir params.result
    echo false
    tag { refName }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.merge_consensus.todo == 1

    input:
    tuple(val(refName), path("varcall/cons/${refName}/*"))
    
    output:
    tuple val(refName), path("varcall/${refName}_consensus.fna")
    
    script:
    //refId = retrieveMetadata(params.merge_consensus["accid2refid"], refName, "refId")
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8

    cat varcall/cons/${refName}/*.fna > varcall/${refName}_consensus.fna
    """
}

process merge_blast {

    // Merge blast matches generated for each sample

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { refName }

    when:
    params.merge_blast.todo == 1

    input:
    tuple(val(refName), path("varcall/blast/${refName}/*"))
    
    output:
    tuple val(refName), path("varcall/${refName}_blast.tsv")
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8

    awk -f ${params.nfpath}/script/rsv/merge_tsv.awk varcall/blast/${refName}/* > "varcall/${refName}_blast.tsv"
    """
}

process gentsv_nextclade {

    // Generate a nextclade report for a given fasta file

    label 'nextclade'
    storeDir params.result
    echo false
    tag { refName }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.gentsv_nextclade.todo == 1

    input:
    tuple(val(refName), path(consensus), path(dataset))
    
    output:
    tuple val(refName), path("${refName}_nextclade.tsv")
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8
    mkdir -p work

    if [[ ${params.workflow} == "denovo" && ${params.gentsv_nextclade["region"]} == "" ]]; then
        echo "in progress"
        cat ${consensus} | perl -i.bkp -pe 's/\\r\$//' | awk '/^>/ {printf("%s%s\\t",(N>0?"\\n":""),\$0);N++;next;} {printf("%s",\$0);} END {printf("\\n");}' | tr "\\t" "\\n" | grep -v "^--\$" > consensus.fasta
    else
        cat ${consensus} | perl -i.bkp -pe 's/\\r\$//' | awk '/^>/ {printf("%s%s\\t",(N>0?"\\n":""),\$0);N++;next;} {printf("%s",\$0);} END {printf("\\n");}' | tr "\\t" "\\n" | grep '|${params.gentsv_nextclade["region"]}' -A 1 | sed 's/|${params.gentsv_nextclade["region"]}//' | grep -v "^--\$" > consensus.fasta
    fi
    
    if [[ ! -s consensus.fasta ]]; then
        touch "${refName}_nextclade.tsv"

        exit 0
    fi

    nextclade run --input-dataset ${dataset} --include-reference --in-order --output-tsv "${refName}_nextclade.tsv" consensus.fasta

    cut -f 2- ${refName}_nextclade.tsv > work/${refName}_nextclade_tmp.tsv
    mv work/${refName}_nextclade_tmp.tsv ${refName}_nextclade.tsv

    rm -rf work
    """
}

process merge_covsummary {

    // Merge coverage files generated for the same reference

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { refName }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.merge_covsummary.todo == 1

    input:
    tuple(val(refName), path("*"))

    output:
    tuple val(refName), path("varcall/${refName}_cov.cov")

    script:
    """
    mkdir -p varcall/
    cat *.cov > varcall/${refName}_cov.cov
    """
}

process merge_fastqstat {

    // Merge fastqstat files generated from fastqc

    label 'varcall'
    storeDir (params.result)
    echo false
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.merge_fastqstat.todo == 1

    input:
    file 'fastqstat/*'

    output:
    file "${runprefix}_qc_fastqstat.tsv"

    script:
    runprefix = (params.prefix =~ /([^\_]+)_(.+)/)[0][1]
    """
    mkdir -p work
    cat fastqstat/*.tsv > ${runprefix}_qc_fastqstat.tsv
    """
}

process merge_bamstat {

    // Merge bamstat files generated for the same reference

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { refName }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.merge_bamstat.todo == 1

    input:
    tuple(val(refName), path("*"))

    output:
    tuple val(refName), path("varcall/${refName}_bamstat.tsv")

    script:
    """
    mkdir -p varcall/
    cat *.tsv > varcall/${refName}_bamstat.tsv
    """
}

process merge_countsummary {

    // Merge varcount files generated for the same reference

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { refName }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.merge_countsummary.todo == 1

    input:
    tuple(val(refName), path("*"))

    output:
    tuple val(refName), path("varcall/${refName}_varcount.tsv")

    script:
    """
    mkdir -p varcall/
    cat *.tsv > varcall/${refName}_varcount.tsv
    """
}

process merge_coinfsummary {

    // Merge coinf files generated for the same reference

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { refName }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.merge_coinfsummary.todo == 1

    input:
    tuple(val(refName), path("*"))

    output:
    tuple val(refName), path("varcall/${refName}_coinf.tsv")

    script:
    """
    mkdir -p varcall/
    cat *.tsv > varcall/${refName}_coinf.tsv
    """
}

process merge_poscsummary {

    // Merge premapping positive controls coverage files generated for each sample

    label 'varcall'
    storeDir (params.result)
    echo false
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.merge_poscsummary.todo == 1

    input:
    file 'posc/*'
    
    output:
    file "${runprefix}_premap_posc.tsv"
    
    script:
    runprefix = (params.prefix =~ /([^\_]+)_(.+)/)[0][1]
    """
    python2.7 ${params.nfpath}/script/rsv/merge_tables.py posc/*.tsv > ${runprefix}_premap_posc.tsv
    """
}

process merge_veriftable {

    //

    label 'varcall'
    storeDir (params.result)
    echo false
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.merge_veriftable.todo == 1

    input:
    file 'verif/*'
    
    output:
    file "${runprefix}_premap_verif.tsv"
    
    script:
    runprefix = (params.prefix =~ /([^\_]+)_(.+)/)[0][1]
    """
    python2.7 ${params.nfpath}/script/rsv/merge_tables.py verif/*.tsv > ${runprefix}_premap_verif.tsv
    """
}

process gentsv_summary {

    // Generate a summary of the run for each reference used

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { refName }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.gentsv_summary.todo == 1

    input:
    tuple val(refName), path("*")

    output:
    tuple val(refName), path("varcall/summary/${refName}_summary.tsv")

    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8
    mkdir -p work varcall/summary/
    
    CONSENSUS=\$(find ./ -name *_consensus.fna)
    COV=\$(find ./ -name *_cov.cov)
    FASTQSTAT=\$(find ./ -name *_fastqstat.tsv)
    BAMSTAT=\$(find ./ -name *_bamstat.tsv)
    VARCOUNT=\$(find ./ -name *_varcount.tsv)
    COINF=\$(find ./ -name *_coinf.tsv)
    POSC=\$(find ./ -name *_posc.tsv)
    CONTATABLE=\$(find ./ -name *_conta.tsv)
    VERIFTABLE=\$(find ./ -name *_verif.tsv)
    BLASTTABLE=\$(find ./ -name *_blast.tsv)

    ARG="--output varcall/summary/${refName}_summary.tsv --refname ${refName} --runname ${params.prefix} --reference ${refName}.fna --cons \${CONSENSUS} --coverage \${COV}"
    if [[ -s \${POSC} ]]; then ARG+=" --posc \${POSC}"; fi
    if [[ -s \${CONTATABLE} ]]; then ARG+=" --conta \${CONTATABLE}"; fi
    if [[ -s \${VERIFTABLE} ]]; then ARG+=" --veriftable \${VERIFTABLE}"; fi
    if [[ -s \${FASTQSTAT} ]]; then ARG+=" --fastqstat-table \${FASTQSTAT}"; fi
    if [[ -s \${BAMSTAT} ]]; then ARG+=" --bamstat-table \${BAMSTAT}"; fi
    if [[ -s \${VARCOUNT} ]]; then ARG+=" --count-table \${VARCOUNT}"; fi
    if [[ -s \${COINF} ]]; then ARG+=" --coinf-table \${COINF}"; fi
    if [[ -s \${BLASTTABLE=} ]]; then ARG+=" --blast \${BLASTTABLE=}"; fi

    echo \$ARG
    /usr/bin/Rscript ${params.nfpath}/script/rsv/make_summary.R \$ARG
    """
}

process merge_summary {

    // Merge all summaries of the run

    label 'varcall'
    storeDir (params.result)
    echo false
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.merge_summary.todo == 1

    input:
    path("work/*")

    output:
    path("${runprefix}_summary.tsv")

    script:
    runprefix = (params.prefix =~ /([^\_]+)_(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8

    awk -f ${params.nfpath}/script/rsv/merge_tsv.awk work/*.tsv > work/${runprefix}_summary.tsv
    python3 ${params.nfpath}/script/rsv/recalc_summary.py -i work/${runprefix}_summary.tsv -o work/${runprefix}_summary_recalc.tsv
    head -n 1 work/${runprefix}_summary_recalc.tsv > "${runprefix}_summary.tsv"
    tail -n+2 "work/${runprefix}_summary_recalc.tsv" | sort -t\$'\\t' -k1 >> "${runprefix}_summary.tsv"

    #rm -rf work
    """
}

process gentsv_validation {

    // Generate a validation report of the run

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { refName }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.gentsv_validation.todo == 1

    input:
    tuple val(refName), path("*"), val(expectedmatrix), val(likelymatrix), val(gff), val(poitable)

    output:
    path("${runprefix}_${type}_fastfinder.csv"), emit: fastfinder
    path("${runprefix}_${type}_indel.tsv"), emit: indel
    path("${runprefix}_${type}_validation.tsv"), emit: validation
    path("${runprefix}_${type}_shortval.tsv"), emit: shortval
    path("${runprefix}_${type}_control.tsv"), emit: control
    path("${runprefix}_${type}_metric.tsv"), emit: metric
    path("${runprefix}_${type}_runsummary.tsv"), emit: runsummary

    script:
    runprefix = (params.prefix =~ /([^\_]+)_(.+)/)[0][1]
    type = (expectedmatrix =~ /([^\_]+)_([^\_]+)_(.+)/)[0][2]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8
    mkdir -p work

    SUMMARY=\$(find ./ -name *_summary.tsv)
    NEXTCLADE=\$(find ./ -name *_nextclade.tsv)
    PANGOLIN=\$(find ./ -name *_pangolin.csv)
    VARTABLE=\$(find ./ -name *_table.tsv)
    ARG="--outdir \$PWD --mode ${params.gentsv_validation["mode"]} --runprefix $runprefix --varcount_threshold ${params.gentsv_validation["varcount_threshold"]} --dp_threshold ${params.gentsv_validation["dp_threshold"]} --coinf_threshold ${params.gentsv_validation["coinf_threshold"]} --cov_minok ${params.gentsv_validation["cov_minok"]} --cov_maxneg ${params.gentsv_validation["cov_maxneg"]}"
    if [[ -s \${SUMMARY} ]]; then ARG+=" --summary \${SUMMARY}"; fi
    if [[ -s \${NEXTCLADE} ]]; then ARG+=" --nextclade \${NEXTCLADE}"; fi
    if [[ -s \${PANGOLIN} ]]; then ARG+=" --pangolin \${PANGOLIN}"; fi
    if [[ -s \${VARTABLE} ]]; then ARG+=" --vartable \${VARTABLE}"; fi
    if [[ "$expectedmatrix" != "" ]]; then ARG+=" --expectedmatrix $expectedmatrix"; fi
    if [[ "$likelymatrix" != "" ]]; then ARG+=" --likelymatrix $likelymatrix"; fi
    if [[ "$gff" != "" ]]; then ARG+=" --gff $gff"; fi
    if [[ "$poitable" != "" ]]; then ARG+=" --poitable $poitable"; fi
    echo \$ARG

    python3 ${params.nfpath}/script/rsv/gen_techval.py --nextcladeversion "${params.gentsv_validation["nextclade_version"]}" \$ARG

    """
}

process merge_fastfinder {

    // Merge fastfinder files to keep validated result

    label 'varcall'
    storeDir (params.result)
    echo false
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.merge_fastfinder.todo == 1

    input:
    path("*")

    output:
    path("${runprefix}_fastfinder.csv")

    script:
    runprefix = (params.prefix =~ /([^\_]+)_(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8
    mkdir -p work

    VAL_FILES=\$(find $PWD -name "*_validation.tsv" | tr '\\n' ',' | sed 's/,\$//g')
    echo "\$VAL_FILES"

    /usr/bin/Rscript ${params.nfpath}/script/rsv/merge_fastfinder.R --val_input "\$VAL_FILES" --output ${runprefix}_fastfinder.csv
    if [[ ${params.gentsv_validation["mode"]} == "rsv" ]]; then sed -i 's/,GABILL,/,VRSMETA,/g' ${runprefix}_fastfinder.csv; sed -i 's/,COMGRAB,/,COMVRS,/g' ${runprefix}_fastfinder.csv ; fi
    

    """
}

process gentree_fasttree {

    // Generate a phylogenetic tree from all sequences generated for the same reference, adding a set of reference sequences
    // TODO Find correct column index, in case columns change

    label 'highcpu'
    storeDir (params.result)
    echo false
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.gentree_fasttree.todo == 1

    input:
    tuple(val(refName), path(summary), path(consensus), path(phylogeny))

    output:
    path("${refName}_validated.nwk")

    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8

    cat "${consensus}" "phylogeny/*" >> toalign.fasta
    mafft --thread $task.cpus --auto --op 5.0 --preservecase --reorder toalign.fasta > aligned.fasta
    fasttree -nt aligned.fasta > ${refName}_validated.nwk
    """
}

process qc_covplot {

    // Generate a coverage plot of a sample mapped on a given reference from a bed file

    label 'varcall'
    storeDir params.result
    echo false
    tag { sampleId }

    when :
    params.qc_covplot.todo == 1 && bga.size() > 0

    input :
    tuple(val(refName), val(sampleId), path(bga))

    output : 
    tuple val(refName), path("varcall/plot/${refName}/${sampleId}_covplot.RDS")

    script :
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p varcall/plot/${refName}

    /usr/bin/Rscript ${params.nfpath}/script/rsv/gen_covplot.R --bga_file ${bga} --sample_name ${sampleId} --ref_name ${refName} --output varcall/plot/${refName}/${sampleId}_covplot.RDS
    """
}

process merge_covplot {

    // Merge coverage plots in a single file by reference

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { refName }

    when : 
    params.merge_covplot.todo == 1 

    input :
    tuple(val(refName), path("work/*"))

    output :
    path("varcall/${refName}_covplot.pdf")

    script :
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p varcall

    /usr/bin/Rscript ${params.nfpath}/script/rsv/merge_covplot.R --output "varcall/${refName}_covplot.pdf" --plot_folder work --chunk ${params.merge_covplot.chunk}
    """
}

process qc_multiqc {
    
    // Generate QC report of the run from different sources

    label 'varcall'
    storeDir params.result
    echo false

    when :
    params.qc_multiqc.todo == 1

    input : 
    path("work/*")    

    output :
    path("${runprefix}_multiqc.html")

    script :
    runprefix = (params.prefix =~ /([^\_]+)_(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p work/

    multiqc work/

    mv multiqc_report.html "${runprefix}_multiqc.html"
    """
}

process plot_tpos {

    // Realign soft-clipped reads, trying to align all reads of a given region in the same way, permit to improve detection of indels
    // TODO update script for new summary heads

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { refName }

    when:
    params.plot_tpos.todo == 1

    input:
    tuple(val(refName), path(summary))
    
    output:
    path("${runprefix}_TposEvol.pdf")

    script:
    runprefix = (params.prefix =~ /([^\_]+)_(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8
    mkdir -p work
    echo ${summary}
    /usr/bin/Rscript ${params.nfpath}/script/rsv/plot_seqplate.R \
        --summary ${summary} \
        --output ${runprefix} 
    """
}
