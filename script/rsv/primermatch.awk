#!/bin/awk -f
#v0.0.1
#USAGE: while read line; do i=$(($i+1)); awk -v loop="${i}" -v prefix="IDTRVP_RSV" -v seq="${line}" -f primermatch.awk /srv/scratch/iai/seqmet/db/reference/rsv/EPIISL1653999/EPIISL1653999.fna >> primer_rsvb_f.bed;done<primer_rsvb_f.txt

BEGIN {
    OFS="\t"
    sep = "|"
    gsub("R", ".", seq)
    gsub("Y", ".", seq)
    gsub("S", ".", seq)
    gsub("W", ".", seq)
    gsub("K", ".", seq)
    gsub("M", ".", seq)
    gsub("B", ".", seq)
    gsub("D", ".", seq)
    gsub("H", ".", seq)
    gsub("V", ".", seq)
    gsub("N", ".", seq)
    regexp = toupper(seq)
    for (i=1; i<=length(seq); i++) {
        regexp = regexp sep substr(seq,1,i-1) "." substr(seq,i+1)
    }
    split(regexp, a, "|")
    for (j=1; j<=length(a); j++) {
        for (i=1; i<=length(seq); i++) {
            regexp = regexp sep substr(a[j],1,i-1) "." substr(a[j],i+1)
        }
    }
    split(regexp, a, "|")
    for (j=1; j<=length(a); j++) {
        for (i=1; i<=length(seq); i++) {
            regexp = regexp sep substr(a[j],1,i-1) "." substr(a[j],i+1)
        }
    }
    split(regexp, a, "|")
    seqrc = ""
    for(i=length(seq); i!=0 ;i--) {
        seqrc = seqrc substr(seq,i,1);
    }
    gsub("A", "t", seqrc)
    gsub("T", "a", seqrc)
    gsub("G", "c", seqrc)
    gsub("C", "g", seqrc)
    regexp = toupper(seqrc)
    for (i=1; i<=length(seqrc); i++) {
        regexp = regexp sep substr(seqrc,1,i-1) "." substr(seqrc,i+1)
    }
    split(regexp, b, "|")
    for (j=1; j<=length(b); j++) {
        for (i=1; i<=length(seqrc); i++) {
            regexp = regexp sep substr(b[j],1,i-1) "." substr(b[j],i+1)
        }
    }
    split(regexp, b, "|")
    for (j=1; j<=length(b); j++) {
        for (i=1; i<=length(seqrc); i++) {
            regexp = regexp sep substr(b[j],1,i-1) "." substr(b[j],i+1)
        }
    }
    split(regexp, b, "|")
}
{
    if (substr($0, 1, 1) == ">") {
        seqname=$0
        next
    }
    for (i=1; i<=length(a); i++) {
        if (match($0,a[i])) {
            print substr(seqname, 2), (match($0,a[i]))-1, (match($0,a[i]))+length(a[i])-1, prefix"_"loop"_F", prefix, "+"
            exit 0
        }
    }
    for (i=1; i<=length(b); i++) {
        if (match($0,b[i])) {
            print substr(seqname, 2), (match($0,b[i]))-1, (match($0,b[i]))+length(b[i])-1, prefix"_"loop"_R", prefix, "-"
            exit 0
        }
    }
    print substr(seqname, 2), "NA", "NA", prefix"_"loop"_NA", a[1], "NA"
}
