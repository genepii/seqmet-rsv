#!/bin/awk -f
#v0.0.1

BEGIN {
   OFS="\t"
   FPAT="[^\t]*|\042[^\042]+\042"
   # keep track of the original argument count
   argc_start=ARGC
}

# Read header and process
# header names are stored as array index in the array "header"
# header order is stored in the array header_order
#    header_order[field_index] = header_name
(FNR == 1) && (ARGIND < argc_start) {
    for(i=1;i<=NF;++i) if (!($i in header)) { header[$i]; header_order[++nf_out]=$i } 
    # add file to end of argument list to be reprocessed
    ARGV[ARGC++] = FILENAME
    # process the next file
    nextfile
}

# Print headers in output file
(FNR == 1) && (ARGIND == argc_start) {
    for(i=1;i<=nf_out;++i) printf header_order[i] (i==nf_out ? ORS : OFS)
}

# Use array h to keep track of the column_name and corresponding field_index
# h[column_name] = field_index
(FNR == 1) { delete h; for(i=1;i<=NF;++i) h[$i]=i; next }

# print record
{
    # process all fields
    for(i=1;i<=nf_out;++i) {
        # get field index using h
        j = h[header_order[i]]+0
        # if field index is zero, print empty field
        printf (j == 0 ? "" : $j) (i==nf_out ? ORS : OFS)
    }
}