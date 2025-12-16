#!/usr/bin/env bash

# shellcheck disable=SC2154
INFILE="${snakemake_input[0]}"
OUTFILE="${snakemake_output[0]}"
LOG="${snakemake_log[0]}"

docx2txt "${INFILE}" | tail -n +21 | tr "\n" "," \
    | sed -e 's/,,/","/g' \
    | awk -F"\",\"" '{for (i=1; i<=NF; i++) \
        {printf("%s,", $i); if (i % 6 == 0) print("")}}' \
    | sed -e "s/,$//g" | sed -e "s/^ //g" \
    | sed -e 's/ ,/",/g' | sed -e 's/ssp. puniceus/"ssp. puniceus/g' \
    | qsv rename -n "Sample,Taxon,Latitude,Longitude,% Reads aligned,Seq. Depth" >"${OUTFILE}" 2>"${LOG}"
