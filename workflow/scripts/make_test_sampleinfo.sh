#!/usr/bin/env bash

# shellcheck disable=SC2154
RUNINFO="${snakemake_input[0]}"
SAMPLES="${snakemake_input[1]}"
OUTFILE="${snakemake_output[0]}"
LOG="${snakemake_log[0]}"

# shellcheck disable=SC2016
csvtk cut -f "Sample,Run,ScientificName,SampleName" "${RUNINFO}" \
    | csvtk mutate - -n AuthorSample -f SampleName -p "^([A-Z]+-.+[A-Z])[0-9]+-2016" \
    | csvtk mutate - -n SampleAlias -f AuthorSample \
    | csvtk replace -f SampleAlias,AuthorSample -p "(^[^A].+)_[0-9]+$" -r '$1' \
    | csvtk replace -f SampleAlias,AuthorSample -p "(^A[A-Z]+-T[0-9]+)_[0-9]+$" -r '$1' \
    | csvtk replace -f SampleAlias,AuthorSample -p "CLV" -r "CLV_" \
    | csvtk replace -f SampleAlias,AuthorSample -p "-GH1" -r "GH" \
    | csvtk replace -f AuthorSample -p "^[A-Z]+-" -r "" \
    | csvtk replace -f SampleAlias -p "PUN-(BCRD|INJ|LO|PCT|POTR)" -r 'PUN-Y-$1' \
    | csvtk replace -f SampleAlias -p "PUN-(ELF|JMC|LH|MT|UCSD)" -r 'PUN-R-$1' >sraruninfo.csv 2>"${LOG}"

csvtk join sraruninfo.csv "${SAMPLES}" -f "AuthorSample;Sample" >"${OUTFILE}" 2>>"${LOG}"
rm -f sraruninfo.csv
