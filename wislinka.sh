#!/bin/bash

set -e

source functions.sh

DATA=/data
MAPPING=${DATA}/20151203-mapping.txt
SAMPLES=$(tail -n +2 ${MAPPING} | cut -f 1)
RESULTS=$(pwd)/wislinka-$(date +%F-%H%M%S)

JOBS=8
CPUS=24

function analyze_otus {
	rarefaction=26099
	log 'Core diversity'
	core_diversity_analyses.py -i ${OTU_TABLE} -o ${RESULTS}/core_output -m ${MAPPING} -c EXTRACTION_DATE,KIT,BUFFER -t ${TREE} -e ${rarefaction} -a -O ${CPUS}

	zcat ${RESULTS}/core_output/table_even${rarefaction}.biom.gz > ${RESULTS}/core_output/table_even${rarefaction}.biom
	for level in 2 5; do
		for category in KIT; do
			summarize_taxa.py \
				-i ${RESULTS}/core_output/table_even${rarefaction}.biom \
				-L ${level} \
				-o ${RESULTS}/taxa_summary_L${level}
			group_significance.py \
				-i ${RESULTS}/taxa_summary_L${level}/table_even${rarefaction}_L${level}.biom \
				-o ${RESULTS}/otu_significance_L${level}_kruskal_wallis_${category} \
				-m ${RESULTS}/mapping.txt \
				-s kruskal_wallis \
				-c ${category}
		done
	done
}

main "$@"
