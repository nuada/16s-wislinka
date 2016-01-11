#!/bin/bash

TEMP=/tmp/pipeline-$(date +%F-%H%M%S)
OTU_REFERENCE=/usr/local/lib/python2.7/dist-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta
REMOVE_CHIMERAS=yes

function log {
	echo "*** $@"
}

function trim {
	sample=$1
	# Trim primers discarding reads without primers
	opts="--match-read-wildcards -e 0.15 -O 8"
	cutadapt ${opts} -g ^CCTACGGGNGGCWGCAG     --discard-untrimmed -o ${TEMP}/${sample}.01.fastq.gz -p ${TEMP}/${sample}.02.fastq.gz ${DATA}/${sample}.r{1,2}.fastq.gz > ${RESULTS}/${sample}-trimming.log
	cutadapt ${opts} -g ^GACTACHVGGGTATCTAATCC --discard-untrimmed -o ${TEMP}/${sample}.12.fastq.gz -p ${TEMP}/${sample}.11.fastq.gz ${TEMP}/${sample}.0{2,1}.fastq.gz >> ${RESULTS}/${sample}-trimming.log
	# Trim primer rc and adapters
	opts="--match-read-wildcards -e 0.15 -O 10 -m 200"
	cutadapt ${opts} -a GGATTAGATACCCBDGTAGTCCTGTCTCTTATACACATCT -o ${TEMP}/${sample}.21.fastq.gz -p ${TEMP}/${sample}.22.fastq.gz ${TEMP}/${sample}.1{1,2}.fastq.gz >> ${RESULTS}/${sample}-trimming.log
	cutadapt ${opts} -a CTGCWGCCNCCCGTAGGCTGTCTCTTATACACATCT -o ${TEMP}/${sample}.r2.fastq.gz -p ${TEMP}/${sample}.r1.fastq.gz ${TEMP}/${sample}.2{2,1}.fastq.gz >> ${RESULTS}/${sample}-trimming.log
}

function preprocess_samples {
	sample=$1
	cpus=$((${CPUS}/${JOBS}))

	log "Trim PCR primers ${sample}"
	trim ${sample}

	log "Join reads ${sample}"
	join_paired_ends.py \
		-f ${TEMP}/${sample}.r1.fastq.gz \
		-r ${TEMP}/${sample}.r2.fastq.gz \
		-o ${TEMP}/${sample}-joined

	# Per sample mapping file
	( head -n 1 ${MAPPING}; grep "^${sample}" ${MAPPING} ) > ${TEMP}/${sample}-map

	log "Convert FASTQs to seqs.fna ${sample}"
	if [[ -s ${TEMP}/${sample}-joined/fastqjoin.join.fastq ]]; then
		split_libraries_fastq.py \
			-i ${TEMP}/${sample}-joined/fastqjoin.join.fastq,${TEMP}/${sample}-joined/fastqjoin.un1.fastq \
			--sample_id ${sample},${sample} \
			-o ${TEMP}/split-${sample} \
			-m ${TEMP}/${sample}-map \
			-q 20 \
			--barcode_type 'not-barcoded'
	else
		# No reads were joined!
		split_libraries_fastq.py \
			-i ${TEMP}/${sample}-joined/fastqjoin.un1.fastq \
			--sample_id ${sample} \
			-o ${TEMP}/split-${sample} \
			-m ${TEMP}/${sample}-map \
			-q 20 \
			--barcode_type 'not-barcoded'
	fi

	log "Done ${sample}"
}

function prepare_seqs {
	./count_reads_per_sample.py $(for sample in ${SAMPLES}; do echo -n "${DATA}/${sample}.r1.fastq.gz "; done) > ${RESULTS}/reads_per_sample.1.raw &

	for sample in ${SAMPLES}; do
		preprocess_samples ${sample} &
		while [[ $(jobs -p | wc --lines) -ge ${JOBS} ]]; do sleep 1; done
	done
	wait

	./count_reads_per_sample.py $(for sample in ${SAMPLES}; do echo -n "${TEMP}/${sample}.r1.fastq.gz "; done) > ${RESULTS}/reads_per_sample.2.trimmed &
	for sample in ${SAMPLES}; do
		ln -s ${TEMP}/${sample}-joined/fastqjoin.join.fastq ${TEMP}/${sample}.joined.fastq
	done
	./count_reads_per_sample.py $(for sample in ${SAMPLES}; do echo -n "${TEMP}/${sample}.joined.fastq "; done) > ${RESULTS}/reads_per_sample.3.joined &

	log 'Concatenate all samples'
	cat ${TEMP}/split-*/seqs.fna > ${RESULTS}/seqs.raw.fna
	./count_seqs_per_sample.py ${RESULTS}/seqs.raw.fna > ${RESULTS}/reads_per_sample.4.seqs &

	if [[ ${REMOVE_CHIMERAS} == "yes" ]]; then
		log 'Remove chimeric sequences'
		identify_chimeric_seqs.py -i ${RESULTS}/seqs.raw.fna -m usearch61 -o ${RESULTS}/usearch_checked_chimeras -r ${OTU_REFERENCE}
		filter_fasta.py -f ${RESULTS}/seqs.raw.fna -o ${RESULTS}/seqs.chimeras_filtered.fna -s ${RESULTS}/usearch_checked_chimeras/chimeras.txt -n
		ln -s ${RESULTS}/seqs.chimeras_filtered.fna ${RESULTS}/seqs.fna
		./count_seqs_per_sample.py ${RESULTS}/seqs.fna > ${RESULTS}/reads_per_sample.5.no_chimeras &
	else
		ln -s ${RESULTS}/seqs.raw.fna ${SEQS}
	fi

	wait
}

function build_otu_table {
	log 'Pick OTUs'
	pick_open_reference_otus.py -i ${SEQS} -o ${RESULTS}/otus -m usearch61 -a -O ${CPUS}
	( cd ${RESULTS}/otus && ln -s otu_table_mc2_w_tax_no_pynast_failures.biom otu_table.biom )

	log 'Summarize OTU table'
	biom summarize-table -i ${OTU_TABLE} -o ${RESULTS}/otus/otu_table_summary.txt
}

function analyze_otus {
	log 'Not implemented!!!'
}

function filter {
	out_dir=$1
	filter="$2"

	mkdir ${out_dir}
	filter_samples_from_otu_table.py \
		-i ${OTU_TABLE} \
		-o ${out_dir}/otu_table.biom \
		-m ${MAPPING} \
		--output_mapping_fp ${out_dir}/mapping.txt \
		-s ${filter}
}

function comparison {
	out_dir=$1
	category=$2

	core_diversity_analyses.py \
		-i ${out_dir}/otu_table.biom \
		-o ${out_dir}/core_output \
		-m ${out_dir}/mapping.txt \
		-t ${TREE} \
		-c ${category} \
		-e ${rarefaction} \
		-a -O ${CPUS}

	zcat ${out_dir}/core_output/table_even${rarefaction}.biom.gz > ${out_dir}/core_output/table_even${rarefaction}.biom
	for ctg in $(echo ${category} | tr ',' ' '); do
		for level in 5 6; do
			summarize_taxa.py \
				-i ${out_dir}/core_output/table_even${rarefaction}.biom \
				-L ${level} \
				-o ${out_dir}/taxa_summary_L${level}
			group_significance.py \
				-i ${out_dir}/taxa_summary_L${level}/table_even${rarefaction}_L${level}.biom \
				-o ${out_dir}/otu_significance_L${level}_kruskal_wallis_${category} \
				-m ${out_dir}/mapping.txt \
				-s kruskal_wallis \
				-c ${ctg} &
		done
	done
	wait
}

function comparison_deseq2 {
	out_dir=$1
	category=$2
	c1=$3
	c2=$4
	suffix=$5

	for level in 5 6 7; do
		if [[ ! -d ${out_dir}/diff_otus_L${level} ]]; then
			summarize_taxa.py \
				-i ${out_dir}/otu_table.biom \
				-L ${level} \
				-a \
				-o ${out_dir}/diff_otus_L${level}
		fi
		differential_abundance.py \
			-i ${out_dir}/diff_otus_L${level}/otu_table_L${level}.biom \
			-o ${out_dir}/diff_otus_L${level}/diff_otus${suffix}_DESeq2_nbinom.txt \
			-m ${out_dir}/mapping.txt \
			-a DESeq2_nbinom \
			-c ${category} \
			-x ${c1} \
			-y ${c2} \
			-d &
	done
	wait
}

function correlation {
	out_dir=$1
	shift 1

	filter_samples_from_otu_table.py \
		-i ${out_dir}/otu_table.biom \
		-o ${out_dir}/table_mc${rarefaction}.biom \
		-n ${rarefaction}

	single_rarefaction.py \
		-i ${out_dir}/table_mc${rarefaction}.biom \
		-o ${out_dir}/table_even${rarefaction}.biom \
		-d ${rarefaction}

	for category in $*; do
		log "Correlation with ${category}"
		observation_metadata_correlation.py \
			-i ${out_dir}/table_even${rarefaction}.biom \
			-m ${out_dir}/mapping.txt \
			-c ${category} \
			-s spearman \
			-o ${out_dir}/correlation_spearman_${category}.txt &
	done
	wait
}

function main {
	# Reastart failed analysis, reuse results directory
	if [[ $1 == '--resume' ]]; then
		RESULTS=${2}
		shift 2
	fi
	
	[[ -d ${TEMP} ]] || mkdir -p ${TEMP}
	[[ -d ${RESULTS} ]] || mkdir -p ${RESULTS}
	
	cp ${MAPPING} ${RESULTS}/mapping.txt

	# If given existing results directory, only analyze OTUs
	if [[ -d $1 ]]; then
		MAPPING=${RESULTS}/mapping.txt
		OTU_TABLE=${1}/otus/otu_table.biom
		TREE=${1}/otus/rep_set.tre
	else
		MAPPING=${RESULTS}/mapping.txt
		SEQS=${RESULTS}/seqs.fna
		OTU_TABLE=${RESULTS}/otus/otu_table.biom
		TREE=${RESULTS}/otus/rep_set.tre
		[[ -f ${SEQS} ]] || prepare_seqs
		[[ -d $(dirname ${OTU_TABLE}) ]] || build_otu_table
	fi
	analyze_otus

	cd $(dirname ${RESULTS})
	zip -r ${RESULTS}.zip $(basename ${RESULTS})
	cp -v ${RESULTS}.zip ${DATA}/

	rm -rf ${TEMP}
}
