#!/bin/bash

GATK="/data/Resources/Software/Javas/GenomeAnalysisTK.jar"
REFERENCE="hg38.bwa"
CORES=8
BED="/data/Resources/BEDs/nexterarapidcapture_exome_targetedregions_v1.2_hg38.bed"
ANNO="/data/Resources/Software/annovar/"

## args
VCF=NULL
PED=NULL
TYPE="GENE"
PATHWAY=NULL
MODE=NULL
TESTS="all"
OUTPUT=NULL
AF_ref=0.01
AF_all=0.05
AF_case=0.2
AF_cont=0.2
MINGQ=30
MINQ=30
MINDP=10
MISS=0.95
CONSEQS="nonframeshift_deletion,nonframeshift_insertion,frameshift_deletion,stopgain,frameshift_insertion,splicing,nonsynonymous_SNV"
PROJECT="default_project"

if [[ $# -eq 0 ]]; then
        echo -e "\n[ASSOCIATION TESTS][MAIN] Error: You need to provide at least SOME arguments! Try using -h / --help for documentation and examples!\n"
        exit
fi

for arg in "$@"; do
        if [[ "$arg" == "-h" ]] || [[ "$arg" == "--help" ]]; then
		echo -e "
	-i  --input
	-c  --cases
	-t  --type
	-p  --pathway
	-nt --threads
	-m  --mode
	-a  --analysis_tests
	-gf --global_MAF
	-if --internal_AF
	-caf --case_AF
	-ctf --control_AF
	-GQ --min_GQ
	-Q --min_QUAL
	-DP --min_DP
	-miss --max_missing
	-csq --conseqs   nonframeshift_deletion,nonframeshift_insertion,frameshift_deletion,stopgain,frameshift_insertion,splicing,nonsynonymous_SNV,synonymous_SNV,stoploss

            "
                echo -e "\n"
                exit
        fi
done

##arugement parsing block
while [[ $# > 1 ]]
	do 
	key="$1"
	case $key in
    		-i|--input)
    		VCF=$2
    		shift
    		;;
    		-c|--cases)
    		PED=$2
    		shift
    		;;
    		-t|--type)
    		TYPE=$2
    		shift
    		;;
    		-p|--pathway)
    		PATHWAY=$2
    		shift
    		;;
    		-nt|--threads)
    		CORES=$2
    		shift
    		;;
    		-m|--mode)
    		MODE=$2
    		shift
    		;;
    		-a|--analysis_tests)
    		TESTS=$2
    		shift
    		;;
		-o|--output)
                OUTPUT=$2
                shift
                ;;
		-pro|--project)
                PROJECT=$2
                shift
                ;;
		-gf|--global_MAF)
                AF_ref=$2
                shift
                ;;
		-if|--internal_AF)
                AF_all=$2
                shift
                ;;
		-caf|--case_AF)
                AF_case=$2
                shift
                ;;
		-ctf|--control_AF)
                AF_cont=$2
                shift
                ;;
		-GQ|--min_GQ)
                MINGQ=$2
                shift
                ;;
		-Q|--min_QUAL)
                MINQ=$2
                shift
                ;;
		-DP|--min_DP)
                MINDP=$2
                shift
                ;;
		-miss|--max_missing)
                MISS=$2
                shift
                ;;
		-csq|--conseqs)
                CONSEQS=$2
                shift
                ;;
	esac
	shift
done


if [[ "$VCF" == "NULL" ]]; then
	echo -e "\n[ASSOCIATION TESTS][MAIN] Error: No input vcf given, provided by -i / --input\nExiting Now"
	exit
fi

if [[ "$PED" == "NULL" ]]; then
        echo -e "\n[ASSOCIATION TESTS][MAIN] Error: No case/control file given, provided by -c / --cases\nExiting Now"
        exit
fi

if [[ "$TYPE" == "PATHWAY" ]] && [[ "$PATHWAY" == "NULL" ]]; then
	echo -e "\n[ASSOCIATION TESTS][MAIN] Error: Pathway analysis selected but no pathway file provided by -p / --pathway\nExiting Now"
        exit
fi

if [[ "$OUTPUT" == "NULL" ]]; then
        echo -e "\n[ASSOCIATION TESTS][MAIN] Error: No output folder given, provided by -o / --output\nExiting Now"
        exit
fi

echo -e "[ASSOCIATION TESTS][MAIN] Input VCF provided - ${VCF}"
echo -e "[ASSOCIATION TESTS][MAIN] Case/Control file provided - ${PED}"

if [[ ! -d ${OUTPUT}${PROJECT} ]]; then
	mkdir ${OUTPUT}${PROJECT}
fi

echo -e "[ASSOCIATION TESTS][MAIN] Analysis start - $date" > assoc_log.txt
echo -e "[ASSOCIATION TESTS][MAIN] Input VCF -${VCF}" >> assoc_log.txt
echo -e "[ASSOCIATION TESTS][MAIN] Case/Control file - ${PED}" >> assoc_log.txt
echo -e "[ASSOCIATION TESTS][MAIN] Analysis type - ${TYPE}" >> assoc_log.txt
echo -e "[ASSOCIATION TESTS][MAIN] Pathway file - ${PATHWAY}" >> assoc_log.txt
echo -e "[ASSOCIATION TESTS][MAIN] Script mode - ${MODE}" >> assoc_log.txt
echo -e "[ASSOCIATION TESTS][MAIN] Tests - ${TESTS}" >> assoc_log.txt
echo -e "[ASSOCIATION TESTS][MAIN] Output folder - ${OUTPUT}" >> assoc_log.txt
echo -e "[ASSOCIATION TESTS][MAIN] Global max_MAF - ${AF_ref}" >> assoc_log.txt
echo -e "[ASSOCIATION TESTS][MAIN] internal max_MAF - ${AF_all}" >> assoc_log.txt
echo -e "[ASSOCIATION TESTS][MAIN] case max_MAF - ${AF_case}" >> assoc_log.txt
echo -e "[ASSOCIATION TESTS][MAIN] control max_MAF ${AF_cont}" >> assoc_log.txt
echo -e "[ASSOCIATION TESTS][MAIN] min_GQ - ${MINGQ}" >> assoc_log.txt
echo -e "[ASSOCIATION TESTS][MAIN] min_qual - ${MINQ}" >> assoc_log.txt
echo -e "[ASSOCIATION TESTS][MAIN] min_read_depth - ${MINDP}" >> assoc_log.txt
echo -e "[ASSOCIATION TESTS][MAIN] max_missing - ${MISS}" >> assoc_log.txt
echo -e "[ASSOCIATION TESTS][MAIN] Mutational consequences - ${CONSEQS}" >> assoc_log.txt
echo -e "[ASSOCIATION TESTS][MAIN] Project name - ${PROJECT}" >> assoc_log.txt

if [[ "$MODE" != "no_prep" ]]; then
	echo -e "[ASSOCIATION TESTS][MAIN] Copying input VCF ${VCF} to analysis folder" | tee -a assoc_log.txt
	cp tests/* ${OUTPUT}${PROJECT}/
	cp ${VCF} ${OUTPUT}${PROJECT}/input_file.vcf
	mv assoc_log.txt ${OUTPUT}${PROJECT}/
	cd ${OUTPUT}${PROJECT}/
	VCF="input_file"

	## Further vcf filtering removing sites with >10% missing particularly
	echo -e "[ASSOCIATION TESTS][MAIN] Filtering VCF - minGQ 30 minDP 10 max-missing 0.8 minQ 30" | tee -a assoc_log.txt
	vcftools --vcf ${VCF}.vcf --hwe 0.05 --non-ref-ac-any 1 --minGQ ${MINGQ} --minDP ${MINDP} --max-missing ${MISS} --minQ ${MINQ} --recode --out ${VCF}
	mv ${VCF}.recode.vcf ${VCF}.filt.vcf
	
	echo -e "[ASSOCIATION TESTS][MAIN] Left aligning and spliting multiallelic sites" | tee -a assoc_log.txt
	## Spliting multiallelics 
	java -Xmx30g -jar ${GATK} \
		-T LeftAlignAndTrimVariants \
		-R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa \
		--variant ${VCF}.filt.vcf \
		-o ${VCF}.filt.bi.vcf \
		--splitMultiallelics

	## Generate list of chr files for
	echo -e "[ASSOCIATION TESTS][MAIN] Spliting VCF by chromosome"  | tee -a assoc_log.txt
	sed -n '/#CHROM/,$p' ${VCF}.filt.bi.vcf | cut -f1 | sort -u | grep -v '#CHROM' > chr_list
	cat chr_list | xargs -n1 -P${CORES} -I {} mkdir temp_{}
	cat chr_list | xargs -n1 -P${CORES} -I {} vcftools --vcf ${VCF}.filt.bi.vcf --chr {} --recode --recode-INFO-all --out temp_{}/${VCF}.splt.{}
	
	echo -e "[ASSOCIATION TESTS][MAIN] Performing ANNOVAR annotation"  | tee -a assoc_log.txt
	cat chr_list | xargs -n1 -P${CORES} -I {} ${ANNO}table_annovar.pl temp_{}/${VCF}.splt.{}.recode.vcf ${ANNO}humandb/ -vcfinput -buildver hg38 -out temp_{}/${VCF}_{} -remove -protocol refGene,exac03,dbnsfp30a -operation g,f,f -nastring .
	
	echo -e "[ASSOCIATION TESTS][MAIN] Concatenating annotated VCF files" | tee -a assoc_log.txt
	cat chr_list | xargs -n1 -P${CORES} -I {} mv temp_{}/${VCF}_{}.hg38_multianno.vcf ${VCF}_{}.hg38_multianno.vcf
	rm -r temp_*
	sed '/#CHROM/,$d' ${VCF}.filt.bi.vcf > header

	for i in `cat chr_list`; do
		cat header ${VCF}_${i}.hg38_multianno.vcf > ${VCF}_${i}.hg38_multianno.header.vcf
	done

	ls *.header.vcf > vcf_list
	bcftools concat -o ${VCF}.FINAL.vcf -Ov -f vcf_list
	rm *_multianno*.vcf

	## Remove header
	sed -i -n '/#CHROM/,$p' ${VCF}.FINAL.vcf
	echo -e "[ASSOCIATION TESTS][MAIN] Cleaning up intermediate files" | tee -a assoc_log.txt
	rm chr_list
	rm header
	rm vcf_list
	rm ${VCF}.filt*
	rm *.log
	rm *.idx
else
	## Run vcf processing script
	echo -e "[ASSOCIATION TESTS][MAIN] Skipping VCF processing (-m ${MODE})" | tee -a assoc_log.txt
	mv assoc_log.txt ${OUTPUT}${PROJECT}/
	cd ${OUTPUT}${PROJECT}/
	VCF="input_file"
fi

echo -e "[ASSOCIATION TESTS][MAIN] Generating association RData file for analysis" | tee -a assoc_log.txt
Rscript vcf_prep.R ${VCF}.FINAL.vcf ${PED} ${CORES} ${AF_ref} ${AF_all} ${AF_case} ${AF_cont} ${CONSEQS}

if [[ "$TESTS" == "all" ]]; then
	echo -e "[ASSOCIATION TESTS][MAIN] Performing ALL statistical tests (-a ${TESTS})" | tee -a assoc_log.txt
fi

if echo ${TESTS} | grep -q "fishers" || [[ "$TESTS" == "all" ]]; then
	echo -e "[ASSOCIATION TESTS][MAIN] Starting Fisher's exact testing" | tee -a assoc_log.txt
	Rscript fishers_tests.R
	echo -e "[ASSOCIATION TESTS][MAIN] Fisher's exact testing finished" | tee -a assoc_log.txt
fi

if echo ${TESTS} | grep -q "skat" || [[ "$TESTS" == "all" ]]; then
 	echo -e "[ASSOCIATION TESTS][MAIN] Starting SKAT burden testing" | tee -a assoc_log.txt
	Rscript SKAT_test.R
	echo -e "[ASSOCIATION TESTS][MAIN] SKAT burden testing finished" | tee -a assoc_log.txt
fi

if echo ${TESTS} | grep -q "bevimed" || [[ "$TESTS" == "all" ]]; then
	echo -e "[ASSOCIATION TESTS][MAIN] Starting Bevimed baysian testing" | tee -a assoc_log.txt
	Rscript bevimed.R ${TYPE} ${PATHWAY}
	echo -e "[ASSOCIATION TESTS][MAIN] Bevimed baysian testing finished" | tee -a assoc_log.txt
fi
