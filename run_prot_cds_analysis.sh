#!/bin/bash
# -*- coding: utf-8 -*-
#
#
#  run_prot_cds_analysis.sh
#  
#  Copyright 2019 Nikos Vakirlis <nvakirlis@nvakirlis-X399-DESIGNARE-EX>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

#This script generates the data that is used as an input for the classifier
#It uses a number of small scripts found in the PARSE_SCRIPT_DIR directory
#You need to provide the names of the directories for the specific tools used
#and also the header files with the column names.
#The tools used are : iupred, tango, codonw, segmasker, psipred, CAIJava, tmhmm
#the tango and codonw executables should be in your $PATH
#sometimes codonw crashes when run from the script, see line 337 if that happens
#Biopython, emboss and possibly other libraries need to be installed in the system
#Make sure to have removed descriptions from your fasta headers (leave only the gene identifier) before running this script, 
#as they will surely result to errors from some tools
#The script takes 4 arguments : the input fasta file, the output directory
#where all the files are written, the min sequence length of the sequences to consider (put 0 to consider all the sequences)
#and whether to delete duplicate sequences from the final file (yes/no)
#the file final_dataframe.csv contains all the info but the names of the columns are confusing
#the file final_dataframe_SMALL.csv contains only the features relevant to the classifier. It is created using the
#short script reformat_final_df.r which should work but as I have not tested it, it might crash (just debug if it does, it's a simple script).

#example : run_prot_cds_analysis.sh input_file.fsa ./ 30 yes

#change the paths here as needed
MINPARAMS=4
CAIJAVA_DIR=/Users/nikos/Documents/tools/CAIJava/
PSIPRED_DIR=/Users/nikos/Documents/tools/psipred/
TMHMM_PATH=/Users/nikos/Documents/tools/tmhmm-2.0c/bin/
PARSE_SCRIPT_DIR=./Parsers_etc/
IUPRED_DIR=/Users/nikos/Documents/tools/iupred/
CAI_HEADER=$PARSE_SCRIPT_DIR"header.txt"
FIN_HEADER=$PARSE_SCRIPT_DIR"fin_header.txt"

if [ $# -lt "$MINPARAMS" ]
then
  echo "This script needs at least $MINPARAMS command-line arguments"
  exit 0
fi

CDS=$(readlink -f $1)
OUT_DIR=$(readlink -f $2)
CWD=$(pwd)
MIN_SEQ_LEN=$3
DEL_DUPS=$4

echo "Checking for missing sequences..."

REM_SEQ=$(python ${PARSE_SCRIPT_DIR}remove_small_orfs.py $CDS $MIN_SEQ_LEN ${1%\.*}"_no_missing.fasta")

if [[ ${REM_SEQ:0:2} != "no" ]]
then
		echo "continuing with "${1%\.*}"_no_missing.fasta"
		CDS=${CDS%\.*}"_no_missing.fasta" 
else
		rm ${1%\.*}"_no_missing.fasta"
fi	

JUST_NAME=${1%\.*}
JUST_NAME=${JUST_NAME##*/}

if [ ! -d ${OUT_DIR} ]
then
	mkdir ${OUT_DIR} ;
	echo "Directory "${OUT_DIR}" created" ;
else 
	echo "Directory "${OUT_DIR}" exists" ;
fi

NAMES_FILE=${OUT_DIR}/${JUST_NAME}"_gene_names.txt"
grep ">" ${CDS} | cut -f 1 | tr -d '>' > ${NAMES_FILE}
	
NO_GENES=$(grep ">" ${CDS} | wc -l)

echo "Found : "${NO_GENES}" sequences in file : "${CDS}

PROT=${CDS%\.*}.prot

echo "Translating sequences"

transeq -sequence ${CDS} -outseq ${PROT} &> /dev/null

if [[ $? == 1 ]]
then
	echo "Problem with translation"
	exit 1
fi

echo "Calculating CAI with CAIJava..."

CAI_OUTPUT_FILE=${OUT_DIR}/${1%\.*}"_CAI_output.txt"

cd ${CAIJAVA_DIR}/ :
java -cp biojava-1.4.jar:CAI.jar CAI.CAIJava ${CDS} -s -f ${CAI_OUTPUT_FILE} -i 15 -k 3 -g &> /dev/null

if [[ $? == 1 ]]
then
	echo "Problem with CAIJava"
	exit 1
else
	echo "Done, reformatting header.."
fi	

cd ${CWD}/ :

awk 'NR >=110' ${CAI_OUTPUT_FILE} > .temp
rm ${CAI_OUTPUT_FILE}
cat ${CAI_HEADER} .temp >> ${CAI_OUTPUT_FILE}
rm .temp

echo "Calculating transmembrane properties with TMHMM..."

TMHMM_DIR=${OUT_DIR}/TMHMM/

if [ ! -d ${TMHMM_DIR} ]
then
	mkdir ${TMHMM_DIR} ;
	echo "Directory "${TMHMM_DIR}" created" ;
else 
	echo "Directory "${TMHMM_DIR}" exists" ;
fi

cd ${TMHMM_DIR} :

TMHMM_OUTPUT=${1%\.*}"_tmhmm_output.txt"

${TMHMM_PATH}tmhmm ${PROT} > ${TMHMM_OUTPUT}

if [[ $? == 1 ]]
then
	echo "Problem with tmhmm"
	exit 1
else
	echo "Done, parsing file"
fi	

grep "Exp number of AAs in TMHs" ${TMHMM_OUTPUT} | cut -f 2,9 -d ' ' > ${OUT_DIR}/${1%\.*}"_tmhmm_output_PARSED.txt"

ONE_SEQ_PER_FILE_DIR=${OUT_DIR}/"split_fasta"

if [ ! -d ${ONE_SEQ_PER_FILE_DIR} ]
then
	mkdir ${ONE_SEQ_PER_FILE_DIR} ;
	echo "Directory "${ONE_SEQ_PER_FILE_DIR}" created" ;
else 
	echo "Directory "${ONE_SEQ_PER_FILE_DIR}" exists" ;
fi

cd ${ONE_SEQ_PER_FILE_DIR} :

csplit -zsf "Seq" ${PROT} '/>/' {*}

ls | xargs -i mv {} {}.fasta

echo "Calculating secondary structure with psipred"

PSIPRED_OUT=${OUT_DIR}/PSIPRED_OUTPUT/

if [ ! -d ${PSIPRED_OUT} ]
then
	mkdir ${PSIPRED_OUT} ;
	echo "Directory "${PSIPRED_OUT}" created" ;
else 
	echo "Directory "${PSIPRED_OUT}" exists" ;
fi

cd ${PSIPRED_DIR} :

parallel --GNU './bin/seq2mtx {} > {.}.mtx' ::: ${ONE_SEQ_PER_FILE_DIR}/*fasta &> /dev/null

if [[ $? == 1 ]]
then
	echo "Problem with psipred : seq2mtx"
	exit 1
fi


parallel --GNU 'psipred {} ./data/weights.dat ./data/weights.dat2 ./data/weights.dat3 > {.}.ss' ::: ${ONE_SEQ_PER_FILE_DIR}/*.mtx  &> /dev/null

if [[ $? == 1 ]]
then
	echo "Problem with psipred : psipred"
	exit 1
fi

parallel --GNU './bin/psipass2 ./data/weights_p2.dat 1 1.0 1.0 {.}.ss2 {} > {.}.horiz' ::: ${ONE_SEQ_PER_FILE_DIR}/*ss  &> /dev/null

if [[ $? == 1 ]]
then
	echo "Problem with psipred : psipass2"
	exit 1
fi

cd ${ONE_SEQ_PER_FILE_DIR} :

for i in *fasta ; do NAME=$(grep ">" ${i} | tr -d '>') ; echo ${NAME} > temp ; cat temp ${i%\.*}".ss2" >> temp2 ; mv temp2 ${i%\.*}".ss2" ; rm temp ; done

cd ${OUT_DIR} :
echo "Tidying up psipred files..."
mkdir ${PSIPRED_OUT}"ss"
mkdir ${PSIPRED_OUT}"ss2"
mkdir ${PSIPRED_OUT}"horiz"
mkdir ${PSIPRED_OUT}"mtx"

mv ${ONE_SEQ_PER_FILE_DIR}/*mtx ${PSIPRED_OUT}"mtx"
mv ${ONE_SEQ_PER_FILE_DIR}/*horiz ${PSIPRED_OUT}"horiz"
mv ${ONE_SEQ_PER_FILE_DIR}/*ss ${PSIPRED_OUT}"ss"
mv ${ONE_SEQ_PER_FILE_DIR}/*ss2 ${PSIPRED_OUT}"ss2"

SS2_DIR=${PSIPRED_OUT}"ss2/"

export PARSE_SCRIPT_DIR

parallel --GNU --env PARSE_SCRIPT_DIR 'python ${PARSE_SCRIPT_DIR}parse_ss2.py {}' ::: PSIPRED_OUTPUT/ss2/*ss2 > Secondary_structure_results.txt

echo "Predicting disordered regions with iupred..."

IUPRED_OUT=${OUT_DIR}/IUPRED_OUTPUT/

if [ ! -d ${IUPRED_OUT} ]
then
	mkdir ${IUPRED_OUT} ;
	echo "Directory "${IUPRED_OUT}" created" ;
else 
	echo "Directory "${IUPRED_OUT}" exists" ;
fi

cd ${IUPRED_DIR} :

parallel --GNU './iupred {} glob > {}.disord 2> /dev/null' ::: ${ONE_SEQ_PER_FILE_DIR}/*

for i in ${ONE_SEQ_PER_FILE_DIR}/*disord ; do csplit -szf ${i}_ ${i} '/>/' {*} ; done

rm ${ONE_SEQ_PER_FILE_DIR}/*00

cat ${ONE_SEQ_PER_FILE_DIR}/*01 >> ${IUPRED_OUT}/All_iupred_fasta.fas

rm ${ONE_SEQ_PER_FILE_DIR}/*01

mkdir ${IUPRED_OUT}/output_files/

mv ${ONE_SEQ_PER_FILE_DIR}/*disord ${IUPRED_OUT}/output_files/

python ${PARSE_SCRIPT_DIR}reformat_iupred_result_file.py ${IUPRED_OUT}/All_iupred_fasta.fas > ${OUT_DIR}/Disordered_results.txt

echo "Done"

cd ${OUT_DIR} :

echo "Calculating low complexiry with segmasker..."

SEGMASKER_OUT=${OUT_DIR}/SEGMASKER_OUTPUT/

if [ ! -d ${SEGMASKER_OUT} ]
then
	mkdir ${SEGMASKER_OUT} ;
	echo "Directory "${SEGMASKER_OUT}" created" ;
else 
	echo "Directory "${SEGMASKER_OUT}" exists" ;
fi

segmasker -in ${PROT} -outfmt fasta -out ${SEGMASKER_OUT}/${JUST_NAME}_segmasker_output.fasta

python ${PARSE_SCRIPT_DIR}reformat_iupred_result_file.py ${SEGMASKER_OUT}/${JUST_NAME}_segmasker_output.fasta > ${OUT_DIR}/Low_comp_results.txt

echo "Done"

echo "Calculating biosynthetic cost..."

python ${PARSE_SCRIPT_DIR}calculate_biosynthetic_cost.py ${PROT} > ${OUT_DIR}/Biosynthetic_cost_results.txt

echo "Done"

echo "Calculating aggregation with tango..."

TANGO_OUT=${OUT_DIR}/TANGO_OUTPUT/

if [ ! -d ${TANGO_OUT} ]
then
	mkdir ${TANGO_OUT} ; mkdir ${TANGO_OUT}/TANGO_SPLIT_FILES/ ; mkdir ${TANGO_OUT}/TANGO_SINGLE_FILES/ ;
	echo "Directory "${TANGO_OUT}" created" ;
else 
	echo "Directory "${TANGO_OUT}" exists" ;
fi

cd ${TANGO_OUT}/TANGO_SPLIT_FILES/ :

python ${PARSE_SCRIPT_DIR}format_fasta_for_tango.py ${PROT} ${JUST_NAME}

cd ${TANGO_OUT}/TANGO_SINGLE_FILES/ :

parallel --GNU 'tango_aggr -inputfile={} &> /dev/null' ::: ../TANGO_SPLIT_FILES/*txt

if [[ $? == 1 ]]
then
	echo "Problem with tango, is the tango executable name correct? It could be tango_aggr instead"
	exit 1
fi

cat *aggregation* | sed '/Sequence/d' >> ${OUT_DIR}/Aggregation_results.txt

echo "Done"

echo "Calculating additional properties with codonW"

cd ${OUT_DIR} :
CODONW_OUTFILE=${OUT_DIR}/${JUST_NAME}"_codonw_output.txt"

CW_OUT=$(codonw ${CDS} ${CODONW_OUTFILE} ${OUT_DIR}/blk -nomenu -nowarn -silent -all_indices 2> /dev/null)

if [[ ${CW_OUT:1:2} != "W" ]]
then 
	echo "Problem with codonw, try running it manually"
	exit 1
fi

#if codonw fails, run it manually using the following line and replace the CODONW_OUTFILE name variable in the next line
#codonw {INFILE} {OUTFILE} ./blk -nomenu -nowarn -silent -all_indices 2> /dev/null
#CODONW_OUTFILE={OUTFILE}

echo "Done"

python ${PARSE_SCRIPT_DIR}Parse_names_build_dataframe_cod.py \
${NAMES_FILE} \
${CAI_OUTPUT_FILE} \
${CODONW_OUTFILE} \
${OUT_DIR}/Aggregation_results.txt \
${OUT_DIR}/${1%\.*}_tmhmm_output_PARSED.txt \
${OUT_DIR}/Low_comp_results.txt \
${OUT_DIR}/Biosynthetic_cost_results.txt \
${OUT_DIR}/Disordered_results.txt \
${OUT_DIR}/Secondary_structure_results.txt \
${OUT_DIR}/ALL_FEATURES_DATAFRAME.csv

echo "File compiled"

if [[ $DEL_DUPS == "no" ]]
then
	cat ${OUT_DIR}/ALL_FEATURES_DATAFRAME.csv | sed 's/\*\**/NA/g' |\
	cut --complement -f 86,87,110,125,128,130,132,134,136,140  |\
	sort -o ${OUT_DIR}/ALL_FEATURES_DATAFRAME_cut.csv
else 
	cat ${OUT_DIR}/ALL_FEATURES_DATAFRAME.csv | sed 's/\*\**/NA/g' |\
	cut --complement -f 86,87,110,125,128,130,132,134,136,140  |\
	sort -u -o ${OUT_DIR}/ALL_FEATURES_DATAFRAME_cut.csv
fi



cat $FIN_HEADER ${OUT_DIR}/ALL_FEATURES_DATAFRAME_cut.csv | sed 's/\*\**/NA/g' > ${OUT_DIR}/final_dataframe.csv

Rscript $PARSE_SCRIPT_DIR"reformat_final_df.r" ${OUT_DIR}/final_dataframe.csv ${OUT_DIR}/final_dataframe_SMALL.csv


