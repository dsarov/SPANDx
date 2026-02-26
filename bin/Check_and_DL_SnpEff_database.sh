#!/bin/bash

#########################################################################
# The following script will check to see if the SnpEff database exists and if not it will attempt to download and install
#
# Written by D. Sarovich
# dereksarovich@gmail.com
#
#########################################################################
database=$1
baseDir=$2
ref=$3
SNPEFF_DATA_DIR=$4 
SNPEFF_CONFIG=$5

echo "We are using snpEff that is located here"
which snpEff


echo "Running: snpEff databases -dataDir ${SNPEFF_DATA_DIR} -c ${SNPEFF_CONFIG} | grep -w ${database} > /dev/null"

snpEff databases -dataDir "${SNPEFF_DATA_DIR}" -c "${SNPEFF_CONFIG}" | grep -w "${database}" > /dev/null
status=$?
if [ ! "$status" == 0 ]; then
    echo "Couldn't find the snpeEff database in the cloud. Is this a standard database? If not you will need to manually add this to snpEff. And change the location for the snpeff_cache variable in the nextflow.config file."
	#exit 1
else
    echo -e "SPANDx found the reference file in the SnpEff database list"
    echo -e "Looking for the resource file"
	echo "Running: snpEff dump -dataDir ${SNPEFF_DATA_DIR} -c ${SNPEFF_CONFIG} ${database} > /dev/null"
	snpEff dump -dataDir "${SNPEFF_DATA_DIR}" -c ${SNPEFF_CONFIG} "${database}" > /dev/null 
	status=$?
    if [ ! "$status" == 0 ]; then
	    echo -e "SPANDx couldn't find the snpEff database, attempting automatic download"
		$(which snpEff) download -dataDir "${SNPEFF_DATA_DIR}" -v "${database}"
		status=$?
        if [ ! "$status" == 0 ]; then
		  echo "Couldn't download and install snpEff database. Please do so manually"
		  exit 1
		else
		  echo "The snpEff database has been successfully installed"
		fi
	else
	  echo "SPANDx found the reference in the default location"
    fi
fi

chr_name=$(snpEff dump -dataDir "${SNPEFF_DATA_DIR}" -c ${SNPEFF_CONFIG} "${database}" | grep -A1 'Chromosomes' | tail -n1 | awk '{print $2}'|sed "s/'//g")

chr_name=$(snpEff dump "${database}" -c ${SNPEFF_CONFIG} | grep -A1 'Chromosomes' | tail -n1 | awk '{print $2}'|sed "s/'//g")
ref_chr_name=$(head -n1 "$ref" | sed 's/>//')
if [ "$chr_name" == "$ref_chr_name" ]; then
	echo -e "Chromosome names in the SnpEff database match the reference chromosome names, good\n"
else
    echo -e "Chromosome names in the SnpEff database DON'T match the reference chromosome names.\n"
	echo -e "Please change the names of the reference file to match those in the SnpEff database.\n"
	echo -e "If you are unsure what these are, run: snpEff dump ${database}\n"
	echo -e "The first chromosome name is $chr_name.\n\n"
    echo -e "If you choose to continue the annotation component of SPANDx may fail.\n"
    #exit 1
fi

exit 0
