#!/bin/bash

#########################################################################
# The following script will check to see if the SnpEff database exists and if not it will attempt to download and install
#
# Written by D. Sarovich
# dsarovich@usc.edu.au
#
#########################################################################

database=$1
baseDir=$2
ref=$3

#determine snpEff location
snpeff_dir=$(dirname $(readlink -f $(which snpEff)))

snpEff databases | grep -w "${database}" > /dev/null
status=$?
if [ ! "$status" == 0 ]; then
    echo "SPANDx couldn't find the annotation database in the SnpEff"
	echo "The name of the annotated reference genome specified with the --database and must match a reference genome in the SnpEff database"
    echo "Does the SnpEff.config file contain the reference specified with the --database switch?"
	echo "Is the SnpEff.config file in the location specified by SPANDx.config?"
	echo "If both of these parameters are correct please refer to the SnpEff manual for further details on the setup of SnpEff"
	exit 1
else
    echo -e "SPANDx found the reference file in the SnpEff database list"
    echo -e "Looking for the resource file"
	snpEff dump "${database}" > /dev/null
	status=$?
    if [ ! "$status" == 0 ]; then
	    echo -e "SPANDx couldn't find the snpEff database in the default location and failed automatic download"
		echo "Attempting manual download"
		dl_link=$(snpEff databases | grep -w "${database}" | awk '{print $4}')
		wget ${dl_link}
		dl_link_basename=$(basename ${dl_link})
		unzip $dl_link_basename
		mv ./home/pcingola/snpEff/data/* ${snpeff_dir}/data
		snpEff dump "${database}" > /dev/null
		status=$?
        if [ ! "$status" == 0 ]; then
		  echo "Couldn't manually download and install snpEff database. Please do so manually"
		  echo "snpEff directory is here ${snpeff_dir}"
		  echo "The database I'm looking for is here ${database}"
		  echo "The download link I'm attempting to grab is ${dl_link}"
		else
		  echo "The snpEff database has been successfully installed"
		fi		  
	else
	  echo "SPANDx found the reference in the default location"
    fi
fi

chr_name=$(snpEff dump "${database}" | grep -A1 'Chromosomes names' | tail -n1 | awk '{print $2}'|sed "s/'//g")
ref_chr_name=$(head -n1 "$ref".fasta | sed 's/>//')
if [ "$chr_name" == "$ref_chr_name" ]; then
	echo -e "Chromosome names in the SnpEff database match the reference chromosome names, good\n"
else
    echo -e "Chromosome names in the SnpEff database DON'T match the reference chromosome names.\n"
	echo -e "Please change the names of the reference file to match those in the SnpEff database.\n"
	echo -e "If you are unsure what these are, run: snpEff dump ${database}\n"
	echo -e "The first chromosome name is $chr_name.\n\n"
    echo -e "If you choose to continue the annotation component of SPANDx may fail.\n"
    exit 1
fi

exit 0