## This script will take two input files - indel_matrix.nex and Ortho_SNP_matrix.nex - and construct a indel and SNP matrix for input into PAUP



if [ ! $PBS_O_WORKDIR ]
    then
        PBS_O_WORKDIR="$PWD"
fi


cd $PBS_O_WORKDIR

log_eval()
{
  cd $1
  echo -e "\nIn $1\n"
  echo "Running: $2"
  eval "$2"
  status=$?

  if [ ! $status == 0 ]; then
    echo "Previous command returned error: $status"
    exit 1
  fi
}


#test for bedcov and if present run P/A matrix creation


if [ ! -s "indel_matrix.nex" ]; then
    echo -e "\nScript must be supplied with indel_matrix.nex. Please check the directory and analysis and run again\n"
	exit 1
	else
	echo -e "Found indel matrix\n"
fi

if [ ! -s "Ortho_SNP_matrix.nex" ]; then
	echo -e "Script must be supplied with Ortho_SNP_matrix.nex\n"
	echo -e "Please check the directory and analysis and run again\n"
	exit 1
	else
	echo -e "Found Ortho_SNP_matrix.nex\n"
fi

#chomp files into just SNPs and positions
tail -n +8 Ortho_SNP_matrix.nex | head -n -2 > SNP_matrix_tmp

awk ' { for (i=3; i<=NF; i++) {if ($i == $2) $i=0; else $i=1}};  {print $0} ' SNP_matrix_tmp > SNP01.tmp

#indels
tail -n +8 indel_matrix.nex | head -n -2 > indel_matrix_tmp


awk ' { for (i=3; i<=NF; i++) {if ($i == $2) $i=0; else $i=1}};  {print $0} ' indel_matrix_tmp > indel01.tmp



cut -d " " -f 3- SNP01.tmp > SNP01.tmp2
awk '{ print $1 }' SNP01.tmp > SNP.loc

sed -i 's/$/ 0/g' SNP.loc 



cut -d " " -f 3- indel01.tmp > indel01.tmp2
awk '{ print $1 }' indel01.tmp > indel.loc

sed -i 's/$/ 0/g' indel.loc

paste -d ' ' indel.loc indel01.tmp2 > indel.mrg

paste -d ' ' SNP.loc SNP01.tmp2 > SNP.mrg

#cat files and create header

cat SNP.mrg indel.mrg > SNPindel.mrg
x=`cat SNPindel.mrg | wc -l`
y=`head -n 1 SNPindel.mrg | awk '{print NF}'`
z=$((y - 1))
taxa=`tail -n +6 Ortho_SNP_matrix.nex | head -n 1 |cut -d ' ' -f 2-`
grid=`cat SNPindel.mrg`
echo -e "\n#nexus\nbegin data;\ndimensions ntax=$z nchar=$x;\nformat symbols=\"01\" gap=. datatype=standard transpose;\ntaxlabels $taxa\nmatrix\n$grid\n;\nend;" > indel_SNP_matrix.nex

cleanup ()
{

rm SNP_matrix_tmp SNP01.tmp indel_matrix_tmp indel01.tmp SNP01.tmp2 SNP.loc
rm indel01.tmp2 indel.loc indel.mrg SNP.mrg SNPindel.mrg

}

cleanup



exit 0