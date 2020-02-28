#!/bin/bash

#rename bedcov because I'm too studpid and can't figure out how to keep file names when collecting multiple outputs into a single process and the Nextflow documentation sucks for this too

for f in *.bedcov; do 
  new_name=$(basename $(readlink -f $f))
  mv "$f" "$new_name"
done

array=($(find *.bedcov -printf "%f "))
array2=("${array[@]/.bedcov/}")
array3=("${array[@]/.bedcov/.cut}")
n=${#array2[@]}
for (( i=0; i<n; i++ )); do
  cat ${array2[i]}.bedcov | sort -k1,1 -k2,2n > ${array2[i]}.sort
  cat ${array2[i]}.sort | cut -f7 > ${array2[i]}.cut
  sed -i '1s/^/'${array2[i]}'\n/' ${array2[i]}.cut
done

cat "${array2[1]}".sort | cut -f1-3 > head.bed
sed -i '1s/^/\t\t\n/' head.bed
paste head.bed ${array3[*]} > Bedcov_merge.txt


exit 0
  