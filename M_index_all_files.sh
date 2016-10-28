#!/bin/bash
indexes=(1 2 3 4 5 6 7 8 9)
for INDEX in "${indexes[@]}"
do
	echo -e "\n\e[7m$INDEX\e[0m"
	time kallisto index -i func_${INDEX}_index meta_nuc_${INDEX}c.fasta
done
