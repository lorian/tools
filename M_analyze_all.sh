#!/bin/bash
SWITCH=$1
if [ $SWITCH = "missing" ]; then
	outputs=(i100_nolisteria_abundance.txt labeled_i100_nolisteria.kraken i100_nobacillus_abundance.txt labeled_i100_nobacillus.kraken)
elif [ $SWITCH = "martin" ]; then
	outputs=(i100_Martin_ralstonia.clark i100_Martin_ralstonia_abundance.txt labeled_i100_Martin_ralstonia.kraken)
elif [ $SWITCH = "i100" ]; then
	outputs=(gasic_output_100.txt testJ2_results.xprs i100_ralstoniaplasmid.clark ralstoniaplasmid_abundance.txt labeled_i100_ralstoniaplasmid.kraken)
else
	outputs=(gasic_output_100.txt testJ2_results.xprs i100_ralstoniaplasmid.clark ralstoniaplasmid_abundance.txt labeled_i100_ralstoniaplasmid.kraken i100_Martin_ralstonia.clark i100_Martin_ralstonia_abundance.txt labeled_i100_Martin_ralstonia.kraken)
fi
for OUTPUT in "${outputs[@]}"
do
	echo -e "\n\e[7m$OUTPUT\e[0m"
	python ~/compbio_tools/M_compare_metagenome_results_to_truth.py -s $OUTPUT
done
