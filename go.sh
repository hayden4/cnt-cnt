#!/bin/bash

clean.sh

python build_systems.py 0 0

for f in sims/*; do
	cd $f
		qsub potential*.pbs
	cd ../..
done

while [ $(qselect -u hayden4 | wc -l) -gt 0 ]; do
	echo "Number of jobs left: $(qselect -u hayden4 | wc -l)"
	for i in {1..60}; do
		echo -ne "\r\033[KRemaining before next check: $((60-i))"
		sleep 1s
	done
	echo -ne "\r\033[K"
done
echo " "

python gen_table.py
