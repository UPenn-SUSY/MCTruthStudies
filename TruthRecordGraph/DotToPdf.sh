#!/bin/bash

dot_file=$1
eps_file=$(echo $1 | sed 's#\.dot#.eps#g')
pdf_file=$(echo $1 | sed 's#\.dot#.pdf#g')

echo "Converting ${dot_file} to ${eps_file} then to ${pdf_file}"

dot -Tps2 $dot_file -o $eps_file
ps2pdf $eps_file $pdf_file

