#!/bin/bash

echo "Number range start: "
read start

echo "Number range end: "
read end

echo "Move? (1 for yes)"
read move

echo "Delete? (1 for yes)"
read delete

if [ $move -eq 1 ];
then
for ((i=$start; i<=$end; i++))
do

    mkdir ./gaussian_filter${i}_reduced
    mv ./gaussian_filter${i}/*0000.npy ./gaussian_filter${i}_reduced/
    mv ./gaussian_filter${i}/*5000.npy ./gaussian_filter${i}_reduced/
    mv ./gaussian_filter${i}/gaussian_filter${i}_plots/ ./gaussian_filter${i}_reduced/

done
fi


if [ $delete -eq 1 ];
then
for ((i=$start; i<=$end; i++))
do
    rm -r ./gaussian_filter${i}
done
fi