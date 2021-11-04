#!/bin/bash

for a in `cat $1`; do
    grep $a athV.ann.vcf  > geneToPat/$a
done
