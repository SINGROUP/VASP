#!/bin/bash

str='00 01 02 03 04'

for i in $str;
do
    echo $i
    povray neb_tmp_$i.ini
done
