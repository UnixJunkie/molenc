#!/bin/bash

~/src/molenc/kbe -i all_uniq_std.txt -k 64 > test_64_1xCPU.txt
sort test_64_1xCPU.txt -o test_64_1xCPU.txt

~/src/molenc/kbe -np 16 -i all_uniq_std.txt -k 64 > test_64_16xCPU.txt
sort test_64_16xCPU.txt -o test_64_16xCPU.txt

diff test_64_1xCPU.txt test_64_16xCPU.txt
