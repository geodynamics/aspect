#!/bin/bash

python testdata.py
cat header.txt > testdata.txt
cat data.txt >> testdata.txt

rm data.txt
