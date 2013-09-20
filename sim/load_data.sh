#!/bin/bash

cat ../data/temp_series_degreesC_0_65MYA_by_0.001MY.csv | tr "\n" "," | sed 's/^/double temps[]={/' | sed 's/,$/};/' > data.cpp
