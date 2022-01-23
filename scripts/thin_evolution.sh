#!/bin/bash

awk -v lastt=0 '{thist=$1; if(thist > lastt + 1e-3 || thist > lastt * 10.0**(4e-4)) {lastt = thist; print;};}' $1
