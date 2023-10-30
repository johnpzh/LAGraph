#!/bin/bash

# LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
# SPDX-License-Identifier: BSD-2-Clause
# See additional acknowledgments in the LICENSE file,
# or contact permission@sei.cmu.edu for the full terms.

# do_gap_tc: run the GAP benchmarks for Triangle Counting

# Usage:
#
#   To run using *.mtx format, with the files in ../../../GAP:
#
#       ./do_gap_tc > myoutput.txt
#
#   To run using binary *.grb format, with the files in ../../../GAP:
#
#       ./do_gap_tc grb > myoutput.txt
#
#   To run using binary *.mtx format, with the files in /my/stuff/GAP
#
#       ./do_gap_tc mtx /my/stuff/GAP > myoutput.txt
#
#   To run using binary *.mtx format, with the files in /my/stuff/GAP
#
#       ./do_gap_tc grb /my/stuff/GAP > myoutput.txt

echo " "
echo "======================================================================"
echo "GAP benchmarks using LAGraph+GraphBLAS: Mask SpGEMM"
echo "======================================================================"

#if [ -z "$1" ]; then KIND="mtx" ; else KIND=$1 ; fi
#echo "Matrix input file format: " $KIND
#
#if [ -z "$2" ]; then GAP="../../../GAP" ; else GAP=$2 ; fi
#echo "GAP matrices located in:  " $GAP
#
##../../build/src/benchmark/tc_demo $GAP/GAP-kron/GAP-kron.$KIND
##../../build/src/benchmark/tc_demo $GAP/GAP-urand/GAP-urand.$KIND
##../../build/src/benchmark/tc_demo $GAP/GAP-twitter/GAP-twitter.$KIND
##../../build/src/benchmark/tc_demo $GAP/GAP-web/GAP-web.$KIND
#../../build/src/benchmark/mxm_demo $GAP/GAP-road/GAP-road.$KIND

#../../build/src/benchmark/mxm_demo \
#  /qfs/people/peng599/data/GAP/GAP-road/GAP-road.mtx \
#  /qfs/people/peng599/pppp/clion/comet_ReACT-fusion-IT_mac/local_stuff/masked_spgemm/data/bcsstk17.zeros.mtx
#../../build/src/benchmark/mxm_demo \
#  /qfs/people/peng599/pppp/clion/react-eval_mac/matrices/bcsstk17/bcsstk17.mtx \
#  /qfs/people/peng599/pppp/clion/react-eval_mac/matrices/bcsstk17/bcsstk17.mtx
#../../build/src/benchmark/mxm_demo \
#  /qfs/people/peng599/pppp/clion/react-eval_mac/matrices/bcsstk17/bcsstk17.mtx \
#  /qfs/people/peng599/pppp/clion/comet_ReACT-fusion-IT_mac/local_stuff/masked_spgemm/data/bcsstk17.zeros.mtx

matrices=("com-LiveJournal" "com-Orkut" "GAP-web")

for mtx in "${matrices[@]}"; do

  ../../build/src/benchmark/mxm_demo \
    "/qfs/people/peng599/data/${mtx}/${mtx}.mtx"

#  ../../build/src/benchmark/mxm_demo \
#    "/qfs/people/peng599/data/${mtx}/${mtx}.mtx" \
#    "/qfs/people/peng599/data/${mtx}/${mtx}.mtx"

#  ../../build/src/benchmark/mxm_demo \
#    "/qfs/people/peng599/pppp/clion/react-eval_mac/matrices/${mtx}/${mtx}.mtx" \
#    "/qfs/people/peng599/pppp/clion/comet_ReACT-fusion-IT_mac/local_stuff/masked_spgemm/data/${mtx}.zeros.mtx"
done

#../../build/src/benchmark/mxm_demo \
#  /qfs/people/peng599/pppp/clion/react-eval_mac/matrices/bcsstk29/bcsstk29.mtx \
#  /qfs/people/peng599/pppp/clion/comet_ReACT-fusion-IT_mac/local_stuff/masked_spgemm/data/bcsstk29.zeros.mtx
