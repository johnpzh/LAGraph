
bin="../cmake-build-debug/src/benchmark/mxm_validate_demo"

matrix_root="/Users/peng599/local/react-eval/matrices"

mask_root="/Users/peng599/pppp/CLion/COMET_masking/local_stuff/masked_spgemm/data"



#data="bcsstk29"
#data="shipsec1"
data="cant"

eval ${bin} "${matrix_root}/${data}/${data}.mtx" "${mask_root}/${data}.ones.mtx"

