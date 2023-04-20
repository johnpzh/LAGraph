
matrices=("bcsstk17" \
          "pdb1HYS" \
          "rma10" \
          "cant" \
          "consph" \
          "pwtk" \
          "shipsec1" \
          "cop20k_A" \
          "scircuit" \
          "mac_econ_fwd500")

for mtx in "${matrices[@]}"; do
  echo ""
  echo "Running ${mtx}..."
  ../../cmake-build-debug/src/benchmark/bfs_burble_demo ~/local/react-eval/matrices/${mtx}/${mtx}.mtx &> burble.fixed_saxpy.${mtx}.txt
done

echo "Done."