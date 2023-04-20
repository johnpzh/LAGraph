import os
import sys
import pandas as pd


def main():
    if len(sys.argv) != 2:
        print(F"Usage: python {sys.argv[0]} <LAGraph.output>")
        exit()

    filename = sys.argv[1]
    output = filename + ".trim.out"
    matrices = []
    burkhardt_rt = []
    cohen_rt = []
    sandia_LL_rt = []
    sandia_UU_rt = []
    sandia_LUT_rt = []
    sandia_ULT_rt = []
    burkhardt_nomask_rt = []
    cohen_nomask_rt = []
    sandia_LL_nomask_rt = []
    sandia_UU_nomask_rt = []
    sandia_LUT_nomask_rt = []
    sandia_ULT_nomask_rt = []
    method_index = 0
    with open(filename) as fin, \
            open(output, "w") as fout :
        for line in fin:
            if line.startswith("threads to test:"):
                ## Read the matrix name
                line = fin.readline()
                mtx_file = line.split()[-1]
                matrices.append(os.path.basename(mtx_file))
                method_index = 0
            if line.endswith("presort: 0\n"):
                ## Read runtime
                runtime = float(line.split()[3])
                if method_index == 0:
                    burkhardt_rt.append(runtime)
                    method_index += 1
                elif method_index == 1:
                    cohen_rt.append(runtime)
                    method_index += 1
                elif method_index == 2:
                    sandia_LL_rt.append(runtime)
                    method_index += 1
                elif method_index == 3:
                    sandia_UU_rt.append(runtime)
                    method_index += 1
                elif method_index == 4:
                    sandia_LUT_rt.append(runtime)
                    method_index += 1
                elif method_index == 5:
                    sandia_ULT_rt.append(runtime)
                    method_index += 1
                elif method_index == 6:
                    burkhardt_nomask_rt.append(runtime)
                    method_index += 1
                elif method_index == 7:
                    cohen_nomask_rt.append(runtime)
                    method_index += 1
                elif method_index == 8:
                    sandia_LL_nomask_rt.append(runtime)
                    method_index += 1
                elif method_index == 9:
                    sandia_UU_nomask_rt.append(runtime)
                    method_index += 1
                elif method_index == 10:
                    sandia_LUT_nomask_rt.append(runtime)
                    method_index += 1
                elif method_index == 11:
                    sandia_ULT_nomask_rt.append(runtime)
                    method_index += 1

        ## Save the results
        columns = {
            'matrices': matrices,
            'burkhardt_rt(sec)': burkhardt_rt,
            'cohen_rt(sec)': cohen_rt,
            'sandia_LL_rt(sec)': sandia_LL_rt,
            'sandia_UU_rt(sec)': sandia_UU_rt,
            'sandia_LUT_rt(sec)': sandia_LUT_rt,
            'sandia_ULT_rt(sec)': sandia_ULT_rt,
            'burkhardt_nomask_rt(sec)': burkhardt_nomask_rt,
            'cohen_nomask_rt(sec)': cohen_nomask_rt,
            'sandia_LL_nomask_rt(sec)': sandia_LL_nomask_rt,
            'sandia_UU_nomask_rt(sec)': sandia_UU_nomask_rt,
            'sandia_LUT_nomask_rt(sec)': sandia_LUT_nomask_rt,
            'sandia_ULT_nomask_rt(sec)': sandia_ULT_nomask_rt
        }
        df = pd.DataFrame(data=columns)
        df.to_csv(output, sep='\t', index=False)



if __name__ == "__main__":
    main()
