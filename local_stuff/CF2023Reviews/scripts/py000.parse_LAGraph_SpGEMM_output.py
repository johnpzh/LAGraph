import os
import sys
import pandas as pd


def main():
    if len(sys.argv) != 2:
        print(F"Usage: python {sys.argv[0]} <LAGraph.output>")
        exit()

    filename = sys.argv[1]
    output = filename + ".trim.out"
    masking_time = []
    no_masking_time = []
    matrices = []
    new_mtx = True
    with open(filename) as fin, \
            open(output, "w") as fout :
        for line in fin:
            if line.startswith("threads to test:"):
                ## Read the matrix name
                line = fin.readline()
                mtx_file = line.split()[-1]
                matrices.append(os.path.basename(mtx_file))
                new_mtx = True

            elif line.endswith("presort: 0\n"):
                # print(F'endswith line: {line}')
                ## Read the runtime
                runtime = float(line.split()[3])
                if new_mtx:
                    no_masking_time.append(runtime)
                    new_mtx = False
                else:
                    masking_time.append(runtime)
        
        ## Save the results
        # print(matrices)
        # print(no_masking_time)
        # print(masking_time)
        columns = {
             'matrices': matrices,
             'no_masking_runtime(sec)': no_masking_time,
             'masking_runtime(sec)': masking_time
        }
        # fout.write(pd.DataFrame(data=columns))
        df = pd.DataFrame(data=columns)
        df.to_csv(output, sep='\t', index=False)


if __name__ == "__main__":
    main()
