import os
import sys


def main():
    if len(sys.argv) != 2:
        print(F"Usage: python {sys.argv[0]} <LAGraph.output>")
        exit()

    filename = sys.argv[1]
    output = filename + ".trim.out"
    method_list = []
    time_list = []
    new_method = False
    mtx_name = ''
    with open(filename) as fin, \
            open(output, "w") as fout :
        for line in fin:
            if line.startswith("Reading matrix market file:"):

                ## Save the last method output
                if '' != mtx_name:
                    fout.write(F'\nMatrix: {mtx_name}\n')
                    for mt in method_list:
                        fout.write(F'{mt}\t')
                    fout.write('\n')
                    for t in time_list:
                        fout.write(F'{t}\t')
                    fout.write('\n')
                    method_list.clear()
                    time_list.clear()

                mtx_file = line.split()[-1]
                mtx_name = os.path.basename(mtx_file)

            elif line.startswith("Method:"):
                method = line.split()[1]
                method_list.append(method)
                new_method = True

            elif new_method and line.startswith("nthreads:"):
                time = line.split()[3]
                time_list.append(time)
                new_method = False

        ## Save the last one
        fout.write(F'\nMatrix: {mtx_name}\n')
        for mt in method_list:
            fout.write(F'{mt}\t')
        fout.write('\n')
        for t in time_list:
            fout.write(F'{t}\t')
        fout.write('\n')
        method_list.clear()
        time_list.clear()


if __name__ == "__main__":
    main()
