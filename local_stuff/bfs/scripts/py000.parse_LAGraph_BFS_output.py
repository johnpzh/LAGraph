import os
import sys


def main():
    if len(sys.argv) != 2:
        print(F"Usage: python {sys.argv[0]} <LAGraph.output>")
        exit()

    filename = sys.argv[1]
    output = filename + ".trim.out"
    with open(filename) as fin, \
            open(output, "w") as fout :
        for line in fin:
            if line.startswith("Avg:"):
                columns = line.split()
                time = columns[7]
                mtx_name = os.path.basename(columns[9])
                fout.write(F'{mtx_name}\t{time}\n')


if __name__ == "__main__":
    main()
