from pathlib import Path
import numpy as np
import csv
import re

def getTime(rtString: str, rtIdx: int) -> float:
    return float(rtString[:rtIdx].strip())

folder_path = Path("../data/")

# Write header to .csv
with open('../data/csv/PerfOutput_Kernel.csv', 'w', newline = '') as csvfile:
    csv_writer = csv.writer(csvfile)
    headerRow = ["N", "Nx", "FP_ARITH", "WALLTIME_MEAN", "WALLTIME_STDDEV"]
    csv_writer.writerow(headerRow)

Npattern = r"_N(\d+)_"
Nxpattern = r"_Nx(\d+)_"

for file in folder_path.iterdir():
    csvArray = []
    if file.is_file():
        print("Processing " + str(file))
        Nmatch = re.search(Npattern, str(file))
        if Nmatch:
            N = Nmatch.group(1)
            csvArray.append(N)
        Nxmatch = re.search(Nxpattern, str(file))
        if Nxmatch:
            Nx = Nxmatch.group(1)
            csvArray.append(Nx)
        # Nidx = str(file).find("N")
        # Nxidx = str(file).find("Nx")
        # N = str(file)[Nidx]
        # Nx = str(file)[Nxidx]
        # csvArray.append(N) 
        # csvArray.append(Nx)
        with open(file, 'r') as fp:
            rtArray = []
            foundArith = False
            while True:
                line = fp.readline()
                if (line == ''):
                    break
                rtIdx = line.find("seconds time elapsed") # magic string because of perf stat output
                arithIdx = line.find("r5301c7") # " 
                if (rtIdx != -1):
                    runtime = getTime(line, rtIdx)
                    rtArray.append(runtime)
                elif (arithIdx != -1 and not foundArith): # number of arithmetic operations performed is always the same for a given point in parameter space
                    # remove all the extra formatting to get only the amount of performed arithmetic
                    base_string = line[:rtIdx].strip() 
                    cleaned_string = base_string.replace(" ", "")
                    numArith = cleaned_string.split("r5301c7")[0] 
                    foundArith = True
                else: 
                    continue
        meanWalltime = np.mean(np.array(rtArray))
        stddevWalltime = np.std(np.array(rtArray))
        csvArray.append(numArith)
        csvArray.append(meanWalltime)
        csvArray.append(stddevWalltime)
        with open('../data/csv/PerfOutput_Kernel.csv', 'a', newline = '') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(csvArray)

            