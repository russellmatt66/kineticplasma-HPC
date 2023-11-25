import subprocess
import os
import sys
import math

def execute_bash_script(bash_script_path: str, executable_path: str, num_runs: str, output_location: str) -> None:
    try:
        subprocess.run(['bash', bash_script_path, executable_path, num_runs, output_location], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing Bash script: {e}")

# Modify N and Nx in *.inp
def modify_inp_file(inp_file_path: str, parameter_values: tuple) -> int:
    # Logic goes here
    print("Modifying input file: " + inp_file_path)
    N_new = parameter_values[0]
    Nx_new = parameter_values[1]
    # Parse file: Change N and Nx, get Nt for providing correct output file name
    with open(inp_file_path, 'r') as fp:
        lines = fp.readlines()

    for i, line in enumerate(lines):
        if line.startswith('N='):
            lines[i] = f'N={N_new}\n'
        elif line.startswith('Nx='):
            lines[i] = f'Nx={Nx_new}\n'
        elif line.startswith('Nt='):
            Nt = int(line.split('=')[1])
    
    with open(inp_file_path, 'w') as fp:
        fp.writelines(lines)

    return Nt


# Scan through N and Nx by running ./executable a set number of times via calling perf_data.sh, then changing the values in *.inp
def main():
    #  Logic goes here
    executable_path = sys.argv[1]
    bash_script_path = sys.argv[2]
    input_file_path = sys.argv[3]
    output_folder = sys.argv[4]

    N_min = int(sys.argv[5])
    N_max = int(sys.argv[6])
    Nx_min = int(sys.argv[7])
    Nx_max = int(sys.argv[8])

    # Should look like [(N_min, Nx_min), 2*(N_min, Nx_min), 4*(N_min, Nx_min), ... , (N_max, Nx_max)]
    parameter_list = [(2**i, 2**j) for i in range(int(math.log2(N_min)), int(math.log2(N_max))+1)
                      for j in range(int(math.log2(Nx_min)), int(math.log2(Nx_max))+1)]

    num_runs = '25'

    for (N,Nx) in parameter_list:
        Nt_num = modify_inp_file(input_file_path, (N,Nx))
        output_location = output_folder + "PerfStat_N" + str(N) + "_Nx" + str(Nx) + "_Nt" + str(Nt_num) + ".txt" 
        execute_bash_script(bash_script_path, executable_path, num_runs, output_location)

if __name__ == "__main__":
    main()
