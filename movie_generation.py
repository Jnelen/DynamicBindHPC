import os
import sys
import subprocess
import time
def do(cmd, get=False, show=True):
    if get:
        out = subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True).communicate()[0].decode()
        if show:
            print(out, end="")
        return out
    else:
        return subprocess.Popen(cmd, shell=True).wait()

file_path = os.path.realpath(__file__)
script_folder = os.path.dirname(file_path)
# print(file_path, script_folder)


import argparse
parser = argparse.ArgumentParser(description="/mnt/nas/glx-share-cache/InfraDev/glx-schrodinger/envs/relax/bin/python movie_generation.py results/s21_3/index0_idx_0/ 1+2 --device 3")

parser.add_argument('prediction_result_path', type=str, default='results/test/index0_idx_0', help='informative name used to name result folder')
parser.add_argument('--rank', type=str, default="1", help='specify the sample to generate movie.\
             (the samples are sorted by their confidence, with rank 1 being considered the best prediction by the model, rank 40 the worst), \
             could give multiple. for example 1+2+3')
parser.add_argument('--device', type=int, default=0, help='CUDA_VISIBLE_DEVICES')
parser.add_argument('--inference_steps', type=int, default=20, help='num of coordinate updates. (movie frames)')
parser.add_argument('--remove_hs', action='store_true', default=False, help='Remove the hydrogens in the final output structures')
parser.add_argument('--debug', action='store_true', default=False, help='Be more verbose about the commands that are executed')
parser.add_argument('--gpu', action='store_true', default=False, help='Accelerate movie generation by using a GPU')
parser.add_argument('--num_workers', default=4, help='How many CPU cores to use. Putting a negative number enables usage of all cores')

args = parser.parse_args()

startTime = time.time()

if args.remove_hs:
    remove_hs_arg = "--remove_hs"
else:
    remove_hs_arg = f""

if args.num_workers is None or args.num_workers < 0:
    args.num_workers = os.cpu_count()
    
if args.gpu:
    gpu_arg_singularity = "--nv"
    gpu_arg = "--gpu"
    
    ## Make sure the correct GPU is being used
    os.environ["CUDA_VISIBLE_DEVICES"]= str(args.device)
else:
    gpu_arg_singularity = ""
    gpu_arg = ""


## Getting the poses from the pkl file
for rank in args.rank.split("+"):
    cmd = f"singularity exec {gpu_arg_singularity} --bind $PWD singularity/DynamicBindHPC.sif python3 -u {script_folder}/save_reverseprocess.py --pklFile {args.prediction_result_path}/rank{rank}_reverseprocess_data_list.pkl {remove_hs_arg}"
    print(f"Extracting poses from pkl ({args.prediction_result_path}/rank{rank}_reverseprocess_data_list.pkl)")
    if args.debug:
        print(cmd)
    do(cmd)

## Check if relaxed final ligand is present, if not first run relax final on final ligand-protein complex before running relax_vis.py
if not "_relaxed.sdf" in os.listdir(args.prediction_result_path):
    cmd = f"singularity exec {gpu_arg_singularity} --bind $PWD singularity/DynamicBindHPC.sif python3 -u relax_final.py --num_workers {args.num_workers} {gpu_arg} --input_paths {args.prediction_result_path}"
    print("performing relaxation of final frame")
    if args.debug:
        print(cmd)
    do(cmd)
    
## Add singularity to command
cmd = f"singularity exec {gpu_arg_singularity} --bind $PWD singularity/DynamicBindHPC.sif python3 -u {script_folder}/relax_vis.py --rank {args.rank} --prediction_result_path {args.prediction_result_path} --inference_steps {args.inference_steps}"
print("Processing movie frames..")
if args.debug:
    print(cmd)
do(cmd)

print(f"Saved the results in {args.prediction_result_path}")
print(f"Finished after {time.time()-startTime:.2f} seconds")