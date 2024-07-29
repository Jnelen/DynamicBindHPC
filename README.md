# DynamicBindHPC
DynamicBindHPC is a fork of [DynamicBind](https://github.com/luwei0917/DynamicBind), which adds enables you to run DynamicBind on HPC systems using Singularity and Slurm.    
For more details about DynamicBind itself, we refer to the [DynamicBind Github](https://github.com/luwei0917/DynamicBind) and the original [publication](https://www.nature.com/articles/s41467-024-45461-2).

### Requirements:
* Singularity 
* Slurm (There is a --no_slurm mode, but using Slurm is highly recommended)

### Installation instructions:
1. Clone the repository and navigate to it
    ```
    git clone https://github.com/Jnelen/DynamicBindHPC
    ```
   ```
   cd DynamicBindHPC
   ```
   
2. Run a test example to be prompted to automatically download the Singularity image (~5 GB) and again to download and install the model weights (~450 MB). Additionally, required cache look-up tables for SO(2) and SO(3) distributions will have to be generated. (This only needs to happen once and usually takes around 15 minutes).  
   The `--no_slurm` flag is optional here, but makes it easier to track the progress.   
   ```
   python inferenceVS.py -p data/origin-1qg8.pdb -l data/ -out TEST -j 1 --no_slurm
   ```  
   Or if you have access to a GPU, you can also add the -gpu tag like this:  
   ```
   python inferenceVS.py -p data/origin-1qg8.pdb -l data/ -out TEST -j 1 -gpu --no_slurm
   ```  
You can also download the Singularity image manually:
   ```
   wget --no-check-certificate -r "https://drive.usercontent.google.com/download?id=1QFfsqvEdDJVLUvh0kVhksG8sIsrGenWO&confirm=t" -O singularity/DynamicBindHPC.sif
   ```
   
   alternatively, you can build the singularity image yourself using:
   ```
   singularity build singularity/DynamicBindHPC.sif singularity/DynamicBindHPC.def
   ```
### Options

The main file to use is `inferenceVS.py`. It has the following options/flags:  

- `-p`, `-r`, `--protein_path`: 
  Path to the protein/receptor `.pdb` file.

- `-l`, `--ligand`: 
  The path to the directory of (separate) `mol2`/`sdf` ligand files.

- `-o`, `--out`, `--out_dir`: 
  Directory where the output structures will be saved to.

- `-j`, `--jobs`: 
  Number of jobs to use.

- `-qu`, `--queue`: 
  On which node to launch the slurm jobs. The default value is the default queue for the user. Might need to be specified if there is no default queue configured.

- `-m`, `--mem`: 
  How much memory to use for each job. The default value is `4GB`.

- `-gpu`, `--gpu`: 
  Use GPU resources. This will accelerate docking calculations if a compatible GPU is available.

- `-c`, `--cores`: 
  How many cores to use for each job. The default value is `1` when used with the GPU option enabled, otherwise it defaults to `4` cores.

- `-n`, `--num_outputs`, `--samples_per_complex`: 
  How many structures to output per compound. The default value is `1`.
  
- `--model`: 
  DynamicBind supports two models: `ema_inference_epoch314_model.pt` and `pro_ema_inference_epoch138_model.pt`.  The default model is the same as the one used in the paper (`ema_inference_epoch314_model.pt`)
  
- `--remove_hs`: 
  Remove the hydrogens in the final output structures.
  
- `--no_slurm`: 
  Don't use slurm to handle the resources. This will run all samples in interactive mode. The `--gpu` and `-c` options will still work to use a gpu and set the number of CPU cores. However, other Slurm arguments such as the amount memory, time limit, ... will be ignored.
  
- `--no_clean`: 
  Don't clean the input protein structure. Not recommended unless you properly prepared the input protein structure (removed ligands, waters, ...)
  
- `--save_visualisation`: 
  Save the output ligand and protein files. These files can be used to generate an animation with the `movie_generation.py` script.
  
- `-h`, `--help`: 
  Show the help message and exit.

### Example
To run a DynamicBindHPC calculation
```
python inferenceVS.py -p data/origin-1qg8.pdb -l data/ -out TEST -j 1 -gpu --no_slurm --save_visualisation
```

And to generate the different steps from the diffusion process you can do:

```
python movie_generation.py VS_DB_TEST_.../complexes/1opj_STI_A/
```

It is also possible to run this with a gpu, which speeds up the calculations significantly. You can specify a device as well (which sets the `CUDA_VISIBLE_DEVICES` to your input):  
    
```
python movie_generation.py outputDir/complexes/complex_of_interest/ --device 1 --gpu  
```

A `rank1_receptor_reverseprocess_relaxed.pdb` and `rank1_ligand_reverseprocess_relaxed.sdf` will be output in the same directory as the original input. These contain the states of the different diffusion steps.  

Note: I found that often in the first steps of the animation, the poses can clash with the protein backbone. However, I run into the exact same problems when running native DynamicBind with the same inputs. 
I raised an issue about this on their GitHub, and if a fix is made, I will also try to update DynamicBindHPC accordingly.

## License
MIT

