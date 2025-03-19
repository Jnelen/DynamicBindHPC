#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Relaunch Failed Compounds for DynamicBindHPC

This script checks for failed molecules that did not successfully dock during a DynamicBindHPC job.
It scans the output directory for missing `.sdf` files (screen mode) or missing complex directories (complex mode),
and allows users to relaunch failed jobs.

Usage:
    python relaunchFailedCompounds.py <DynamicBindHPC_run_directory>

Example:
    python relaunchFailedCompounds.py VS_DB_.../

Author:
    Jochem Nelen (jnelen@ucam.edu)
"""

import glob
import os
import re
import shutil
import subprocess
import sys
from typing import List, Dict


def split_list(items: List[str], num_splits: int) -> List[List[str]]:
    """
    Splits a list into `num_splits` roughly equal sublists.

    Args:
        items (List[str]): List of items to split.
        num_splits (int): Number of sublists to create.

    Returns:
        List[List[str]]: A list of lists, each containing a subset of the input items.
    """
    if num_splits > len(items):
        print("More jobs than files, launching 1 job per file")
        return [[item] for item in items]

    k, m = divmod(len(items), num_splits)
    return [
        items[i * k + min(i, m) : (i + 1) * k + min(i + 1, m)]
        for i in range(num_splits)
    ]


def get_all_molecules(input_dir: str) -> Dict[str, str]:
    """
    Reads CSV files from the input directory to determine which molecules were used as input.

    Args:
        input_dir (str): Path to the DynamicBindHPC run directory.

    Returns:
        Dict[str, str]: A dictionary where keys are molecule names and values are their paths.
    """
    path_dict = {}
    csv_files = glob.glob(f"{input_dir}/csvs/*.csv")

    for path in csv_files:
        with open(path, "r") as file:
            lines = file.readlines()
            for line in lines[1:]:  # Skip header
                line_split = line.strip().split(";")
                mol_name, mol_path = line_split[0], line_split[-1]
                path_dict[mol_name] = mol_path

    return path_dict


def get_finished_items(input_dir: str) -> List[str]:
    """
    Identifies successfully processed molecules, either from 'molecules/' (screen mode)
    or 'complexes/' (complex mode), based on the existing directory structure.

    Args:
        input_dir (str): Path to the DynamicBindHPC run directory.

    Returns:
        List[str]: A list of successfully processed molecule names.

    Raises:
        SystemExit: If neither 'molecules/' nor 'complexes/' directories exist.
    """
    molecules_path = f"{input_dir}/molecules/"
    complexes_path = f"{input_dir}/complexes/"

    if os.path.isdir(molecules_path):
        print("Detected screen mode: Checking SDF files in 'molecules/'")
        return [
            os.path.basename(path).split("VS_DB_")[-1].split("_rank")[0]
            for path in glob.glob(f"{molecules_path}/*.sdf")
        ]

    elif os.path.isdir(complexes_path):
        print("Detected complex mode: Checking directories in 'complexes/'")
        return [
            os.path.basename(path.rstrip("/"))  # Remove trailing slash
            for path in glob.glob(f"{complexes_path}/*")  # Get directory names
            if os.path.isdir(path)
        ]  # Ensure it's a directory

    else:
        sys.exit(
            "ERROR: Neither 'molecules/' nor 'complexes/' directories exist. Invalid DynamicBindHPC run."
        )


def clean_output_directory(redo_dir: str) -> None:
    """Checks if the redo directory exists and prompts the user for cleanup."""
    if os.path.isdir(redo_dir):
        print(
            f"The directory {redo_dir} already exists. To continue, delete it or choose another output directory."
        )
        answer = input("Do you want to remove it? (y/n) ").lower()

        while answer not in ("y", "n", "yes", "no"):
            print("Invalid input. Please enter y(es) or n(o).")
            answer = input("Do you want to remove or overwrite it? (y/n) ").lower()

        if answer in ("y", "yes"):
            shutil.rmtree(redo_dir, ignore_errors=False)
        else:
            sys.exit("Exiting without relaunching jobs.")


def create_redo_directory(redo_dir: str) -> None:
    """
    Creates necessary subdirectories for relaunching jobs.

    Args:
        redo_dir (str): Path to the redo directory.
    """
    os.mkdir(redo_dir)
    os.mkdir(f"{redo_dir}/csvs/")
    os.mkdir(f"{redo_dir}/jobs_out/")
    os.mkdir(f"{redo_dir}/jobs/")


def relaunch_jobs(
    input_dir: str, ligand_paths_split: List[List[str]], redo_dir: str
) -> None:
    """Relaunches failed docking jobs using the original job settings."""
    job_paths = glob.glob(f"{input_dir}/jobs/job_*.sh")

    if not job_paths:
        sys.exit("No job files found in the original directory.")

    with open(job_paths[0], "r") as job_file:
        job_template = job_file.readlines()[1]

    protein_path = glob.glob(f"{input_dir}/*.pdb")[0]

    for i, job_ligands in enumerate(ligand_paths_split):
        csv_file_path = f"{redo_dir}/csvs/job_csv_{i+1}.csv"

        with open(csv_file_path, "w") as job_csv:
            job_csv.write("name;protein_path;ligand\n")
            for ligand in job_ligands:
                complex_name = os.path.basename(ligand).split(".")[0]
                job_csv.write(f"{complex_name};{protein_path};{ligand}\n")

        job_cmd = re.sub(
            r"/jobs_out/job_.*\.out",
            f"/redo/jobs_out/redo_job_{i}_%j.out",
            job_template,
        )
        job_cmd = re.sub(
            r"(--protein_ligand_csv\s+)[^ ]+(\s+--samples_per_complex)",
            rf"\1{csv_file_path}\2",
            job_cmd,
        )

        job_file_path = f"{redo_dir}/jobs/redo_job_{i+1}.sh"
        with open(job_file_path, "w") as job_file:
            job_file.write("#!/usr/bin/env bash\n")
            job_file.write(job_cmd)

        subprocess.run(job_cmd, shell=True)


def main() -> None:
    """Main function to process failed molecules and relaunch jobs."""
    if len(sys.argv) < 2:
        sys.exit("You must provide a DynamicBindHPC run directory as an argument.")

    input_path = sys.argv[1]
    if not os.path.isdir(input_path):
        sys.exit("The input path is not a valid directory.")

    path_dict = get_all_molecules(input_path)
    finished_list = get_finished_items(input_path)

    total_count = len(path_dict)
    success_count = len(finished_list)

    for finished in finished_list:
        path_dict.pop(finished, None)  # removing completed molecules

    failed_count = len(path_dict)

    print(f"Total compounds: {total_count}")
    print(
        f"Successfully processed: {success_count}/{total_count} ({(success_count / total_count) * 100:.1f}%)"
    )
    print(
        f"Failed to process: {failed_count}/{total_count} ({(failed_count / total_count) * 100:.1f}%)"
    )

    if failed_count == 0:
        print("All compounds processed successfully! No jobs to relaunch.")
        return

    while True:
        answer = input("How many jobs should be launched? ").strip()
        if answer.isdigit() and int(answer) > 0:
            job_number = int(answer)
            break
        elif answer.lower() in ("n", "no"):
            print("Exiting without launching any jobs.")
            sys.exit(0)

        else:
            print("Invalid input. Please enter a positive integer greater than 0.")

    ligand_paths_split = split_list(list(path_dict.values()), job_number)
    redo_directory = f"{input_path}/redo/"
    clean_output_directory(redo_directory)
    create_redo_directory(redo_directory)

    relaunch_jobs(input_path, ligand_paths_split, redo_directory)


if __name__ == "__main__":
    main()
