import csv
import glob
import os
import sys

# input should be the main VS_DB directory

args = sys.argv

inputDir = args[1]

# Get all relevant File Paths (different cases if visualization is saved or not)
if os.path.isdir(f"{inputDir}/complexes"):	
	filePaths = glob.glob(f"{inputDir}/complexes/*/*_lddt*_affinity*.sdf")
else:
	filePaths = glob.glob(f"{inputDir}/molecules/VS_DB_*_lddt*_affinity*.sdf")

finalData = []

# For each filepath, get the filename, lddt and affinity
for filePath in filePaths:
	fileName = os.path.basename(filePath).replace("VS_DB_","").replace("_ligand","").split("_rank1")[0]
	lddtScore = filePath.split("_lddt")[-1].split("_affinity")[0]
	affinityScore = filePath.split("_affinity")[-1].split(".sdf")[0]
	finalData.append([fileName,lddtScore,affinityScore,filePath])

# Sort the finalData based on the affinityScore
finalData = sorted(finalData, key=lambda x:(x[2], x[1]),reverse=True)

# Write the sorted dataset to a csv file
with open(f'{inputDir}/summary_results.csv', 'w') as csvfile:
	writer = csv.writer(csvfile, delimiter=';')

	writer.writerow(["Compound_Name","lddt_score","affinity_score","file_path"])

	# Write the data rows
	writer.writerows(finalData)
	
print("Finished summarizing results")
