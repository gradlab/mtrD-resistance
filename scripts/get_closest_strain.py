#!/usr/bin/env python

from collections import defaultdict
import argparse

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Find closest relative to mtrD mutants")
    parser.add_argument("matrix", help="Distances of all mtrD mutants to other strains")
    parser.add_argument("strains", help="TSV with each strain and its mtrD mutation")
    parser.add_argument("metadata", help="metadata file with MICs")
    parser.add_argument("antibiotic", help="anbiotic name for output file")
    parser.add_argument("mic_column", help="column number (0 index) for MIC in metadata file", type = int)
    return parser.parse_args()

def read_distances(matrix_file):
    """Reads a distance matrix created with snp-dists and filtered so that rows are removed if not of interest"""
    with open(matrix_file, "r") as infile:
        distances = defaultdict(dict)
        for i,line in enumerate(infile):
            if i == 0:
                strains = line.strip().split('\t')[1:]
            else:
                line = line.strip().split()
                distances[line[0]] = dict(zip(strains, [int(x) for x in line[1:]]))
    return distances

def get_mtrD_strains(mtrD_file):
    """Reads a TSV list strain names and mtrD mutation"""
    mtrD_strains = {}
    with open(mtrD_file, "r") as infile:
        for line in infile:
            line = line.strip().split('\t')
            mtrD_strains[line[0]] = line[1]
    return mtrD_strains

def get_mic(metadata_file, mic_column):
    """Reads metadata file (specifically strain-table-filtered.csv)"""
    mic = {}
    with open(metadata_file, "r") as infile:
        for line in infile:
            line = line.strip().split(',')
            mic[line[1]] = line[mic_column]
    return mic

def find_closest(distances, mtrD_strains, mic, antibiotic):
    """Finds closest strain that does not also have mtrD mutation and has azithromycin MIC"""
    with open(f"distances/paired_by_distance_{antibiotic}.txt", "w") as outfile:
        outfile.write("mtrD_strain\tmtrD_mutation\tmtrD_mic\twt_strain\twt_mic\tSNP_distance\n")
        for strain in distances:
            distance_dict = distances[strain]
            found = False
            while found == False:
                closest_strain = min(distance_dict, key=distance_dict.get)
                try:
                    if closest_strain not in mtrD_strains.keys() and mic[closest_strain] not in ["NA", ""]:
                        found = True
                    else:
                        del distance_dict[closest_strain]
                except KeyError:
                    del distance_dict[closest_strain]
            outfile.write(f"{strain}\t{mtrD_strains[strain]}\t{mic[strain]}\t{closest_strain}\t{mic[closest_strain]}\t{distance_dict[closest_strain]}\n")

def main():
    args = get_args()
    distances = read_distances(args.matrix)
    mtrD_strains = get_mtrD_strains(args.strains)
    mic = get_mic(args.metadata, args.mic_column)
    find_closest(distances, mtrD_strains, mic, args.antibiotic)

if __name__ == "__main__":
    main()
