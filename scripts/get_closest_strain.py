#!/usr/bin/env python

from collections import defaultdict
import argparse

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Find closest relative to mtrD mutants")
    parser.add_argument("matrix", help="Distances of all mtrD mutants to other strains")
    parser.add_argument("strains", help="TSV with each strain and its mtrD mutation")
    parser.add_argument("metadata", help="metadata file with MICs")
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

def get_azi(metadata_file):
    """Reads metadata file (specifically strain-table-filtered.csv)"""
    azi_mic = {}
    with open(metadata_file, "r") as infile:
        for line in infile:
            line = line.strip().split(',')
            azi_mic[line[1]] = line[9]
    return azi_mic

def find_closest(distances, mtrD_strains, azi_mic):
    """Finds closest strain that does not also have mtrD mutation and has azithromycin MIC"""
    with open("distances/paired_by_distance.txt", "w") as outfile:
        outfile.write("mtrD_strain\tmtrD_mutation\tmtrD_azi_mic\twt_strain\twt_azi_mic\tSNP_distance\n")
        for strain in distances:
            distance_dict = distances[strain]
            found = False
            while found == False:
                closest_strain = min(distance_dict, key=distance_dict.get)
                try:
                    if closest_strain not in mtrD_strains.keys() and azi_mic[closest_strain] not in ["NA", ""]:
                        found = True
                    else:
                        del distance_dict[closest_strain]
                except KeyError:
                    del distance_dict[closest_strain]
            outfile.write(f"{strain}\t{mtrD_strains[strain]}\t{azi_mic[strain]}\t{closest_strain}\t{azi_mic[closest_strain]}\t{distance_dict[closest_strain]}\n")

def main():
    args = get_args()
    distances = read_distances(args.matrix)
    mtrD_strains = get_mtrD_strains(args.strains)
    azi_mic = get_azi(args.metadata)
    find_closest(distances, mtrD_strains, azi_mic)

if __name__ == "__main__":
    main()
