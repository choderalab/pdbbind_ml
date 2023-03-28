"""
Cluster a group of ligands. Input CSV file must at minimum have the columns "ligand_id"
and either "str_smiles" or "iso_smiles". "str_smiles" will be preferred over
"iso_smiles" if both are present. Output CSV file contains the columns "ligand_id" and
"cluster". "ligand_id" corresponds to the same ids in the input CSV file, and "cluster"
gives the index of the cluster that the given ligand was placed in.
"""

import argparse
import multiprocessing as mp
import numpy as np
from openeye import oechem, oegraphsim
import os
import pandas
from rdkit.ML.Cluster import Butina


def build_fp_db(smiles_list, ligand_ids):
    """
    Build an OE fingerprint database containing all the passed SMILES.

    Parameters
    ----------
    smiles_list : List[str]
        List of SMILES strings
    ligand_ids : List[str]
        List of ligand ids corresponding to SMILES

    Returns
    -------
    oegraphsim.OEFPDatabase
        OE fingerprint database containing all the given SMILES
    Dict[str, int]
        Mapping from ligand id to index in fingerprint database
    """
    fpdb = oegraphsim.OEFPDatabase(oegraphsim.OEFPType_MACCS166)

    ligand_idx_dict = {}
    for smi, lig in zip(smiles_list, ligand_ids):
        # Ligand ids not necessarily unique, so make sure we haven't already seen it
        if lig in ligand_idx_dict:
            continue

        # Create the molecule object from SMILES
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)

        # Add to db
        lig_idx = fpdb.AddFP(mol)
        if lig_idx == -1:
            print(f"Error adding ligand {lig} to database.", flush=True)
        ligand_idx_dict[lig] = lig_idx

    return fpdb, ligand_idx_dict


def mp_func(fp, fpdb):
    return [s.GetScore() for s in fpdb.GetScores(fp)]


################################################################################
def get_args():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-i", "--in_file", required=True, help="Input CSV file.")
    parser.add_argument("-o", "--out_file", required=True, help="Output CSV file.")

    parser.add_argument(
        "-c", "--cache_file", help="Optional file to save the full distance matrix to."
    )

    parser.add_argument(
        "-n",
        "--num_procs",
        type=int,
        default=1,
        help="Number of concurrent processes to run.",
    )

    parser.add_argument(
        "-t",
        "--thresh",
        type=float,
        default=0.5,
        help="Cutoff threshold for clustering.",
    )

    return parser.parse_args()


def main():
    args = get_args()

    # Load input
    df = pandas.read_csv(args.in_file, index_col=0)
    if "str_smiles" in df.columns:
        smi_col = "str_smiles"
    elif "iso_smiles" in df.columns:
        smi_col = "iso_smiles"
    else:
        raise RuntimeError("Couldn't find valid SMILES column.")

    fpdb, ligand_idx_dict = build_fp_db(df[smi_col], df["ligand_id"])
    print(fpdb.NumFingerPrints(), len(ligand_idx_dict), flush=True)

    if args.cache_file and os.path.exists(args.cache_file):
        print("Loading similarities from cache file.", flush=True)
        all_sims = np.loadtxt(args.cache_file)
    else:
        print("Calculating similarities.", flush=True)
        # mp_args = [(fp, fpdb) for fp in fpdb.GetFingerPrints()]
        # n_procs = min(args.num_procs, len(mp_args), mp.cpu_count())
        # with mp.Pool(processes=n_procs) as pool:
        # all_sims = np.asarray(pool.starmap(mp_func, mp_args))
        # can redo this at some point to not calculate unneccessary values
        all_sims = np.asarray(
            [
                [s.GetScore() for s in fpdb.GetScores(fp)]
                for fp in fpdb.GetFingerPrints()
            ]
        )

        if args.cache_file:
            np.savetxt(args.cache_file, all_sims)

    print("Clustering", flush=True)
    all_dists = 1 - all_sims
    n_mols = all_dists.shape[0]
    # Get the lower triangular entries, not including the diagonal (k=-1)
    all_dists = all_dists[np.tril_indices_from(all_dists, k=-1)]
    all_clusters = Butina.ClusterData(all_dists, n_mols, args.thresh, isDistData=True)

    lig_clusters = {}
    idx_lig_dict = {idx: lig for lig, idx in ligand_idx_dict.items() if idx != -1}
    for i, cl in enumerate(all_clusters):
        for lig_idx in cl:
            lig_clusters[idx_lig_dict[lig_idx]] = i

    out_df = pandas.DataFrame(
        {"ligand_id": lig_clusters.keys(), "cluster": lig_clusters.values()}
    )
    out_df.to_csv(args.out_file)


if __name__ == "__main__":
    main()
