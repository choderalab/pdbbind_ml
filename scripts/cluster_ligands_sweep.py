"""
Sweep over a range of clustering thresholds to see how many clusters you get.
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas
import re
import seaborn as sns

from cluster_ligands import cluster_mols, load_fp_db


def cluster_sweep(sim_arr, cutoffs, ligand_idx_dict, cluster_cache=None):
    all_n_clusters = []
    print("Clustering", flush=True)
    for t in cutoffs:
        print(f"\tCutoff: {t:0.2f}", flush=True)
        lig_clusters = cluster_mols(sim_arr, t, ligand_idx_dict)
        all_n_clusters += [len(np.unique(list(lig_clusters.values())))]
        if cluster_cache:
            # Replace decimal with underscore to make the filename better
            t = re.sub("\.", "_", f"{t:0.2f}")
            pandas.DataFrame(
                {"ligand_id": lig_clusters.keys(), "cluster": lig_clusters.values()}
            ).to_csv(os.path.join(cluster_cache, f"{t}.csv"))

    return all_n_clusters


def load_cluster_sweep(cluster_cache, cutoffs):
    """
    Load clusterings from directory.

    Parameters
    ----------
    cluster_cache : str
        Directory where clustering CSV files are stored
    cutoffs : List[float]
        List of cutoffs to check for

    Returns
    -------
    List[int]
        List of numbers of clusters
    """
    all_n_clusters = []
    for t in cutoffs:
        # Replace decimal with underscore to make the filename better
        t_fn = re.sub("\.", "_", f"{t:0.2f}")
        fn = os.path.join(cluster_cache, f"{t}.csv")
        if os.path.isfile(fn):
            df = pandas.read_csv(fn, index_col=0)
            all_n_clusters += [df.shape[0]]
        else:
            print(f"No file found for cutoff {t:0.2f}", flush=True)
            # Already know we're going to have to redo the clustering, so no need to
            #  keep checking files
            break

    return all_n_clusters


################################################################################
def get_args():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-i", "--in_file", required=True, help="Input CSV file.")
    parser.add_argument("-o", "--out_file", required=True, help="Output plot file.")
    parser.add_argument(
        "-c", "--cache_file", help="Optional file to save the full distance matrix to."
    )
    parser.add_argument(
        "-cl_o", "--cluster_out", help="Output directory to store clusterings."
    )

    parser.add_argument(
        "-ts",
        "--thresh_start",
        type=float,
        default=0.2,
        help="Starting threshold for the sweep.",
    )
    parser.add_argument(
        "-te",
        "--thresh_end",
        type=float,
        default=1,
        help="Ending threshold for the sweep.",
    )
    parser.add_argument(
        "-s", "--step", type=float, default=0.05, help="Step size for the sweep."
    )

    return parser.parse_args()


def main():
    args = get_args()

    # Load input
    fpdb, ligand_idx_dict = load_fp_db(args.in_file)
    print(fpdb.NumFingerPrints(), len(ligand_idx_dict), flush=True)

    # Load/calculate similarities
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

    # add step to the end to make it inclusive
    all_thresh = np.arange(args.thresh_start, args.thresh_end + args.step, args.step)
    print(all_thresh, flush=True)

    if args.cluster_out:
        if os.path.isdir(args.cluster_out):
            # Just load clusterings
            all_n_clusters = load_cluster_sweep(args.cluster_out, all_thresh)
            if len(all_n_clusters) != len(all_thresh):
                print(
                    (
                        f"Found {len(all_n_clusters)} clusterings but "
                        f"{len(all_thresh)} cutoffs, will recalculate clusters."
                    ),
                    flush=True,
                )
        else:
            # Create directory for saving clusterings
            os.makedirs(args.cluster_out)
            all_n_clusters = []

    # Clustering sweep
    if (not args.cluster_out) or (len(all_n_clusters) != len(all_thresh)):
        all_n_clusters = cluster_sweep(
            all_sims, all_thresh, ligand_idx_dict, args.cluster_out
        )

    # Plot
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.scatterplot(x=all_thresh, y=all_n_clusters, ax=ax)
    ax.set_yscale("log")
    ax.set_xlabel("Clustering Cutoff")
    ax.set_ylabel("Number of Clusters")
    fig.savefig(args.out_file, bbox_inches="tight", dpi=200)


if __name__ == "__main__":
    main()
