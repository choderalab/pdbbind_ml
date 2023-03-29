import argparse
import matplotlib.pyplot as plt
import pandas
import seaborn as sns


################################################################################
def get_args():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-i", "--in_file", required=True, help="Input CSV file.")
    parser.add_argument(
        "-lo",
        "--ligands_out",
        required=True,
        help="Output filename for counts of each ligand.",
    )
    parser.add_argument(
        "-co",
        "--counts_out",
        required=True,
        help="Output base filename for counts of counts.",
    )

    return parser.parse_args()


def main():
    args = get_args()

    # Load input
    df = pandas.read_csv(args.in_file, index_col=0)

    # Count occurrences of each ligand id
    lig_counts = df["ligand_id"].value_counts().to_frame().reset_index()
    lig_counts = lig_counts.rename(columns={"index": "ligand_id", "ligand_id": "count"})
    lig_counts.to_csv(args.ligands_out, index=False)

    # Count occurrences of each count
    count_counts = lig_counts["count"].value_counts().to_frame().reset_index()
    count_counts = count_counts.rename(
        columns={"index": "occurrences", "count": "n_ligs"}
    )

    # # Get order for x axis
    # x_order = map(str, sorted(lig_counts["count"].values))

    # # Set to str so bins aren't numerical
    # lig_counts.loc[:, "count"] = lig_counts["count"].astype(str)

    # Plot
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.barplot(count_counts, x="occurrences", y="n_ligs")
    # sns.histplot(x=map(str, sorted(lig_counts["count"].values)))
    ax.set_xlabel("Ligand Occurrences")
    ax.set_ylabel("Number of Ligands")
    ax.set_title("Ligand Distribution in PDBBind")
    fig.savefig(f"{args.counts_out}.png", bbox_inches="tight", dpi=200)

    count_counts.to_csv(f"{args.counts_out}.csv", index=False)


if __name__ == "__main__":
    main()
