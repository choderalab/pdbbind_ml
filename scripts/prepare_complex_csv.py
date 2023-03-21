import argparse
from glob import glob
import os
import pandas


################################################################################
def get_args():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-i", "--in_files", required=True, help="Input glob/directory.")
    parser.add_argument("-o", "--out_file", required=True, help="Output CSV file.")

    return parser.parse_args()


def main():
    args = get_args()

    if os.path.isdir(args.in_files):
        file_glob = os.path.join(args.in_files, "*", "*.sdf")
    else:
        file_glob = args.in_files
    all_sdf_files = glob(file_glob)
    assert len(all_sdf_files) > 0, "No SDF files found."
    print(f"Found {len(all_sdf_files)} SDF files", flush=True)

    csv_cols = ["ligand_id", "du_structure", "docked_file"]
    df = [("", fn.split("/")[-2], fn) for fn in all_sdf_files]
    df = pandas.DataFrame(df, columns=csv_cols)
    df.to_csv(args.out_file)


if __name__ == "__main__":
    main()
