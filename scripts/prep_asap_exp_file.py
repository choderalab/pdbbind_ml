"""
Script to prepare the output of parse_pdbbind_affinities.py for ingestion by the
asapdiscovery pipeline (namely, asapdiscovery.data.utils.cdd_to_schema).
"""

import argparse
import numpy as np
import pandas
from rdkit import Chem


def convert_row(r):
    """
    Convert a row of the input CSV file into a row of output format.

    Parameters
    ----------
    r : pandas.Series
        Row of input DataFrame

    Returns
    -------
    pandas.Series
        Row of output DataFrame
    """

    # Needed output columns:
    #  * suspected_SMILES
    #  * Canonical PostEra ID
    #  * ProteaseAssay_Fluorescence_Dose-Response_Weizmann: Avg pIC50
    #  * ProteaseAssay_Fluorescence_Dose-Response_Weizmann: IC50 CI (Lower) (µM)
    #  * ProteaseAssay_Fluorescence_Dose-Response_Weizmann: IC50 CI (Upper) (µM)

    if ("iso_smiles" in r.index) and (not pandas.isna(r["iso_smiles"])):
        smi = r["iso_smiles"]
    else:
        smi = r["str_smiles"]
    ligand_id = r["ligand_id"]
    # For now, only handle IC50 measurements
    # Ki/Kd conversion will need to use Cheng-Prussof
    pic50 = -np.log10(r["measurement_value"])
    # pic50 = convert_to_pic50(r["measurement_value"], r["measurement_type"])

    return pandas.Series(
        [smi, ligand_id, pic50, np.nan, np.nan],
        index=[
            "suspected_SMILES",
            "Canonical PostEra ID",
            "ProteaseAssay_Fluorescence_Dose-Response_Weizmann: Avg pIC50",
            "ProteaseAssay_Fluorescence_Dose-Response_Weizmann: IC50 CI (Lower) (µM)",
            "ProteaseAssay_Fluorescence_Dose-Response_Weizmann: IC50 CI (Upper) (µM)",
        ],
    )


def convert_to_pic50(measurement_value, measurement_type):
    """
    Convert an experimental measurement to a pIC50 value.

    Parameters
    ----------
    measurement_value : float
        Measurement value
    measurement_type : str
        What type of measurement the input is. Should be one of ["Ki", "IC50", "Kd"]

    Returns
    -------
    float
        Measurement converted to pIC50 value
    """
    pass


################################################################################
def get_args():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-i", "--in_file", required=True, help="Input CSV file.")
    parser.add_argument("-o", "--out_file", required=True, help="Input CSV file.")

    return parser.parse_args()


def main():
    args = get_args()

    df_in = pandas.read_csv(args.in_file, index_col=0)
    df_out = df_in.apply(convert_row, axis=1)

    # Filter out invalid SMILES from RCSB
    keep_idx = []
    for _, r in df_out.iterrows():
        smi = r["suspected_SMILES"]
        keep_idx += [bool(Chem.MolFromSmiles(smi))]
        if not keep_idx[-1]:
            print(r["Canonical PostEra ID"], smi, flush=True)
    print(
        f"Removing {len(keep_idx) - sum(keep_idx)} entries with invalid SMILES",
        flush=True,
    )
    df_out = df_out.loc[keep_idx, :]

    df_out.to_csv(args.out_file, index=False)


if __name__ == "__main__":
    main()
