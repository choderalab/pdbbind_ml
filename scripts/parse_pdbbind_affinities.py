"""
Convert PDBBind index files into a CSV file.
"""
import argparse
from openeye import oechem
import os
import pandas
import re

# Strings to avoid floating point weirdness when saving CSV files
UNITS_MULT = {
    "fM": "e-15",
    "pM": "e-12",
    "nM": "e-9",
    "uM": "e-6",
    "mM": "e-3",
}


def fix_units(measurement):
    # Function to add the correct multiplier for the unit
    val = measurement[:-2]
    units = measurement[-2:]
    return val + UNITS_MULT[units]


def get_smiles(fn):
    # Load molecule
    ifs = oechem.oemolistream()
    ifs.SetFlavor(
        oechem.OEFormat_SDF,
        oechem.OEIFlavor_SDF_Default,
    )
    ifs.open(fn)
    mol = oechem.OEGraphMol()
    oechem.OEReadMolecule(ifs, mol)
    ifs.close()

    # # Delete explicit hydrogens
    # for a in mol.GetAtoms():
    #     if a.GetAtomicNum() == 1:
    #         mol.DeleteAtom(a)

    # Get smiles
    return oechem.OEMolToSmiles(mol)


def parse_index(fn, str_dir):
    df_cols = [
        "PDB_id",
        "ligand_id",
        "resolution",
        "release_year",
        "measurement_type",
        "measurement_relationship",
        "measurement_value",
        "pdb_path",
        "sdf_path",
        "smiles",
    ]
    entries = []
    for line in open(fn).readlines():
        if line[0] == "#":
            continue

        line = line.split()
        pdb_id = line[0]
        ligand_id = line[-1].strip("()")
        entry = [pdb_id, ligand_id] + line[1:3]
        try:
            measurement = line[4]
        except IndexError as e:
            print(line, flush=True)
            raise e

        # Split up the measurement value entry
        split = re.split(r"([<>=~]+)", measurement)
        if len(split) != 3:
            print(
                f"Unable to split measurement value {measurement} for {pdb_id}",
                flush=True,
            )
            continue
        # Adjust units for measurement value
        split[2] = fix_units(split[2])
        # Add to entry
        entry.extend(split)

        # Try and find the protein PDB file and ligand SDF file
        base_fn = os.path.join(str_dir, pdb_id, f"{pdb_id}")
        prot_fn = f"{base_fn}_protein.pdb"
        if not os.path.isfile(prot_fn):
            print(f"PDB file {prot_fn} not found", flush=True)
            continue
        prot_fn = os.path.abspath(prot_fn)
        lig_fn = f"{base_fn}_ligand.sdf"
        if not os.path.isfile(lig_fn):
            print(f"SDF file {lig_fn} not found", flush=True)
            continue
        lig_fn = os.path.abspath(lig_fn)
        # Add file names to entry
        entry.extend([prot_fn, lig_fn])

        # Load SMILES
        smiles = get_smiles(lig_fn)
        entry.append(smiles)

        entries.append(entry)

    return pandas.DataFrame(entries, columns=df_cols)


################################################################################
def get_args():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument(
        "-i", "--in_files", required=True, nargs="+", help="Input index files."
    )
    parser.add_argument(
        "-d", "--str_dir", required=True, help="Directory to look for structure files."
    )
    parser.add_argument("-o", "--out_file", required=True, help="Output CSV file.")

    return parser.parse_args()


def main():
    args = get_args()

    # Parse all index files
    all_dfs = [parse_index(fn, args.str_dir) for fn in args.in_files]
    # Combine and save
    df = pandas.concat(all_dfs)
    df.to_csv(args.out_file)


if __name__ == "__main__":
    main()
