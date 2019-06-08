#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import pdb_utils
import utils
import extract

from operator import attrgetter

RESIDUES_LIST = "data/residues.list"
ATOMS_LIST = "data/atoms.list"


def clean_format(f):
    """CleanFormat is a function used to format different platform formats to unix. Users need to install dos2unix

    .. warning:: dos2unix required

    """
    os.system( "mac2unix -q %s" %f )
    os.system( "dos2unix -q %s" %f )


def normalize_structure(struct, out_file = None, index_file=None, extract_file = None):
    """
    Different prediction methods may give different PDB formats. The standard PDB format considered is the [1992]
    (https://www.rcsb.org/pdb/file_formats/pdb/pdbguide2.2/PDB_format_1992.pdf) format.


    To normalize the structure format, we need to:
    * Correct the residue names and the atom names:
    The name mapping dictionaries are: defined as "RESIDUES_LIST" and "ATOMS_LIST" in the data folder.

    Index file format:

    * select the part of structure we are interested. The part of structure is assigned in "index_file". The selection grammar is: [Chain id]:[Beg of Residue number]:[length of fragment]
        + [Chain id] is single letter character
        + [Beg of Residue number] is the number of the beginning residue in the PDB file. This number is exact the "res. seq no." field of the PDB file without any renumbering.
        + [length of fragment] is the length of the fragment. numbering of the residues in the fragment is not considered.
        + for example, A:1:31,A:33:29 includes two fragments of 31 residue and 29 residues. The first fragment starts from residue 1 and until residue 31. leaving out residue 32, the second fragment starts from residue 33 and ends till 29 residues after. The number of these 29 residues are not considered except the first one.

        Args:
        path (str): The path of the file to wrap
        findex_file (fn) The :class:Y instance to wrap
        extract_file (???): Extract the interesting part. To do this, we need to set the output file for "extract_file" parameter

        Returns:
        None

    """
    pdb_normalizer = pdb_utils.PDBNormalizer( RESIDUES_LIST, ATOMS_LIST )
    ok = pdb_normalizer.parse( struct, out_file )
    if not ok:
        sys.stderr.write("ERROR: structure not normalized!\n")
    else:
        sys.stderr.write("INFO: Normalization succeded!\n")
    if not extract_file is None:
        coords=open(index_file).read()
        extract.extract_PDB(SOLUTION_NORMAL, coords, extract_file)
        sys.stderr.write("INFO:	structure extracted\n")


def interaction_network_fidelity(native_file, native_index, prediction_file, prediction_index, check_seq=False):
    """
    Index file format:

    * select the part of structure we are interested. The part of structure is assigned in "index_file". The selection grammar is: [Chain id]:[Beg of Residue number]:[length of fragment]
        + [Chain id] is single letter character
        + [Beg of Residue number] is the number of the beginning residue in the PDB file. This number is exact the "res. seq no." field of the PDB file without any renumbering.
        + [length of fragment] is the length of the fragment. numbering of the residues in the fragment is not considered.
        + for example, A:1:31,A:33:29 includes two fragments of 31 residue and 29 residues. The first fragment starts from residue 1 and until residue 31. leaving out residue 32, the second fragment starts from residue 33 and ends till 29 residues after. The number of these 29 residues are not considered except the first one.
    """
    res_struct = pdb_utils.PDBStruct()
    res_struct.load(native_file, native_index)
    res_raw_seq = res_struct.raw_sequence()

    sol_struct = pdb_utils.PDBStruct()
    sol_struct.load( prediction_file, prediction_index )
    sol_raw_seq = sol_struct.raw_sequence()

    if sol_raw_seq != res_raw_seq and check_seq:
        sys.stderr.write("ERROR Result sequence != Solution sequence!\n")
        sys.stderr.write("DATA Solution sequence --> '%s'\n" %sol_raw_seq )
        sys.stderr.write("DATA Result sequence	 --> '%s'\n" %res_raw_seq )
        return(-1)

    # computes the RMSD
    comparer = pdb_utils.PDBComparer()
    rmsd = comparer.rmsd(sol_struct, res_struct)
    INF_ALL = comparer.INF(sol_struct, res_struct, type="ALL")
    DI_ALL = rmsd / INF_ALL
    INF_WC = comparer.INF(sol_struct, res_struct, type="PAIR_2D")
    INF_NWC = comparer.INF(sol_struct, res_struct, type="PAIR_3D")
    INF_STACK = comparer.INF(sol_struct, res_struct, type="STACK")
    return (rmsd, DI_ALL, INF_ALL, INF_WC, INF_NWC, INF_STACK)


def calc_RMSD(native_file, native_index, prediction_file, prediction_index, PVALUE = "-"):
    """
    PVALUE set according to Hajdin et al., RNA (7) 16, 2010, either "+" or "-"
    """
    res_struct = pdb_utils.PDBStruct()
    res_struct.load( native_file, native_index )
    res_raw_seq = res_struct.raw_sequence()

    sol_struct = pdb_utils.PDBStruct()
    sol_struct.load( prediction_file, prediction_index )
    sol_raw_seq = sol_struct.raw_sequence()

    if( sol_raw_seq != res_raw_seq ):
        sys.stderr.write("ERROR Result sequence != Solution sequence!\n")
        sys.stderr.write("DATA Solution sequence --> '%s'\n" %sol_raw_seq )
        sys.stderr.write("DATA Result sequence   --> '%s'\n" %res_raw_seq )
        return(-1)
    # computes the RMSD
    comparer = pdb_utils.PDBComparer()
    rmsd = comparer.rmsd( sol_struct, res_struct )
    sys.stderr.write("INFO Partial RMSD --> %f\n" %rmsd )
    pvalue = comparer.pvalue( rmsd, len(sol_raw_seq), PVALUE )
    sys.stderr.write("INFO Partial P-Value --> %e\n" %pvalue )
    return(rmsd, pvalue)
