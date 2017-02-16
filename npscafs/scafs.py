# -*- coding: utf-8 -*-
"""
#####
Scafs
#####

*Created on Tue Feb 14 14:00 2017 by A. Pahl*

Handling Natural Product Scaffolds.
"""

import os
from copy import deepcopy

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
import rdkit.Chem.Descriptors as Desc
import rdkit.Chem.Scaffolds.MurckoScaffold as MurckoScaffold
Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
Draw.DrawingOptions.atomLabelFontSize = 18

from rdkit_ipynb_tools import tools


def sort_by_size(s):
    """Sorts the input elements by decreasing length.
    Returns a list with the sorted elements."""
    return sorted(s, key=len, reverse=True)


def calc_murcko_scaf(mol):
    "Calculate the Murcko scaffold from a molecule and return as Smiles."
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol)


def get_scaffolds_for_mol(mol):
    """Generates the BRICS scaffold from the input mol.
    Returns a set with the scaffolds as Smiles."""
    scaf_set = set()
    frags = Chem.FragmentOnBRICSBonds(mol)
    frag_list = Chem.MolToSmiles(frags).split(".")
    frag_list = [f.replace("*", "H") for f in frag_list]
    for frag in frag_list:
        mol = Chem.MolFromSmiles(frag)
        if Desc.RingCount(mol) > 0:
            murcko = calc_murcko_scaf(mol)
            scaf_set.add(murcko)
    return scaf_set


def get_scaffolds_for_mols(mol_list):
    """Generates the BRICS scaffold from the input mols.
    Returns a set with the scaffolds as Smiles."""
    scaf_set = set()
    for mol in mol_list:
        scaf_set.update(get_scaffolds_for_mol(mol))
    return scaf_set


def pipe_get_scaffolds(stream, scaf_set, summary=None, comp_id="pipe_get_scaffold"):
    "Put the scaffolds of the molecules in the stream into scaf_set."
    rec_counter = 0
    for rec in stream:
        mol = rec["mol"]
        if mol:
            scaf_set.update(get_scaffolds_for_mol(mol))
            rec_counter += 1
            if summary is not None:
                summary[comp_id] = rec_counter
            yield rec


def pipe_count_scaffolds(stream, scaf_ctr, summary=None, comp_id="pipe_count_scaffolds"):
    """Count the scaffolds of the molecules in the stream.
    The result is put in `scaf_ctr` (of type collections.Counter)."""
    rec_counter = 0
    for rec in stream:
        mol = rec["mol"]
        if mol:
            scaf_set = get_scaffolds_for_mol(mol)
            for scaf in scaf_set:
                scaf_ctr[scaf] += 1
            rec_counter += 1
            if summary is not None:
                summary[comp_id] = rec_counter
            yield rec


