"""Choose the ligand residue: explicit mask, residue name, or auto-detect."""

from . import UnleapError

_AMINO = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # protonation / disulfide / alternate names
    "ASH", "GLH", "HID", "HIE", "HIP", "HYP", "LYN", "CYX", "CYM", "CYN", "ARN", "TYM", "AR0",
    "MSE", "SEC", "PYL",
    # caps
    "ACE", "NME", "NHE", "NH2",
}
_NUCLEIC = {
    "DA", "DC", "DG", "DT", "A", "C", "G", "U", "RA", "RC", "RG", "RU",
    "DA5", "DC5", "DG5", "DT5", "DA3", "DC3", "DG3", "DT3", "DAN", "DCN", "DGN", "DTN",
    "A5", "C5", "G5", "U5", "A3", "C3", "G3", "U3", "AN", "CN", "GN", "UN",
    "RA5", "RC5", "RG5", "RU5", "RA3", "RC3", "RG3", "RU3",
}
_SOLVENT = {"WAT", "HOH", "TIP3", "TIP4", "TP3", "TP4", "T3P", "T4P", "SPC", "SPCE", "OPC", "OPC3",
            "SOL", "H2O"}
_IONS = {
    "NA", "Na+", "NA+", "K", "K+", "CL", "Cl-", "CL-", "MG", "MG2", "CA", "CA2", "ZN", "ZN2",
    "Li+", "LI", "Rb+", "RB", "Cs+", "CS", "F-", "Br-", "BR", "I-", "IOD", "SR", "BA", "CU",
    "FE", "FE2", "MN", "CO", "NI", "CD", "HG",
}
STANDARD = _AMINO | _NUCLEIC | _SOLVENT | _IONS


def _is_standard(name):
    if name in STANDARD:
        return True
    # N-/C-terminal amino-acid variants (NALA, CGLY, ...)
    return len(name) == 4 and name[0] in "NC" and name[1:] in _AMINO


def ligand_mask(struct, mask=None, resname=None):
    """Return an Amber mask string selecting the ligand."""
    if mask and resname:
        raise UnleapError("give either --mask or --resname, not both")
    if mask:
        return mask
    if resname:
        if not any(r.name == resname for r in struct.residues):
            raise UnleapError(f"no residue named {resname!r} in the topology")
        return f":{resname}"

    counts = {}
    for res in struct.residues:
        if not _is_standard(res.name):
            counts[res.name] = counts.get(res.name, 0) + 1
    if not counts:
        raise UnleapError("auto-detect found no non-standard residue; use --mask or --resname")
    if len(counts) > 1 or next(iter(counts.values())) > 1:
        found = ", ".join(f"{k} (x{v})" for k, v in sorted(counts.items()))
        raise UnleapError(f"auto-detect is ambiguous; candidates: {found}. Use --mask or --resname")
    return f":{next(iter(counts))}"
