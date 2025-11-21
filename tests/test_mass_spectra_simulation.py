from proteosim.mass_spectra_simulation import (calculate_mol_mass, calculate_mol_mass_collection, calculate_mz_collection, fragment_peptide)

def test_calculate_mol_mass():
    aa_mass_dict = {
        'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
        'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
        'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
        'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
        'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
    }
    test_peptide = "TESTTHISPEPTIDE"
    expected_mass = {'TESTTHISPEPTIDE': 1638.75}
    calculated_mass = calculate_mol_mass(test_peptide, aa_mass_dict)
    assert calculated_mass == expected_mass

test_calculate_mol_mass()

def test_calculate_mol_mass_collection():
    aa_mass_dict = {
        'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
        'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
        'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
        'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
        'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
    }
    peptides = ["PEPTIDE", "TEST"]
    expected = {'PEPTIDE': 781.84, 'TEST': 418.42}
    actual = calculate_mol_mass_collection(peptides,aa_mass_dict)

    assert actual == expected

    # Test for unknown amino acid, test if ValueError is raised
    try:
        calculate_mol_mass_collection(["MÄSS"], aa_mass_dict)
        assert False, "Expected a ValueError to be raised."
    except ValueError as e:
        assert str(e) == "Unknown amino acid 'Ä' in peptide 'MÄSS'."

test_calculate_mol_mass_collection()

def test_calculate_mz_collection():
    peptide_mass_map = {'MATSR': 546.65, 'TYLDK': 620.71, 'DVSLPR': 667.77}
    actual = calculate_mz_collection(peptide_mass_map, charge=2)
    expected = {'MATSR': 274.332, 'TYLDK': 311.362, 'DVSLPR': 334.892}

    assert actual == expected

    # Test for charge=0, test if ValueError is raised
    try:
        calculate_mz_collection({"MÄSS":500}, charge=0)
        assert False, "Expected a positive integer as charge."
    except ValueError as e:
        assert str(e) == "Charge must be a positive integer."

test_calculate_mz_collection()

def test_fragment_peptide():
    peptide = "TESTPEPTIDE"
    expected = ['T', 'TE', 'TES', 'EST', 'ST', 'T']
    actual = fragment_peptide("TEST")

    assert set(actual) == set(expected)

test_fragment_peptide()