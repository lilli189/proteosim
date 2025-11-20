import re

enzyme_cleavage_patterns = {
    #schneide nach K
    'LysC': r'(?<=K)',
    #schneide vor K
    'LysN': r'(?=K)',
    'ArgC': r'(?<=R)',
    #eckige klammer bedeutet entweder oder, schneide nach K oder R, aber nicht wenn danach P kommt
    'Trypsin': r'(?<=[KR])(?!P)',
}

def digest_protein_sequence(protein_seq, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    Add a short description here.

    Parameters
    ----------
    protein_seq : str
        Protein sequence to be digested.
    cleave_pattern : str
        The expression pattern defining where to cleave the protein sequence.
    min_pep_len : int
        default = 5
        Minimum length of peptides to keep.
    max_pep_len : int
        default = 30
        Maximum length of peptides to keep.
    Returns
    -------
    """
    peptides = re.split(cleave_pattern, protein_seq)
    peptides_filtered = [x for x in peptides if min_pep_len <= len(x) <= max_pep_len]
    return peptides_filtered

def digest_protein_collection(protein_map, cleave_pattern, min_pep_len = 5, max_pep_len = 30):
    """
    Function to digest more than one protein sequence based on a cleavage patterns.

    Parameters
    ----------
    protein_map : dict
        Dictionary with protein IDs as keys and sequences as values.
    cleave_pattern : str
        The expression pattern defining where to cleave the protein sequence.
    min_pep_len : int
        Minimum length of peptides to keep.
    max_pep_len : int
        Maximum length of peptides to keep.
    Returns
    -------
    all_peptides : dict
        Dictionary with protein IDs as keys and lists of resulting peptide sequences as values.
    """
    #empty dictionary to store all peptides
    all_peptides = {}
    #iterate over all key-value pairs in the protein_map dictionary
    for protein_id, protein_seq in protein_map.items():
        peptides = re.split(cleave_pattern, protein_seq)
        #filters peptides based on length
        peptides_filtered = [x for x in peptides if min_pep_len <= len(x) <= max_pep_len]
        #stores the filtered peptides in the all_peptides dictionary with protein_id as key
        all_peptides[protein_id] = peptides_filtered

    return all_peptides 

def compute_sequence_coverage(protein_seq, peptides):
    """
    Function that reports the sequence coverage of a protein given a list of peptides.

    Parameters
    -------
    protein_seq : str
        The sequnce of the full protein.
    peptides : list of str
        A list of the digested peptides.
    Returns
    -------
    coverage : float
        The protein coverage (between 0 and 1).
    """
    #set to store covered positions, sets automatically handle duplicates
    covered_positions = set()
    #iterate over each peptide in the list of peptides
    for peptide in peptides:
        #sets the starting position for searching
        start = 0
        #loop to find all positions of the peptide in the protein sequence
        while True:
            #saves the start index of the peptide in the given protein sequence
            start = protein_seq.find(peptide, start)
            #breaks the loop if the peptide is not found
            if start == -1:
                break
            #adds all positions covered by the peptide to the covered_positions set
            for pos in range(start, start + len(peptide)):
                covered_positions.add(pos)
            #moves the start index forward to continue searching for the next occurrence
            start = start + 1
    #calculates the coverage as the ratio of covered positions to the total length of the protein sequence
    coverage = len(covered_positions) / len(protein_seq) * 100
    return coverage

