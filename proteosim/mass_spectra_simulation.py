import numpy as np
import matplotlib.pyplot as plt

amino_acid_mass_dalton = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
}

def calculate_mol_mass(peptide_seq, amino_acid_mass_dict):
    """
    A function that calculates the molecular mass of a peptide sequence.

    Parameters
    ----------
    peptide_seq : str
        The peptide sequence for which to calculate the molecular mass.
    amino_acid_mass_dict : dict
        A dictionary mapping amino acids to their molecular 

    Returns
    -------
    protein_mass : dict
        The molecular mass of the peptide sequence.

    """
    protein_mass = {}
    total_mass = 0
    for latter in peptide_seq:
        total_mass += amino_acid_mass_dict[latter]
    protein_mass[peptide_seq] = total_mass
    
    return protein_mass

def calculate_mol_mass_collection(peptides, amino_acid_mass_dict):
    """
    A function that calculates the molecular mass for a collection of peptide sequences.
    
    Parameters
    ----------
    peptides : list
        A list of peptide sequences.
    amino_acid_mass_dict : dict
        A dictionary mapping amino acids to their molecular.

    Returns
    -------
    total_mass : dict 
        The molecular mass of the peptide sequences.

    """
    mass = {}
    for pep in peptides:
        total_mass = 0
        for latter in pep:
            if latter not in amino_acid_mass_dict:
                raise ValueError(f"Unknown amino acid '{latter}' in peptide '{pep}'.")
            total_mass += amino_acid_mass_dict[latter]
        mass[pep] = total_mass
    
    return mass

def calculate_mz_collection(peptide_mass_map, charge=2, proton_mass=1.007):
    """
    A function that calculates the m/z values for a collection of peptide masses.

    Parameters
    ----------
    peptide_mass_map : dict
        A dictionary mapping peptide sequences to their molecular masses.
    charge : int
        default = 2
        The charge state of the peptide ions.
    proton_mass : float
        default = 1.007
        The mass of a proton.

    Returns
    -------
    mz_map : dict
        A dictionary mapping peptide sequences to their m/z values.

    """
    mz_map = {}
    for peptide, x in peptide_mass_map.items():
        #raise ValueError if charge <= 0:
        if charge <= 0:
            raise ValueError("Charge must be a positive integer.")
        calculation = (x + (charge * proton_mass)) / charge
        mz_map[peptide] = calculation
    return mz_map

def plot_spectrum(mz_values, random_count_range=(0, 30000), seed=42, width=2):
    """
    A function that plots a simulated mass spectrum given m/z values.

    Parameters
    ----------
    mz_values : list
        A list of m/z values.
    random_count_range : tuple
        default = (0, 30000)
        The range for generating random intensity counts.
    seed : int
        default = 42
        The random seed for reproducibility.
    width : int
        default = 2
        The width of the bars in the bar chart.
    
    Returns
    -------
    A bar chart representing the mass spectrum.

    """
    # Set random seed for reproducibility
    np.random.seed(seed)
    #define random counts for each mz value
    counts = np.random.randint(random_count_range[0], random_count_range[1], size=len(mz_values))
    width = width

    plt.bar(mz_values, counts, width=2, color='pink')
    plt.xlabel('m/z')
    plt.ylabel('Intensity')
    plt.title('Simulated Mass Spectrum')
    plt.show()

def fragment_peptide(peptide):
    """
    A function that fragments a peptide into b- and y-ions.

    Parameters
    ----------
    peptide : str
        The peptide sequence to be fragmented.

    Returns
    -------
    fragments : list
        A list of fragment ions (b- and y-ions comnined).
    """
    fragments = []
    length = len(peptide)
    
    # Generate b-ions
    for i in range(1, length):
        # Generate b-ion by taking the 1., 2., ..., (n-1). amino acids and appending to fragments
        b_ion = peptide[:i]
        fragments.append(b_ion)
    
    # Generate y-ions
    for i in range(1, length):
        # Generate y-ion by taking backwards from the end
        y_ion = peptide[i:]
        fragments.append(y_ion)
    
    return fragments