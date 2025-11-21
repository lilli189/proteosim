# Proteosim — Bottom-Up Proteomics Simulation Package
proteosim is a Python package that simulates a complete bottom-up proteomics pipeline.
It reads protein FASTA files, performs in-silico digestion, predicts LC retention times, and simulates MS1/MS2 spectra including peptide masses, m/z values, and fragment ions.
This repository also contains an end-to-end tutorial notebook.

## How to install the package
Clone the repository and install the environment
<pre> ```bash git clone &lt;your-repo-url&gt; cd proteosim ``` </pre>
Install dependencies:
<pre> ```bash pip install -r requirements.txt ``` </pre>
You can now import it in Python:
<pre> ```bash import proteosim as ps ``` </pre>

## What the Package Does
proteosim provides a full workflow for simulating bottom-up proteomics experiments.
It takes protein sequences as input and performs a multi-step transformation:

### Input
- Multi-protein FASTA file
  → Parsed into {protein_id: sequence}
### Processing Pipeline
1. Protein digestion
   - Trypsin and other protease rules included
   - Adjustable minimum/maximum peptide length
   - Sequence coverage computation
2. Chromatography simulation
   - Predict liquid chromatography (LC) retention times
   - Visualize chromatograms
   - Select time windows of interest
3. Mass spectrometry simulation
   - Compute peptide masses
   - Compute m/z values for MS1
   - Generate MS1 spectra
   - Fragment peptides (b/y ions)
   - Simulate MS2 spectra
### Output
- Dictionaries of digested peptides
- Retention-time tables
- Mass/mz mappings
- Simulated MS1 and MS2 spectra (via plotting functions)

## Modules & Functions Overview

This section provides newcomers with a quick orientation of the `proteosim` package structure and the key functions available in each module.

### FASTA Handling

- **`read_fasta(path)`**  
  Reads a FASTA file and returns a dictionary mapping protein IDs to amino-acid sequences.


### Digestion Module

- **`digest_protein_sequence(seq, cleave_pattern, min_pep_len, max_pep_len)`**  
  Digests a single protein sequence based on a cleavage pattern.
- **`digest_protein_collection(protein_map, cleave_pattern, min_pep_len, max_pep_len)`**  
  Digests all proteins in a FASTA dataset.
- **`compute_sequence_coverage(seq, peptide_list)`**  
  Calculates the percentage of the protein covered by resulting peptides.
- **`enzyme_cleavage_patterns`**  
  Provides cleavage patterns for Trypsin, LysC, LysN, and ArgC.

### Chromatography Module

- **`predict_lc_retention_times(peptide_list)`**  
  Predicts approximate LC retention times for each peptide.
- **`plot_retention_time(retention_times, resolution)`**  
  Visualizes a simplified chromatographic profile.
- **`select_retention_time_window(rt_map, lower, upper)`**  
  Filters peptides by a retention-time interval.

---

### Mass Spectrometry Module

- **`calculate_mol_mass(peptide_seq, aa_mass_dict)`**  
  Computes the monoisotopic mass of a peptide.
- **`calculate_mol_mass_collection(peptides, aa_mass_dict)`**  
  Computes masses for a full peptide list.
- **`calculate_mz_collection(mass_map, charge, proton_mass)`**  
  Converts peptide masses into m/z values for MS1.
- **`fragment_peptide(peptide_seq)`**  
  Generates fragment ions (b/y series) for MS2 simulation.
- **`plot_spectrum(mz_values, seed, width)`**  
  Plots an MS1 or MS2 mass spectrum.
- **`amino_acid_mass_dalton`**  
  Dictionary of amino-acid masses used for simulations.


## Full Example Workflow (Tutorial Notebook)

A complete end-to-end example is provided in the notebook:
tutorials/ms_experiment_final.ipynb
