from pyteomics import achrom 
import matplotlib.pyplot as plt

def predict_lc_retention_times(peptides):
    """
    A function that maps a list of peptide sequences to their predicted LC retention times.

    Parameters
    ----------
    peptides : list of str
        A list of peptide sequences.

    Returns
    -------
    dict
        A dictionary mapping each peptide sequence to its predicted retention time (RT).
    """
    RT_of_peptides = {}
    for peptide in peptides:
        RT = achrom.calculate_RT(peptide, achrom.RCs_guo_ph7_0)
        RT = round(float(RT), 2) 
        RT_of_peptides[peptide] = RT

    return RT_of_peptides

def plot_retention_time(retention_times, resolution=30):
    """
    A function that plots a histogram of retention times.

    Parameters
    ----------
    retention_times : list
        A list of retention times.
    resolution : int
        default = 30
        Controls the number of bins in the histogram.

    Returns
    -------
    A histogram plot of retention times.

    """
    plt.figure()
    plt.hist(retention_times, bins=resolution, color="pink", edgecolor='black')
    plt.xlabel('Relative retention time')
    plt.ylabel('Intensities')
    plt.title('Histogram of Retention Times')
    return plt.show()

def select_retention_time_window(peptide_rt_map, lower_ret_time, upper_ret_time):
    """
    A function that selects peptides within a specified retention time window.

    Parameters
    ----------
    peptide_rt_map : dict
        A dictionary mapping peptide sequences to their retention times.
    lower_ret_time : float
        The lower bound of the retention time (in minutes) window. 
    upper_ret_time : float
        The upper bound of the retention time (in minutes) window.

    Returns
    -------
    A dictionary of peptides with retention times within the specified window.

    """
    #create an empty dictionary to store selected peptides
    selected_peptides = {}
    #iterate through the input dictionary
    for peptide, rt in peptide_rt_map.items():
        if lower_ret_time <= rt <= upper_ret_time:
            #save the peptide (key) and its retention time (value) to the new dictionary
            selected_peptides[peptide] = rt
    return selected_peptides