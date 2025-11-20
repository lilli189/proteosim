from proteosim.file_handling import read_fasta

def test_read_faster():
    tmp_fasta_path = 'data/dummy_proteins.fasta'
    protein_map = read_fasta(tmp_fasta_path)
    
    # Replace the strings with your fasta content
    # which you expect to be now available as a dictionary
    assert protein_map["A069384"] == "BWETRSDFVRTBTWWTQDBTRTB"
    assert protein_map["PG356G"] == "EGTZNGHNTEZNETZRZNZNZNT"