from proteosim.protein_digestion import digest_protein_collection, compute_sequence_coverage

def test_digest_protein_collection():
    dummy_proteins = {"Sequnece_ID_1": "AKLLSEBN", "Sequnece_ID_2": "HSJERLSDE"}
    cleavage_patter_trypsin = r'(?<=[KR])(?!P)'
    cleavage_patter_LysC = r'(?<=K)'
    cleavage_patter_LysN = r'(?=K)'
    cleavage_patter_ArgC = r'(?<=R)'
    
    #Trypsin cleaves after K or R
    #only those peptides that are between 5 and 30 amino acids long are kept
    assert digest_protein_collection(dummy_proteins, cleave_pattern = cleavage_patter_trypsin) == {
        "Sequnece_ID_1": ["LLSEBN"],
        "Sequnece_ID_2": ["HSJER"]
    }
    assert digest_protein_collection(dummy_proteins, cleave_pattern = cleavage_patter_LysC) == {
        "Sequnece_ID_1": ["LLSEBN"],
        "Sequnece_ID_2": ["HSJERLSDE"]
    }
    assert digest_protein_collection(dummy_proteins, cleave_pattern = cleavage_patter_LysN) == {
        "Sequnece_ID_1": ["KLLSEBN"],
        "Sequnece_ID_2": ["HSJERLSDE"]
    }
    assert digest_protein_collection(dummy_proteins, cleave_pattern = cleavage_patter_ArgC) == {
        "Sequnece_ID_1": ["AKLLSEBN"],
        "Sequnece_ID_2": ["HSJER"]
    }

test_digest_protein_collection()

def test_compute_sequence_coverage():
    dummy_prot_seq = "APTSWQKLVFTYQK"
    dummy_peps_100 = ["APTSWQK", "LVFTYQK"]
    #hier sollten mehr Tests ergänzt werden

    
    assert compute_sequence_coverage(dummy_prot_seq, dummy_peps_100) == 100

test_compute_sequence_coverage()