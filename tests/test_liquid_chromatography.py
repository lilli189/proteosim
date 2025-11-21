from proteosim.liquid_chromatography import predict_lc_retention_times, select_retention_time_window

def test_predict_lc_retention_times():
    peptides = ["DMMY", "PEPTIDE", "LIST", "TEST"]
    expected = {'DMMY': 14.0, 'PEPTIDE': 7.8, 'LIST': 17.1, 'TEST': -1.2}

    actual = predict_lc_retention_times(peptides)
    assert actual == expected

test_predict_lc_retention_times()

def test_select_retention_time_window():
    peptide_rt_map = {'MATSR': 8.9, 'YEPVAEIGVGAYGTVYK': 42.7, 'DPHSGHFVALK': 29.0}
    selected = select_retention_time_window(peptide_rt_map, lower_ret_time=10, upper_ret_time=30)

    assert selected == {'DPHSGHFVALK': 29.0}

test_select_retention_time_window()