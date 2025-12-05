from __future__ import annotations

from pathlib import Path

from DAJIN2.core.preprocess.structural_variants import sv_handler


def test_add_unique_allele_keys_with_internal_suffix_and_mapping(tmp_path):
    fasta_sv_alleles = {"sv1": "AAA", "sv2": "TTT"}
    fasta_alleles_user = {"control": "AAA", "deletion01": "CCC"}

    internal, name_map = sv_handler.add_unique_allele_keys(
        fasta_sv_alleles, fasta_alleles_user, key="deletion", internal_suffix="uuid123"
    )

    assert set(internal.keys()) == {"deletion001__uuid123", "deletion002__uuid123"}
    assert internal["deletion001__uuid123"] == "AAA"
    assert name_map["deletion001__uuid123"] == "DAJIN_deletion001"

    sv_handler.save_sv_name_map(tmp_path, "sample1", name_map)
    loaded = sv_handler.load_sv_name_map(tmp_path, "sample1")
    assert loaded == name_map
    inverted = sv_handler.invert_sv_name_map(loaded)
    assert inverted["DAJIN_deletion001"] == "deletion001__uuid123"
    assert inverted["DAJIN_deletion002"] == "deletion002__uuid123"
