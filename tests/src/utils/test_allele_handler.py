from __future__ import annotations

import json

from DAJIN2.utils.allele_handler import save_allele_name_map, to_allele_key


def test_to_allele_key_returns_md5_hex():
    assert to_allele_key("control") == "fc5364bf9dbfa34954526becad136d4b"


def test_save_allele_name_map(tmp_path):
    path_output = save_allele_name_map(tmp_path, ["control", "deletion"])
    payload = json.loads(path_output.read_text(encoding="utf-8"))

    assert path_output == tmp_path / "cache" / "allele_name_map.json"
    assert payload["name_to_key"]["control"] == to_allele_key("control")
    assert payload["key_to_name"][to_allele_key("deletion")] == "deletion"
