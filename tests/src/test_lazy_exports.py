import importlib
import sys


def unload_modules(*module_names: str) -> None:
    for module_name in module_names:
        sys.modules.pop(module_name, None)


def test_preprocess_package_defers_export_modules_until_attribute_access():
    unload_modules(
        "DAJIN2.core.preprocess",
        "DAJIN2.core.preprocess.alignment.mapping",
        "DAJIN2.core.preprocess.infrastructure.directory_manager",
    )

    preprocess = importlib.import_module("DAJIN2.core.preprocess")

    assert "DAJIN2.core.preprocess.alignment.mapping" not in sys.modules
    assert "DAJIN2.core.preprocess.infrastructure.directory_manager" not in sys.modules
    assert preprocess.create_temporal_directories.__name__ == "create_temporal_directories"
    assert "DAJIN2.core.preprocess.infrastructure.directory_manager" in sys.modules


def test_report_package_defers_bam_exporter_until_attribute_access():
    unload_modules(
        "DAJIN2.core.report",
        "DAJIN2.core.report.bam_exporter",
        "DAJIN2.core.report.sequence_exporter",
    )

    report = importlib.import_module("DAJIN2.core.report")

    assert "DAJIN2.core.report.bam_exporter" not in sys.modules
    assert "DAJIN2.core.report.sequence_exporter" not in sys.modules
    assert report.__all__ == ["bam_exporter", "sequence_exporter", "vcf_exporter"]
    assert "sequence_exporter" in dir(report)
