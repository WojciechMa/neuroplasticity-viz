import pytest

from fi_analysis.cell_class import classify_ttype


@pytest.mark.parametrize(
    "label, expected",
    [
        ("Pvalb Reln Itm2a", "I"),
        ("Sst Calb2 Pdlim5", "I"),
        ("Vip Rspo4 Rxfp1 Chat", "I"),
        ("Lamp5 Lsp1", "I"),
        ("L2/3 IT VISp Agmat", "E"),
        ("L5 PT VISp C1ql2 Ptgfr", "E"),
        ("Inh L1-2 PAX6 CDH12", "I"),
        ("Exc L2-3 LINC00507 FREM3", "E"),
        ("Astro Aqp4", None),
        ("", None),
        (None, None),
    ],
)
def test_classify_ttype(label, expected):
    assert classify_ttype(label) == expected
