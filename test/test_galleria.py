from sequana.htmltools import galleria


def test_galleria():
    html = galleria('test')
    assert "galleria" in html
