from sequana.sphinxext import snakemakerule



def test_doc():
    res  = snakemakerule.get_rule_doc("dag")
    assert len(res)
