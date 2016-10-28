from sequana.sphinxext import snakemakerule



@skip('!travis')
def test_doc():
    res  = snakemakerule.get_rule_doc("dag")
    res  = snakemakerule.get_rule_doc("fastqc")

    try:
        res  = snakemakerule.get_rule_doc("dummy")
        assert False
    except FileNotFoundError:
        assert True
    except:
        assert False


    from sphinx.application import Sphinx
    app = Sphinx(".", "/home/cokelaer/Work/github/sequana/doc/", ".", ".", "html")
