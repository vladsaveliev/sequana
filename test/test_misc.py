from sequana.misc import *


def test_textwrap():
    res = textwrap("test1test2", width=5, indent=0)
    assert res.split("\n")[0] == "test1"
    assert res.split("\n")[1] == "test2"

def findpos():
    pass

def wget():
    pass
