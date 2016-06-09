from sequana import Coverage


def test_coverage():

    cover = Coverage(L=500, G=300000000)
    cover.__repr__() # test __repr__


    cover.a = 10
    assert cover.N == 6000000.0

    cover.N = 600000
    assert cover.a == 1
   
    cover  # Now that N and a are defined
    # one can cahnge L and G afterwards
    cover.L = 500
    cover.G=3e9 

    df = cover.get_table()

    cover.get_required_coverage(0.01)
    cover.get_required_coverage([1e-2,1e-3])


