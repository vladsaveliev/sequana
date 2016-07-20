from sequana.stats import moving_average, evenness


def test_ma():
    ma = moving_average([1,1,1,1,3,3,3,3], 4)
    assert list(ma) == [ 1. ,  1.5,  2. ,  2.5,  3. ]

def test_evenness():
    assert evenness([1,1,1,1,4,4,4,4]) == 0.75
