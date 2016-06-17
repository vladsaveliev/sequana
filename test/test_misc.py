from sequana.misc import moving_average


def test_moving_average():

    moving_average([1,2,3,4,2,1,2,4,3,4], 3)
