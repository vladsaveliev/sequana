from sequana import sequana_data


def test_sequana_data():

    f1 = sequana_data("Institut_Pasteur.png")
    f2 = sequana_data("Institut_Pasteur.png", "images")
    assert f1 == f2
  
