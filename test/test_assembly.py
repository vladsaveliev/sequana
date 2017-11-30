from sequana.assembly import BUSCO
from sequana import sequana_data


def test_busco():
    filename = sequana_data('test_busco_full_table.tsv')
    b = BUSCO(filename)
    print(b)
    b.pie_plot()
    b.scatter_plot()
    assert b.score > 90
    assert b.score <91
