from sequana.plots.canvasjs_linegraph import CanvasJSLineGraph
from sequana import bedtools, sequana_data


def test_canvasjs_linegraph():
    bed = bedtools.GenomeCov(sequana_data("JB409847.cov.csv"))
    df = bed[0].df
    csv = df.to_csv(columns=['pos', 'cov', 'gc'], index=False,
                    float_format='%.3g')
    # create CanvasJS stuff
    cjs = CanvasJSLineGraph(csv, 'cov', 'pos', ['cov', 'gc'])
    # set options
    cjs.set_options({'zoomEnabled': 'true',
                     'zoomType': 'x',
                     'exportEnabled': 'true'})
    # set title
    cjs.set_title("Genome Coverage")
    # set legend
    cjs.set_legend({'verticalAlign': 'bottom',
                    'horizontalAlign': 'center',
                    'cursor':'pointer'},
                    hide_on_click=True)
    # set axis
    cjs.set_axis_x({'title': "Position (bp)",
                    'labelAngle': 30,
                    'minimum': 0,
                    'maximum': len(df)})
    cjs.set_axis_y({'title': "Coverage (Count)"})
    cjs.set_axis_y2({'title': "GC content (ratio)",
                     'minimum':0,
                     'maximum': 1,
                     'lineColor': '#FFC425',
                     'titleFontColor': '#FFC425',
                     'labelFontColor': '#FFC425'})
    # set datas
    cjs.set_data(index=0, data_dict={'type': 'line',
                                     'name': "Coverage",
                                     'showInLegend': 'true',
                                     'color': '#5BC0DE',
                                     'lineColor': '#5BC0DE'})
    cjs.set_data(index=1, data_dict={'type': 'line',
                                     'axisYType': 'secondary',
                                     'name': "GC content",
                                     'showInLegend': 'true',
                                     'color': '#FFC425',
                                     'lineColor': '#FFC425'})
    # create canvasJS
    cjs.create_canvasjs()
