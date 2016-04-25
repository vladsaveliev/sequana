


def stacked_bar(tag, title, datalist):
    """

    data list should be a list of dictionary formatted as follows

    {"name": "A"
     "data": {
        "R1_mapped": 50,
        "R2_mapped": 50,
        "R1_unmapped": 50,
        "R2_unmapped": 50,
        }
    }
    """
    

    dataitems = ""
    for item in datalist:
        datatext = []
        for k,v in item['data'].items():
            datatext.append('{y:%s,label:"%s"}' % (v,k))
        datatext = ",\n        ".join(datatext)

        params = {
            "name": item['name'], 
            "datatext": datatext}

        dataitems += """
      {
        type: "stackedBar100",
        showInLegend: true,
        name: "%(name)s",
        dataPoints: [
        %(datatext)s
        ]
        },
    """ % params

    metadata = {
        'tag': tag,
        'title': title,
        'dataitems': dataitems}

    script = """
<script type="text/javascript">
  window.onload = function () {
    var chart = new CanvasJS.Chart("chartContainer%(tag)s",
    {
      theme: "theme2",
      title:{
        text: "%(title)s"
      },
      animationEnabled: true,
      axisY:{
        title: "percent"
      },
      legend :{
        horizontalAlign: 'center',
        verticalAlign: 'bottom'
      },
      toolTip: {
        shared: true
      },
      data:[
        %(dataitems)s
      ]

    });

chart.render();
}
</script>
"""
    return script % metadata

#<script type="text/javascript" src="js/canvasjs.min.js"></script></head>
#<body>
#  <div id="chartContainer" style="height: 300px; width: 100%;">
#  </div>
#</body>
#</html>


