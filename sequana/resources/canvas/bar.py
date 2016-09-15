


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


class CanvasBar(object):
    def __init__(self, data, title="", tag="", xlabel="",
        ylabel_max_length=20, **kargs):

        self.title = title
        self.tag = tag

        if isinstance(data, dict):
            dataitems = "["
            for label, value in data.items():
                dataitems += '{y:%s,label:"%s"},' % (value, label)
            dataitems += "]"
        else:
            raise NotImplementedError

        self.metadata = {
            "dataitems": dataitems,
            "title": self.title,
            "xlabel":xlabel,
            "tag": self.tag}

    def to_html(self):
        script = """
    var chart = new CanvasJS.Chart("chartContainer%(tag)s",
    {
      theme: "theme2",
      title:{
        text: "%(title)s"
      },
      animationEnabled: true,
      axisY:{
        title: "%(xlabel)s"
      },
      legend :{
        horizontalAlign: 'center',
        verticalAlign: 'bottom'
      },
      toolTip: {
        shared: true
      },
      data:[
      {
          type: "bar",
          dataPoints: %(dataitems)s
      },
      ]
    });

chart.render();

"""
        return script % self.metadata
#<div id="chartContainer%(tag)s" style="height: 300px; width: 100%%;">

