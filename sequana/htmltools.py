


def galleria(filenames):

    html = """ <div class="galleria">"""
    for filename in filenames:
        html += """<img src="%s" 
                        data-title="%s" 
                        data-description="%s">\n""" % (filename, filename, filename)
    html += "</div>"
    html += """
 <script>
     Galleria.loadTheme('galleria/themes/classic/galleria.classic.min.js');
     Galleria.run('.galleria');
 </script>
    """

    return html

