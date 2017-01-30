from sequana.iotools import YamlDocParser
from sequana import sequana_data
from sequana import snaketools



def test_yamlreader():
    filename = sequana_data("test_gui_generic_config.yaml")
    r = YamlDocParser(filename)
    assert r.sections['N'] == '# example of docstring\n'

    # check docstring is parsed with #### removed
    r = YamlDocParser(snaketools.Module('quality_control').config)
    docstring = r._block2docstring("fastqc")
    assert docstring.startswith("FastQC")

    assert r._block2docstring("dummy") is None

    # check that pipelines can be parsed
    for pipeline in snaketools.pipeline_names:
        filename = snaketools.Module(pipeline).config
        r = YamlDocParser(filename)
