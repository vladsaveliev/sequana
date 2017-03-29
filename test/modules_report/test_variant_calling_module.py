from sequana import sequana_data
from sequana.modules_report.variant_calling import VariantCallingModule
from sequana.utils import config


def test_variant_calling_module(tmpdir):
    directory = tmpdir.mkdir('test_variant_calling_module')
    config.output_dir = str(directory)
    config.sample_name = 'JB409847'
    VariantCallingModule(sequana_data('JB409847.vc.csv'))
