from sequana.vcf_filter import VCF
from sequana import sequana_data
from easydev import TempFile



def test_vcf_filter_freebayes():

    data = sequana_data("test.vcf")
    v = VCF(data)
    v.hist_qual()

def test_vcf_filter():

    data = sequana_data("test_vcf_mpileup_4dot1.vcf")
    v = VCF(data)
    assert v.vcf.version == "4.1"


    # Test the INFO set to {}
    v.vcf.filter_dict['INFO'] = {}
    v.vcf.filter_dict['QUAL'] = 50
    with TempFile() as fh:
        res = v.vcf.filter_vcf(fh.name)
    assert res == {'N': 573, 'filtered': 308, 'unfiltered': 265}


    # Test the & filter
    v.vcf.filter_dict["QUAL"] = 0
    v.vcf.filter_dict["INFO"] = {}
    v.vcf.rewind();
    v.vcf.filter_dict['INFO']['DP'] = ">40&<80"
    with TempFile() as fh:
        res = v.vcf.filter_vcf(fh.name)
    assert res == {'N': 573, 'filtered': 10, 'unfiltered': 563}


    # check the | filter
    v.vcf.filter_dict["QUAL"] = 0
    v.vcf.filter_dict["INFO"] = {}
    v.vcf.rewind();
    v.vcf.filter_dict['INFO']['DP'] = "<40|>80"
    with TempFile() as fh:
        res = v.vcf.filter_vcf(fh.name)
    assert res ==  {'N': 573, 'filtered': 562, 'unfiltered': 11}


    # check the array filter
    v.vcf.filter_dict["QUAL"] = 0
    v.vcf.filter_dict["INFO"] = {}
    v.vcf.rewind();
    v.vcf.filter_dict['INFO']['DP4[0]'] = "<2"
    with TempFile() as fh:
        res = v.vcf.filter_vcf(fh.name)
    assert res ==  {'N': 573, 'filtered': 173, 'unfiltered': 400}

    # Check the array access to one item
    v.vcf.filter_dict["QUAL"] = 0
    v.vcf.filter_dict["INFO"] = {}
    v.vcf.rewind();
    v.vcf.filter_dict['INFO']['DP4[2]'] = "<2"
    with TempFile() as fh:
        res = v.vcf.filter_vcf(fh.name)
    assert res ==  {'N': 573, 'filtered': 199, 'unfiltered': 374}

    # Check the sum of elements in an array
    v.vcf.filter_dict["QUAL"] = 0
    v.vcf.filter_dict["INFO"] = {}
    v.vcf.rewind();
    v.vcf.filter_dict['INFO']['sum(DP4[2], DP4[3])'] = "<4"
    with TempFile() as fh:
        res = v.vcf.filter_vcf(fh.name)
    assert res ==  {'N': 573, 'filtered': 7, 'unfiltered': 566}
