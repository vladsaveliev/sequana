from sequana.modules_report import pacbio_input_bam



def test_pacbio_input_bam():
    try:
        # we need a summary and a bunch of images
        ff = pacbio_input_bam.PacbioInputBAMModule()
    except:
        pass

