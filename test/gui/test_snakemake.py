from sequana.gui import snakemake

def test_snakemakedialog(qtbot):
    form = snakemake.SnakemakeDialog()
    qtbot.addWidget(form)

    form.get_snakemake_local_options()
    form.get_snakemake_general_options()
    form.get_snakemake_cluster_options()

    form.accept()

    form.reject()
    form.get_settings()

