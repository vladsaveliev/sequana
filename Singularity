BootStrap: debootstrap
DistType "debian"
MirrorURL: http://us.archive.ubuntu.com/ubuntu/
OSVersion: xenial

%labels

    AUTHOR Thomas Cokelaer

%post

    apt-get install -y wget
    apt-get install -y bzip2
    apt-get install -y vim
    apt-get install -y libgl1-mesa-glx  # will be required by pyqt
    apt-get install -y fontconfig # for sequanix/qt fonts otherwise no text in menus

    # install anaconda
    if [ ! -d /usr/local/anaconda ]; then
        #wget https://repo.continuum.io/miniconda/Miniconda3-4.3.14-Linux-x86_64.sh\
        # for now, we use 4.2.12 to have python3.5 by default so no need to
        # create a new env saving space in the process. The reason for using 3.5
        # is inherent to the packages used at the moment.
        wget https://repo.continuum.io/miniconda/Miniconda3-4.2.12-Linux-x86_64.sh\
           -O ~/anaconda.sh && \
        bash ~/anaconda.sh -b -p /usr/local/anaconda && \
        rm ~/anaconda.sh
    fi

    # set anaconda path
    export PATH=$PATH:/usr/local/anaconda/bin

    conda config --add channels r
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

    #if [ ! -d /usr/local/anaconda/envs/sequana ]; then
    #    conda create --name sequana python=3.5
    #fi

    conda install --file https://raw.githubusercontent.com/sequana/sequana/master/requirements.txt
    conda install -y bwa fastqc kraken krona cutadapt
    conda install -y bowtie bowtie2 star subread
    conda install -y bcftools bedtools khmer samtools pigz bleach
    conda install -y snpeff freebayes spades multiqc sambamba


    conda clean --packages -y # next requires lots of space
    #conda install -y prokka  # takes a lot!!!

    #conda install -y busco==3.0.2
    conda install -y picard shustring
    conda install -y atropos<=1.1.10
    pip install sequana
    conda clean --all -y # next requires lots of space
    rm -rf /usr/local/anaconda/pkgs

    conda install --override-channels -c conda-forge bzip2
    conda install --override-channels -c bioconda -c conda-forge htslib==1.5.0

    if [ ! -d /data ]; then mkdir /data; fi
    if [ ! -d /scripts ]; then mkdir /scripts; fi
    if [ ! -d /scratch ]; then mkdir /scratch; fi
    if [ ! -d /mounting ]; then mkdir /mounting; fi
    if [ ! -d /pasteur ]; then mkdir /pasteur; fi

%environment
    export PATH=$PATH:/usr/local/anaconda/bin
    export LANG=C   # prevents perl for raising warnings
    export PERL5LIB=/usr/local/anaconda/lib/perl5/5.22.2/
    echo "backend:agg" > matplotlibrc

