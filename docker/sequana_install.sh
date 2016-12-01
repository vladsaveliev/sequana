cd /home/sequana

apt-get update
apt-get install wget -y
apt-get install bzip2 -y
apt-get -y install gcc # required to compile the cython code
apt-get install libxtst6 -y
apt-get install -y python-qt4
apt-get install git -y
apt-get install vim -y

apt-get install libpng-dev -y
apt-get install zlib1g-dev -y
apt-get install make -y
apt-get install xvfb -y
apt-get install firefox -y

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod 755 Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh -b -f -p /home/sequana/miniconda3

source /home/sequana/miniconda3/bin/activate

conda config --add channels r
conda config --add channels bioconda

conda install numpy matplotlib pandas cutadapt pysam pyvcf snpeff -y
conda install snakemake biokit bioservices spades khmer -y
conda install bwa bcftools samtools==1.3.0 bedtools picard freebayes fastqc -y
conda install kraken krona scipy graphviz pigz bowtie2 -y
conda install bowtie bowtie2 fastq-screen subread tophat -y

conda install ipython -y

# Finally, Sequana release itself
pip install sequana==0.1.14

#git clone https://github.com/sequana/sequana

#
cp /home/sequana/miniconda3/lib/python3.5/site-packages/sequana/resources/data/Hm2_GTGAAA_L005_R* .
cp /home/sequana/miniconda3/lib/python3.5/site-packages/sequana/resources/data/virus.bed  .

# change matplotlibrc to switch QtAgg to Agg. Nov 2016 seems that now, we have
# PyQt5 in anaconda by default. Let us keep Qt4 sed as well in case of
cp /home/sequana/miniconda3/lib/python3.5/site-packages/matplotlib/mpl-data/matplotlibrc .
sed -i -e 's/Qt4Agg/Agg/'g matplotlibrc
sed -i -e 's/Qt5Agg/Agg/'g matplotlibrc

