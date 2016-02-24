
Quality
===============

wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.4.zip
unzip fastqc_v0.11.4.zip
cd FastQC
chmod 755 fastqc
cd ..
# rm fastqc_v0.11.4.zip
cd bin
ln -s fastqc ../FastQC/fastqc
cd ..




Adapter Removal
=====================

cutadapt
-------------

Installation::

    pip install cutadapt


btrim
--------
The code is here: http://graphics.med.yale.edu/trim/ . Is it maintained ? 
Seems fast but last version is 2014.


Alignment 
===========

bwa
---------

C compilation but no dependencies so should work out of the box::

    git clone https://github.com/lh3/bwa
    cd bwa 
    make

    cp ./bwa $(BIN)    



















