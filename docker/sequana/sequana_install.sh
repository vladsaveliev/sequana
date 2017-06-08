cd /home/sequana

# Finally, Sequana release itself
pip install sequana==0.3.0
pip install line_profiler


# Some dummy data to play with
cp /home/sequana/miniconda3/lib/python3.5/site-packages/sequana/resources/data/Hm2_GTGAAA_L005_R* .
cp /home/sequana/miniconda3/lib/python3.5/site-packages/sequana/resources/data/virus.bed  .


# This build the .config file once for all and avoid warning message
ipython -c "import sequana; import pylab; pylab.plot([1,1])"
