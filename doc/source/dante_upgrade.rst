sequoia provides snakefiles and in particular, a snakemake file that mimics the
behaviour of DANTE. There are however some differences. Some of which are
summarized here below.


#. **bwa** is run with the -M option possibly
# cutadapt in DANTE is ran on each reads R1 and R2 indepdendently. Reports of
  parsing of the 2 log files is done twice on each log file. In Sequoia, we 
  run the cutadapt on the pair of reads together. The parsing of the cutadapt
  log is therefore different.




