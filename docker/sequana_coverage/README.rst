**sequana_coverage** docker
===============================


This docker creates a single entry point for the sequana_coverage **standalone**. It is a simple layer on top 
of the **sequana** image, which makes available the **sequana_coverage** standalone only. For the full library, please see the main docker image.

#. Make sure you are logged in hub.docker.com:

    docker login

#. Get the docker image and rename it::

    docker pull cokelaer/sequana:sequana_coverage
    
    # Rename it (optional but simplifies the code hereafter)
    docker tag cokelaer/sequana sequana
    sudo docker rmi cokelaer/sequana

#. To get some help::

    docker run --user 1000 -v $(pwd):/home/sequana/data -w /home/sequana/data -it sequana_coverage

#. Analyse a BED file and create a report

Assuming that:

#. you have a BED file in your local directory (pwd)
#. Your bed file is named *JB409847.bed*::

    docker run --user 1000 -v $(pwd):/home/sequana/data -w /home/sequana/data -it sequana_coverage --input *JB409847.bed*

You may want to create an alias::

    alias sequana_cov = "docker run --user 1000 -v $(pwd):/home/sequana/data -w /home/sequana/data -it sequana_coverage"




For developers
-----------------

Create the docker as follows::

    sudo docker build -t="sequana_coverage" .

It is a layer on top of the **sequana** docker.




