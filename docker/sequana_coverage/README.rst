**sequana_coverage** docker
===============================


This docker creates a single entry point for the sequana_coverage **standalone**.



Get the docker image::

    docker pull cokelaer/sequana_coverage



Assuming you have a BED file in your local directory (pwd), type::

    docker run --user 1000 -v $(pwd):/home/sequana/data -w /home/sequana/data -it sequana_coverage

to get some help, or type::

    docker run --user 1000 -v $(pwd):/home/sequana/data -w /home/sequana/data -it sequana_coverage --input *JB409847.bed*

to analyse the local file named *JB409847.bed*



For developers
-----------------

Create the docker as follows::

    sudo docker build -t="sequana_coverage" .

It is a layer on top of the **sequana** docker.




