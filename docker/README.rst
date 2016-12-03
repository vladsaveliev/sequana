Docker containers for **Sequana**
====================================

Docker containers wrap a piece of software in a complete filesystem that contains everything needed to run (code, system libraries). We assume you have installed Docker on your system (see  `Docker <https://www.docker.com>`_ otherwise).

In order to allows anyone to use **Sequana** without needs for complex installation, we provide 
`Docker images <https://hub.docker.com/u/sequana>`, which are synchronized on the *master* 
branch of the source code.


Quick start
----------------
Since images are on `Hub Docker <https://www.hub.docker.com>`_ , you will first need to 
sign up to be able to login as follows::

    docker login

**Sequana** itself (all library, pipelines and standalones) is available and can be downloaded as follows (1.5Gb)::

    # Get the docker image
    docker pull sequana/sequana

Now, you can use it as a developer to access for instance to executables,
pipelines or the Python library itself::

    cd <Directory_with_bed_files>
    docker run -v $PWD:/home/sequana/data -it sequana/sequana

Standalone
----------------

The primary goal of the docker is to make it possible to quickyl test the
standalones. For now, we expose only one but more will come. For specific
standalone, please see links here below:

- sequana_coverage_

.. _sequana_coverage: sequana_coverage/README.rst


Usage
---------------------------

Here below, we provide a quick tutorial that will guide you on using **Sequana**
thanks to the docker. To do so, we will focus on one standalone application
called **sequana_coverage**. In brief, the standalone takes as input a BED file
that contains the genome coverage of a set of mapped DNA reads onto a reference
genome. Then, the standalone creates a report with relevant information about
the coverage (See `Sequana documentation <sequana.readthedocs.org>`_ for 
more information).

Use the **sequana** Docker image
---------------------------------------

Once you downloaded the **sequana** image, you can then enter into the image as follows::

    docker run -it sequana

This opens an interactive shell with sequana v0.2 pre-installed. Have a go with
the test file called virus.bed::

    sequana_coverage --input virus.bed

This should print information and create a report/ directory. This is not very
practical if you have your own files or want to open the HTML page stored in
./reports. So, let us quit::

    exit

and do it the proper way. Go to a working directory and start the docker image again as
follows::

    docker run -v $PWD:/home/sequana/data -it sequana

This should start the docker image again but you should now have a *./data*
directory. **Be aware that if you modify data here (in the image),
you will also modify the data in your local data file.**

Now, you can run sequana_coverage in this directory::

   sequana_coverage --input yourfile.bed

This analyses the data and creates a report/ directory. The container has no
display but you can now go back to your computer in /home/user/mydatapath and
browse the HTML page that was created.


For developers:
------------------

Build the image::

    git clone https://github.com/sequana/sequana
    cd sequana/docker
    sudo docker  build  -t="sequana" .

Run the image::

    sudo docker run -it sequana


To avoid using sudo, check out various forum. Seefor example:  http://askubuntu.com/questions/477551/how-can-i-use-docker-without-sudo
