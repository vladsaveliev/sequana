Docker containers for **Sequana**
====================================

.. warning:: Although we provide a Docker recipes, this method will not be
    maintained after release 0.3.0 of Sequana. However, we keep the version 0.3
    and the following recipes here below for book-keeping and those willing to build
    their own docker of Sequana


`Docker <http://www.docker.com>`_ containers wrap a piece of software in a complete filesystem that contains everything needed to run the software.

In order to allow anyone to use **Sequana** without needs for complex installation, we provide
`Docker images <https://hub.docker.com/u/sequana>`_, which are synchronized on the *master*
branch of the source code.

We assume that:

#. You have installed Docker on your system (see  `Docker <https://www.docker.com>`_ otherwise).
#. You have an account on  `Hub Docker <https://hub.docker.com>`_ .


Quick start
----------------
With your hub.docker account, first login::

    docker login

Then download (pull) a **Sequana** image (all library, pipelines and standalones) as follows (2Gb image in total)::

    docker pull sequana/sequana

Now, you should be ready to try it. To start an interactive session, type::

    cd <Directory_with_data>
    docker run -v $PWD:/home/sequana/data -it sequana/sequana
    

Standalone
----------------

The primary goal of the docker is to make it possible to quickly test the
standalones. For now, we expose only one docker. Please see specific 
documentation following the links here below:

- sequana_coverage: (https://github.com/sequana/sequana/tree/master/docker/sequana_coverage)
- sequana_taxonomy: (https://github.com/sequana/sequana/tree/master/docker/sequana_taxonomy)


More advanced Usage
---------------------------

Here below, we provide a quick tutorial that will guide you on using **Sequana**
thanks to the docker. To do so, we will focus on one standalone application
called **sequana_coverage**. In brief, the standalone takes as input a BED file
that contains the genome coverage of a set of mapped DNA reads onto a reference
genome. Then, the standalone creates a report with relevant information about
the coverage (See `Sequana documentation <http://sequana.readthedocs.org>`_ for
more information).

Use the **sequana** Docker image
---------------------------------------

Once you downloaded the **sequana** image, you can then enter into the image as follows::

    docker run -it sequana/sequana

This opens an interactive shell with latest sequana library pre-installed. For instance, you can
start an IPython shell::

    ipython
 
and import the library::

    import sequana


Or within the unix shell, you can use standalones. For instance there is a test
BED file that can be analysed as follows to get a coverage report::

    sequana_coverage --input virus.bed

This should print information and create a report/ directory. This is not very
practical if you have your own files or want to open the HTML page stored in
./report. So, let us quit the docker::

    exit

and do it the proper way. Go to a working directory (or your computer )and start the 
docker image again as follows::

    docker run -v $PWD:/home/sequana/data -it sequana/sequana

This should start the docker image again but you should now have a *./data*
directory. **Be aware that if you modify data here (in the image),
you will also modify the data in your local data file.**

Now, you can run sequana_coverage in this directory::

    cd data
    sequana_coverage --input yourfile.bed

This analyses the data and creates a report/ directory. The container has no
display but you can now go back to your computer in /home/user/mydatapath and
browse the HTML page that was created.

Each time, we entered in the image but you can also use the images as
executables (see standalone section above).


For developers:
------------------


Build the image::

    git clone https://github.com/sequana/sequana
    cd sequana/docker/sequana_core
    sudo docker  build  -t="sequana/sequana_core" .

Run the image::

    sudo docker run -it sequana/sequana_core


Layers
~~~~~~~~~~~
Here are the layers made available on hub.docker.com/u/sequana organizations.
Each layer is built on top of the previous one

- sequana_core_  (only ubuntu + some packages)
- sequana_conda_core_ (sequana_core + conda + common scientific packages)
- sequana_conda_ngs_ (sequana_conda_core + NGS conda packages)
- sequana_ (sequana_conda_ngs + sequana specific version)
- Standalone Layers:

  - sequana_coverage_ (sequana + sequana_coverage standalone)

.. _sequana_core: https://github.com/sequana/sequana/tree/master/docker/sequana_core
.. _sequana_conda_core: https://github.com/sequana/sequana/tree/master/docker/sequana_conda_core
.. _sequana_conda_ngs: https://github.com/sequana/sequana/tree/master/docker/sequana_conda_ngs
.. _sequana: https://github.com/sequana/sequana/tree/master/docker/sequana
.. _sequana_coverage: https://github.com/sequana/sequana/tree/master/docker/sequana_coverage



Sudo
~~~~~~~~~

To avoid using sudo, check out various forum. See for example:  http://askubuntu.com/questions/477551/how-can-i-use-docker-without-sudo
