

For end-users
----------------

In order to allows anyone to use **Sequana** without needs for complex installation, we provide a
`Docker <https://www.docker.com/>` image. It is synchronized on the *master*
branch on the source code, which means on official releases.


Here below, we provide a quick tutorial that will guide you on using **Sequana**
thanks to the docker. To do so, we will focus on one standalone application
called **sequana_coverage**. In bried, the standalone takes as input a BED file
that contains the genome coverage of a set of mapped DNA reads onto a reference
genome. Then, the standalone creates a report with relevant information about
the coverage (See `Sequana documentation <sequana.readthedocs.org>`_ for details.

Get the docker image
-------------------------

A docker image is provided on `hub.docker <https://hub.docker.com/r/cokelaer/sequana/>`_.

We assume you have install docker on your system.


First, you need a login on `docker <hub.docker.com>`_. Then, in a shell, type::

    docker login

This will allow you to obtain the **Sequana** docker image using::

    docker pull cokelaer/sequana

Let us rename it into *sequana*::

    docker tag cokelaer/sequana sequana
    sudo docker rmi cokelaer/sequana

You can then enter into the image as follows::

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


sudo docker run -v $PWD:/tempdir -it sequana3
exit


.. seealso:: to avoid sudo
    http://askubuntu.com/questions/477551/how-can-i-use-docker-without-sudo
