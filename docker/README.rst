

For end-users
----------------

We will demonstrate the usage of the docker using the sequana_coverage
standalone. This application requires a BED file as input. An example is
available within the source code: https://raw.githubusercontent.com/sequana/sequana/master/sequana/resources/data/JB409847.bed)


A docker image is provided on `hub.docker <https://hub.docker.com/r/cokelaer/sequana/>`_.

We assume you have install docker on your system. First, get the docker image using::

    git pull cokelaer/sequana

You can then enter into the image as follows::

    docker run -it sequana

This opens an interactive shell with sequana v0.2 pre-installed. Have a go with
the test file called virus.bed::

    sequana_coverage --input virus.bed

This should print information and create a report/ directory. This is not very
practical if you have your own files or want to open the HTML page stored in
./reports. So, let us quit::

    exit

and do it the proper way. You need to tell docker where is your data on your
local machine. Let us say your data is in */home/user/mydatapath*. Start the docker as
follows::

    docker run -v /home/user/mydatapath:/home/sequana/data -it sequana

An alternative is to go to your local directory and type::

    cd  /home/user/mydatapath
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
