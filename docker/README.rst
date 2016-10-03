
In progress

For end-users
----------------

Assuming you have nice data (e.g., a BED file) and want to try the
sequana_coverage standalone quickly, you can use a docker image available on
`hub.docker <https://hub.docker.com/r/cokelaer/sequana/>`_

We assume you have install docker on your system.

First, get the docker image using::

    git pull cokelaer/sequana

You can test it as follows::

    docker run -it sequana

This opens an interactive shell with sequana v0.2 pre-installed. Have a go with
the test file called virus.bed::

    sequana_coverage --input virus.bed

This should print information and create a report/ directory. This is not very
practical if you have your own files or want to open the HTML page stored in
./reports. So, type::

    exit

and let us now do it the proper way. You need to tell docker where is you data.
Let us say your data is in */home/user/mydatapath*. Start the docker as
follows::

    docker run -v /home/user/mydatapath:/home/sequana/data -it sequana

An alternative is to go to your local directory and type::

    cd  /home/user/mydatapath
    docker run -v $PWD:/home/sequana/data -it sequana

This should start the docker image again but you should now have a *./data*
directory. Be aware that if you modify data here (in the image),
you will also modify the data in your local data file.

Now, you can run sequana_coverage in this directory::

   sequana_coverage --input yourfile.bed

This analyses the data and creates a report/ directory. You can exit and 
the modifications are persistent. So you can see *./report* in yout local 
directory and in particular your The container has no
display but you can now go back to your computer in /home/user/mydatapath and
browse the HTML page that was created.

You can also exit the docker interactive shell.



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
