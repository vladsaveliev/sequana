Docker **sequana_conda_ngs**
====================================

This is a Docker image on top on sequana_conda_ngs that install conda and all sequana
dependencies. 



For developers:
------------------

Build the image::

    git clone https://github.com/sequana/sequana
    cd sequana/docker/sequana
    sudo docker  build  -t="sequana" .


Tag and push on hub.docker.com::

   docker tag sequana sequana/sequana
   docker push sequana/sequana

Run the local image (not a pulled one in this example)::

    sudo docker run -it sequana


