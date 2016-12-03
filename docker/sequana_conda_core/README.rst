Docker **sequana_conda_core**
====================================

This is a Docker image on top on sequana_conda_core that install conda and all sequana
dependencies. 



For developers:
------------------

Build the image::

    git clone https://github.com/sequana/sequana
    cd sequana/docker/sequana_conda_core
    sudo docker  build  -t="sequana_conda_core" .


Tag and push on hub.docker.com::

   docker tag sequana_conda_core sequana/sequana_conda_core
   docker push sequana/sequana_conda_core

Run the local image (not a pulled one in this example)::

    sudo docker run -it sequana_conda_core


