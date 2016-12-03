Docker **sequana_core**
====================================

This is a Docker image with only the Linux system requirements (no conda, no
sequana package here). This can be used later a the core layer for all Sequana



For developers:
------------------

Build the image::

    git clone https://github.com/sequana/sequana
    cd sequana/docker/sequana_core
    sudo docker  build  -t="sequana_core" .


Tag and push on hub.docker.com::

   docker tag sequana_core sequana/sequana_core
   docker push sequana/sequana_core

Run the local image (not a pulled one in this example)::

    sudo docker run -it sequana_core


