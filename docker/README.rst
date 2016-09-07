
In progress


Build the image::

    git clone https://github.com/sequana/sequana
    cd sequana/docker
    sudo docker  build  -t="sequana" .


Run the image::

    sudo docker run -i -t -entrypoint='/root' sequana

