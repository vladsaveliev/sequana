# instructions to use dreamtools inside a Fedora docker (22)

FROM ubuntu:16.04
MAINTAINER Thomas Cokelaer <thomas.cokelaer@pasteur.fr>

# Setup home environment
RUN useradd sequana
RUN mkdir /home/sequana && chown -R sequana: /home/sequana
RUN echo $HOME


# Instructions after ADD are not cached
ADD sequana_install.sh .
RUN /bin/bash sequana_install.sh

#ENV GOPATH /home/dev/go:$GOPATH
#export PYTHONPATH=$HOME/miniconda3/lib/python3.5/site-packages/
#export PATH=$PATH:$HOME/miniconda3/bin/


WORKDIR /home/sequana


ENV PATH $PATH:/home/sequana/miniconda3/bin
ENV MATPLOTLIBRC=/home/sequana
CMD ["/bin/bash"]




