# instructions to use dreamtools inside a Fedora docker (22)

FROM ubuntu:16.04
MAINTAINER Thomas Cokelaer


#ENTRYPOINT ["/home/sequana"]
#CMD ["sequana --help"]
#EXPOSE 8100






# Just copy the installation scripts and execute it
ADD sequana_install.sh /
RUN /bin/bash /sequana_install.sh

ENV HOME /home/sequana
