# instructions to use dreamtools inside a Fedora docker (22)

FROM ubuntu:16.04
MAINTAINER Thomas Cokelaer


# Instructions after ADD are not cached
ADD sequana_install.sh /
RUN /bin/bash /sequana_install.sh

# the source is done in the previous RUN command
# RUN source conda/bin/activate
CMD ["sequana_coverage"]




