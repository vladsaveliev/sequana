# Sequana docker builder
FROM sequana
MAINTAINER Thomas Cokelaer <thomas.cokelaer@pasteur.fr>


# Make sure the backend works
ENV MPLBACKEND "Agg"

ENTRYPOINT ["sequana_coverage"]
CMD ['--help']
