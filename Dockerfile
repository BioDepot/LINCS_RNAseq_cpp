FROM biodepot/alpine-bwa:3.7-0.7.15 
RUN apk update upgrade --no-cache && apk add bash musl-dev libgcc libgomp libstdc++  boost-filesystem
COPY w384/ /384/
COPY w96/ /96/
ADD scripts/multibwa.sh /usr/local/bin/multibwa.sh
ADD scripts/start.sh start.sh
ENV NWELLS 96
ENTRYPOINT ["/start.sh"]
