# How to build
docker build -t coldfire79/server .

# How to debug
docker run -v /Volumes/userdata/Joonyong/bsf-server-gse/:/tmp/server:consistent -p 127.0.0.1:3333:3333 -it --name server coldfire79/server /bin/bash

# How to run
docker run -v /Volumes/userdata/Joonyong/bsf-server-gse/:/tmp/server:consistent -it --name server coldfire79/bsf-py /bin/bash