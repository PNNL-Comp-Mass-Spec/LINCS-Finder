[![Docker Build Status](https://img.shields.io/docker/build/coldfire79/lincs-finder.svg)](https://hub.docker.com/r/coldfire79/lincs-finder/)

# LINCS-Finder: a Python web app for finding LINCS signatures
This is a web application to quickly find the most interesting LINCS signatures. It's based on Blazing Signature Filter (BSF) behind. Please refer to https://github.com/PNNL-Comp-Mass-Spec/bsf-py/ for details.

# How to run
We provide a docker version for this application. Please refer to these following steps to run this app on your local machine.

1. Install the latest docker in your machine.
Please refer to https://docs.docker.com/engine/installation/.
2. Pull the LINCS Finder image.
```bash
docker pull coldfire79/lincs-finder
```
3. Run the LINCS Finder image.
```bash
docker run -p 127.0.0.1:8081:8081 -it --name lincs-finder coldfire79/lincs-finder
```
3. Go to 'http://localhost:8081/'

# How to use
Please submit your own list of genes into corresponding text boxes to query the similar signatures. Then, click the `Submit` button.

# How to open KEGG pathways
For example, http://localhost:8081/pathway.html?id=hsa05145&nodes=LJP005_A375_24H:G01.
You can try to open your own links to display a KEGG pathway of interest as follows:
localhost:8081/pathway.html?id=\<KEGG ID\>&nodes=\<LINCS ID\>.
