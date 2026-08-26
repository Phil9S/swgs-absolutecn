#!/bin/sh
# argo workflow deployment for swgs-absolutecn
# Requires Docker, kubectl & minikube installed
## Run the following
minikube start
nohup minikube mount /home/psmith/swgs-absolutecn/data/:/mnt/data/ &
kubectl create ns argo
kubectl apply -n argo -f https://github.com/argoproj/argo-workflows/releases/download/v3.5.10/quick-start-minimal.yaml
## For config.yaml - output_dir and config.yaml should be set to the mounted directory in the minikube mount command and available
argo submit -n argo --watch workflow/argo-auto-workflow.yaml \
	-f config/config.yaml \
	-p samples='[{"sampleId":"SAMPLE_1","bam":"/mnt/data/RMNISTHS_0_point_1xdownsample.bam"},
	{"sampleId":"SAMPLE_2","bam":"/mnt/data/RMNISTHS_0_point_1xdownsample.bam"}]'
