minikube start
nohup minikube mount /home/psmith/swgs-absolutecn/data/:/mnt/data/ &
argo submit -n argo --watch workflow/argo-auto-workflow.yaml \
	-f config/config.yaml \
	-p samples='[{"sampleId":"SAMPLE_1","bam":"/mnt/data/RMNISTHS_0_point_1xdownsample.bam"},{"sampleId":"SAMPLE_2","bam":"/mnt/data/RMNISTHS_0_point_1xdownsample.bam"}]'
