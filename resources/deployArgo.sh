minikube start
nohup minikube mount /home/psmith/swgs-absolutecn/data/:/mnt/data/ &
argo submit -n argo --watch workflow/argo-auto-workflow.yaml
