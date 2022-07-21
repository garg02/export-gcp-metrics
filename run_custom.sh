curl -sSO https://dl.google.com/cloudagents/add-google-cloud-ops-agent-repo.sh
sudo bash add-google-cloud-ops-agent-repo.sh --also-install
sudo apt-get install -y stress-ng git golang-go
git clone https://github.com/garg02/gcp_report_batch.git
cd gcp_report_batch
go mod init metrics
go build batch_metrics.go
chmod +x *.sh
