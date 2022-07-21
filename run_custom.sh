curl -sSO https://dl.google.com/cloudagents/add-google-cloud-ops-agent-repo.sh
sudo bash add-google-cloud-ops-agent-repo.sh --also-install
sudo apt-get install -y stress-ng git golang-go
git clone https://github.com/garg02/gcp_report_batch.git
cd gcp_report_batch
go mod init metrics
go build batch_metrics.go
chmod +x *.sh

echo "export CURRTIME=$(date +%s)" > /tmp/curr_time
echo "export PROJECT_ID=$(gcloud info --format='value(config.project)')" > /tmp/gce_meta
echo "export INSTANCE_ID=$(curl http://metadata.google.internal/computeMetadata/v1/instance/id -H Metadata-Flavor:Google)" >> /tmp/gce_meta
echo "export ZONE_ID=$(curl http://metadata.google.internal/computeMetadata/v1/instance/zone -H Metadata-Flavor:Google | rev | cut -d/ -f1 | rev)" >> /tmp/gce_meta
echo -e "*/30 * * * * bash -c ./stress.sh >> stress_stdout.log 2>> stress_stderr.log\n*/1 * * * * bash -c ./batch_metrics.sh >> batch_metrics_stdout.log 2>> batch_metrics_stderr.log" | crontab -u $USER -
