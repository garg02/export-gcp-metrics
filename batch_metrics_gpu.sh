#!/bin/bash
source /tmp/gce_meta
gcp_report_batch/batch_metrics_gpu -batch_gpu=$(date +%s) -instID=$INSTANCE_ID -projID=$PROJECT_ID -zoneID=$ZONE_ID
