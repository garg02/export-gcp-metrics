#!/bin/bash
source /tmp/gce_meta
source /tmp/curr_time_gpu
gcp_report_batch/batch_metrics_gpu -batch_gpu=$CURRTIME -instID=$INSTANCE_ID -projID=$PROJECT_ID -zoneID=$ZONE_ID
