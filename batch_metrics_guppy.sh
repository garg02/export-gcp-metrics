#!/bin/bash
source /tmp/gce_meta
source /tmp/curr_time_guppy
gcp_report_batch/batch_metrics_guppy -batch=$CURRTIME -instID=$INSTANCE_ID -projID=$PROJECT_ID -zoneID=$ZONE_ID
