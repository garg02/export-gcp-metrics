#!/bin/bash
source /tmp/gce_meta
source /tmp/curr_time
gcp_report_batch/batch_metrics -batch=$CURRTIME -instID=$INSTANCE_ID -projID=$PROJECT_ID -zoneID=$ZONE_ID
