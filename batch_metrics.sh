#!/bin/bash
source /tmp/gce_meta
source /tmp/curr_time
batch_metrics -batch=$CURRTIME -instID=$INSTANCE_ID -projID=$PROJECT_ID -zoneID=$ZONE_ID
