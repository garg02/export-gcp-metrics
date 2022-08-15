# Replace the values as needed
PROJECT_ID = "$PROJECT_ID"

# Add/Update the queries for your metrics
MQL_QUERYS = {
"minimap-CPU-Usage-by-Batch":
"""
fetch gce_instance
| { metric 'agent.googleapis.com/processes/cpu_time';
      metric 'custom.googleapis.com/batch_num'
    }
| join
| filter (metric.command == 'minimap2')
| value val(0)
| align rate(1m)
| every 1m
| group_by
    [metric.batch_num, metadata.system.name: metadata.system_labels.name],
    [value_cpu_time_aggregate: aggregate(value.cpu_time)]
| within $START_TIME, $END_TIME
""",
"minimap-Resident-Memory-by-Batch":
"""
fetch gce_instance
| 	{
	metric 'agent.googleapis.com/processes/rss_usage';
	metric 'custom.googleapis.com/batch_num'
	} 
| join 
| filter (metric.command == 'minimap2')
| value val(0)
| group_by 1m, [value_rss_usage_mean: mean(value.rss_usage)]
| every 1m
| group_by [metric.batch_num, metadata.system.name: metadata.system_labels.name],
    [value_rss_usage_mean_aggregate: aggregate(value_rss_usage_mean)]
| within $START_TIME, $END_TIME
""",
"minimap-Virtual-Memory-by-Batch":
"""
fetch gce_instance
| 	{
	metric 'agent.googleapis.com/processes/vm_usage';
	metric 'custom.googleapis.com/batch_num'
	} 
| join 
| filter (metric.command == 'minimap2')
| value val(0)
| group_by 1m, [value_vm_usage_mean: mean(vm_usage)]
| every 1m
| group_by [metric.batch_num, metadata.system.name: metadata.system_labels.name],
    [value_vm_usage_mean_aggregate: aggregate(value_vm_usage_mean)]
| within $START_TIME, $END_TIME
""",
"minimap-Bytes-Read-by-Batch":
"""
fetch gce_instance
| 	{
	metric 'agent.googleapis.com/processes/disk/read_bytes_count';
	metric 'custom.googleapis.com/batch_num'
	} 
| join 
| filter (metric.command == 'minimap2')
| value val(0)
| group_by 1m, [value_read_bytes_mean: mean(read_bytes_count)]
| every 1m
| group_by [metric.batch_num, metadata.system.name: metadata.system_labels.name],
    [value_read_bytes_mean_aggregate: aggregate(value_read_bytes_mean)]
| within $START_TIME, $END_TIME
""",
"minimap-Bytes-Written-by-Batch":
"""
fetch gce_instance
| 	{
	metric 'agent.googleapis.com/processes/disk/write_bytes_count';
	metric 'custom.googleapis.com/batch_num'
	} 
| join 
| filter (metric.command == 'minimap2')
| value val(0)
| group_by 1m, [value_write_bytes_mean: mean(write_bytes_count)]
| every 1m
| group_by [metric.batch_num, metadata.system.name: metadata.system_labels.name],
    [value_write_bytes_mean_aggregate: aggregate(value_write_bytes_mean)]
| within $START_TIME, $END_TIME
""",
"samtools-CPU-Usage-by-Batch":
"""
fetch gce_instance
| { metric 'agent.googleapis.com/processes/cpu_time';
      metric 'custom.googleapis.com/batch_num'
    }
| join
| filter (metric.command == 'samtools')
| value val(0)
| align rate(1m)
| every 1m
| group_by
    [metric.batch_num, metadata.system.name: metadata.system_labels.name],
    [value_cpu_time_aggregate: aggregate(value.cpu_time)]
| within $START_TIME, $END_TIME
""",
"samtools-Resident-Memory-by-Batch":
"""
fetch gce_instance
| 	{
	metric 'agent.googleapis.com/processes/rss_usage';
	metric 'custom.googleapis.com/batch_num'
	} 
| join 
| filter (metric.command == 'samtools')
| value val(0)
| group_by 1m, [value_rss_usage_mean: mean(value.rss_usage)]
| every 1m
| group_by [metric.batch_num, metadata.system.name: metadata.system_labels.name],
    [value_rss_usage_mean_aggregate: aggregate(value_rss_usage_mean)]
| within $START_TIME, $END_TIME
""",
"samtools-Virtual-Memory-by-Batch":
"""
fetch gce_instance
| 	{
	metric 'agent.googleapis.com/processes/vm_usage';
	metric 'custom.googleapis.com/batch_num'
	} 
| join 
| filter (metric.command == 'samtools')
| value val(0)
| group_by 1m, [value_vm_usage_mean: mean(vm_usage)]
| every 1m
| group_by [metric.batch_num, metadata.system.name: metadata.system_labels.name],
    [value_vm_usage_mean_aggregate: aggregate(value_vm_usage_mean)]
| within $START_TIME, $END_TIME
""",
"samtools-Bytes-Read-by-Batch":
"""
fetch gce_instance
| 	{
	metric 'agent.googleapis.com/processes/disk/read_bytes_count';
	metric 'custom.googleapis.com/batch_num'
	} 
| join 
| filter (metric.command == 'samtools')
| value val(0)
| group_by 1m, [value_read_bytes_mean: mean(read_bytes_count)]
| every 1m
| group_by [metric.batch_num, metadata.system.name: metadata.system_labels.name],
    [value_read_bytes_mean_aggregate: aggregate(value_read_bytes_mean)]
| within $START_TIME, $END_TIME
""",
"samtools-Bytes-Written-by-Batch":
"""
fetch gce_instance
| 	{
	metric 'agent.googleapis.com/processes/disk/write_bytes_count';
	metric 'custom.googleapis.com/batch_num'
	} 
| join 
| filter (metric.command == 'samtools')
| value val(0)
| group_by 1m, [value_write_bytes_mean: mean(write_bytes_count)]
| every 1m
| group_by [metric.batch_num, metadata.system.name: metadata.system_labels.name],
    [value_write_bytes_mean_aggregate: aggregate(value_write_bytes_mean)]
| within $START_TIME, $END_TIME
"""
}

BASE_URL = "https://monitoring.googleapis.com/v3/projects"
QUERY_URL = f"{BASE_URL}/{PROJECT_ID}/timeSeries:query"

API_VALUE_MAP = {
    "INT64": "int64Value",
    "BOOL": "booleanValue",
    "DOUBLE": "doubleValue",
    "STRING": "stringValue",
    "DISTRIBUTION": "distributionValue"
}