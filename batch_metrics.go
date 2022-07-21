package main

import (
        "context"
        "fmt"
        "log"
        "time"
        "flag"
        "strconv"

        monitoring "cloud.google.com/go/monitoring/apiv3/v2"
        googlepb "github.com/golang/protobuf/ptypes/timestamp"
        metricpb "google.golang.org/genproto/googleapis/api/metric"
        monitoredrespb "google.golang.org/genproto/googleapis/api/monitoredres"
        monitoringpb "google.golang.org/genproto/googleapis/monitoring/v3"
)

func main() {
        time_ptr := flag.Int64("batch", 0, "int value for batch metric")
        project_ptr :=  flag.String("projID", "", "string value for the project ID of the instance")
        instance_ptr :=  flag.String("instID", "", "string value for the ID of the instance")
        zone_ptr :=  flag.String("zoneID", "", "string value for the zone ID of the instance")
        flag.Parse()

                if *time_ptr == 0 || len(*project_ptr) == 0 || len(*instance_ptr) == 0 || len(*zone_ptr) == 0 {
                        log.Fatalf("arguments may not be empty")
                }
        fmt.Println(*time_ptr)
        ctx := context.Background()

        // Creates a client.
        client, err := monitoring.NewMetricClient(ctx)
        if err != nil {
                log.Fatalf("Failed to create client: %v", err)
        }

        // Prepares an individual data point
        dataPoint := &monitoringpb.Point{
                Interval: &monitoringpb.TimeInterval{
                        EndTime: &googlepb.Timestamp{
                                Seconds: time.Now().Unix(),
                        },
                },
                Value: &monitoringpb.TypedValue{
                        Value: &monitoringpb.TypedValue_Int64Value{
                                Int64Value: 0,
                        },
                },
        }

        // Writes time series data.
        if err := client.CreateTimeSeries(ctx, &monitoringpb.CreateTimeSeriesRequest{
                Name: fmt.Sprintf("projects/%s", *project_ptr),
                TimeSeries: []*monitoringpb.TimeSeries{
                        {
                                Metric: &metricpb.Metric{
                                        Type: "custom.googleapis.com/batch_num",
                                        Labels: map[string]string{
                                                "batch_num": strconv.Itoa(int(*time_ptr)),
                                        },
                                },
                                Resource: &monitoredrespb.MonitoredResource{
                                        Type: "gce_instance",
                                        Labels: map[string]string{
                                                "instance_id": *instance_ptr,
                                                        "zone": *zone_ptr,
                                        },
                                },
                                Points: []*monitoringpb.Point{
                                        dataPoint,
                                },
                        },
                },
        }); err != nil {
                log.Fatalf("Failed to write time series data: %v", err)
        }

        // Closes the client and flushes the data to Stackdriver.
        if err := client.Close(); err != nil {
                log.Fatalf("Failed to close client: %v", err)
        }
}
