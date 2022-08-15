## Export GCP MQL queries to CSV

This solution uses the Google Cloud [Monitoring Query Language (MQL)](https://cloud.google.com/monitoring/mql) to export queries to CSV.


### Before starting
1.  Ensure you have access to GCP Monitoring
2.  Have a project with active resources and metrics (i.e. usage within the last 25 hours). 


### Get Code and Configure MQL Queries
Clone the repository and enter the directory:

        git clone https://github.com/garg02/gcp-custom-metric-export.git
        cd gcp-custom-metric-export

Add and update the MQL queries as needed in the `config-template.py` file. You can find examples in the [MQL documentation](https://cloud.google.com/monitoring/mql/examples). This repository's default queries align CPU/GPU usage with an instance's batch number over time as created in [this project](https://github.com/garg02/gcp-custom-metric).


### Set up environment variables 
Export PROJECT_ID, START_TIME, and END_TIME variables to the shell environment. The start and end times are written in d'YYYY/MM/DD-HH:MM:SSZ' format, where Z is the time offset from UTC. Examples in Pacific time given below.

        export PROJECT_ID=`gcloud config get-value core/project`
        export START_TIME=d'2022/07/26-00:00:00-08:00'
        export END_TIME=d'2022/07/27-00:00:00-08:00'
        
### Create configuration

1.  Update the variables in the configuration file:
       
        envsubst < config-template.py > config.py
    
    If envsubst is not installed, manually open the config-template.py file, update the variables, and save it as config.py.

### Export metrics to csv
1.  Run the script

        python3 export_metric_data.py

