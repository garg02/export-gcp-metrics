import logging
import config
import requests
import subprocess
import pandas as pd
import os

token = None

def get_access_token_from_gcloud(force=False):
    global token
    if token is None or force:
        token = subprocess.check_output(
            ["gcloud", "auth", "application-default", "print-access-token"],
            text=True,
        ).rstrip()
    return token


def get_mql_result(token, query):
    q = f'{{"query": "{query}"}}'

    headers = {"Content-Type": "application/json",
               "Authorization": f"Bearer {token}"}
    return requests.post(config.QUERY_URL, data=q, headers=headers).json()


def build_rows(data):
    logging.debug("build_row")
    dataframe = {}

    labelDescriptors = data["timeSeriesDescriptor"]["labelDescriptors"]
    pointDescriptors = data["timeSeriesDescriptor"]["pointDescriptors"] 
    
    label_descriptions = []
    for i in range(len(labelDescriptors)):
        for label_description in labelDescriptors[i].values():
            label_descriptions.append(label_description)
            dataframe[label_description] = []
    dataframe["start_time"] = []
    dataframe["end_time"] = []
    dataframe["value"] = []

    for timeseries in data["timeSeriesData"]:
        labelValues = timeseries["labelValues"]
        pointData = timeseries["pointData"]
        numPoints = len(pointData)

        for i, label_description in zip(range(len(labelValues)), label_descriptions):
            for label_value in labelValues[i].values():
                dataframe[label_description] += [label_value] * numPoints

        # handle >= 1 points, potentially > 1 returned from Monitoring API call
        for point_idx in range(numPoints):
            dataframe["start_time"].append(pointData[point_idx]["timeInterval"]["startTime"])
            dataframe["end_time"].append(pointData[point_idx]["timeInterval"]["endTime"])
            
            # map the API value types to the BigQuery value types
            value_type = pointDescriptors[0]["valueType"] 
            api_value_type_index = config.API_VALUE_MAP[value_type]
            value = timeseries["pointData"][point_idx]["values"][0][api_value_type_index]

            if value_type == "DISTRIBUTION":
                return # error as this script does not handle Distribution data
            else:
                dataframe["value"].append(value)

    length = len(dataframe["value"])
    for j in range(len(pointDescriptors)):
        for k, v in pointDescriptors[j].items():
            dataframe[k] = [v]*length

    return pd.DataFrame(dataframe)

def save_to_csv(token):
    for metric, query in config.MQL_QUERYS.items():
        metric = metric.replace("/", "-")
        result = get_mql_result(token, query)
        if result.get("timeSeriesDescriptor"):
            df = build_rows(result)
            df.to_csv(metric+".csv", index=False)
            os.makedirs("stats", exist_ok=True)
            df.groupby(['metric.batch_num'])['value'].agg('sum').to_csv("stats"/metric+"_batched_aggr.csv")
            print(metric + ".csv written to file")
        else:
            print("No data retrieved")

if __name__ == "__main__":
    cmd = ["gcloud", "auth", "application-default", "login"]
    subprocess.Popen(cmd).wait()
    token = get_access_token_from_gcloud()
    save_to_csv(token)
