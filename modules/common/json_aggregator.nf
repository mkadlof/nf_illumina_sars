process json_aggregator {
    tag "json_aggregator:${sampleId}"

    input:
    val pathogen
    val pipeline_version

     output:
     file('output.json')

    script:
    """
    echo ${pathogen}
    echo ${pipeline_version}
    python -c '
import json
import sys

output = {"output": {}}
output["pipeline_version"] = "${pipeline_version}"
output["pathogen"] = "${pathogen}"

with open("output.json", "w") as f:
    json.dump(output, f, indent=4)

# Here will go rest of aggregation
'

    """
}