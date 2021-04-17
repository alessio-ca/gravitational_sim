from typing import Dict, Any
import yaml


def load_parameters() -> Dict[str, Any]:
    with open("parameters.yml") as file:
        params = yaml.load(file, Loader=yaml.FullLoader)

    return params
