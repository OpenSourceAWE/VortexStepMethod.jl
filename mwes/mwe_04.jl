# read the vsm_settings.yaml file

import YAML

filename = "vsm_settings.yaml"
data = YAML.load_file(joinpath("data", filename))

