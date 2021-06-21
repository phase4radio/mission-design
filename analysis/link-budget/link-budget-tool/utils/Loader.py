import yaml
import os

# Allows include statements inside YAML files
#
# Taken from:
#   https://stackoverflow.com/questions/528281/how-can-i-include-a-yaml-file-inside-another
class Loader(yaml.SafeLoader):

    def __init__(self, stream):

        self._root = os.path.split(stream.name)[0]

        super(Loader, self).__init__(stream)

    def include(self, node):

        filename = os.path.join(self._root, self.construct_scalar(node))

        with open(filename, 'r') as f:
            return yaml.load(f, Loader)

Loader.add_constructor('!include', Loader.include)


# read in the YAML configuration files supplied in the filenames
def load_configurations(satellite_configuration_file, groundstation_configuration_file, channel_configuration_file):

    satellite_configuration = load_configuration(satellite_configuration_file)
    groundstation_configuration = load_configuration(groundstation_configuration_file)
    channel_configuration = load_configuration(channel_configuration_file)

    return satellite_configuration, groundstation_configuration, channel_configuration


# read in the YAML configuration file
def load_configuration(file):

    with open(file) as file:
        configuration = yaml.load(file, Loader)

    return configuration