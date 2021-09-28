import configparser


def parse_config():
    config = configparser.ConfigParser()
    config.read("config.ini")
    files = config["INPUT FILES"]
    output = config["OUTPUT FILES"]
    params = config['PARAMETERS']

    return files, output, params
