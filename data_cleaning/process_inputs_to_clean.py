import configparser


def parse_config():
    config = configparser.ConfigParser()
    config.read("config.ini")
    files = config["INPUT FILES"]
    output = config["OUTPUT FILE"]
    file_specifications = config['FILE SPECIFICATIONS']
    params = config['PARAMETERS']

    return files, output, file_specifications, params
