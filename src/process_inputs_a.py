import configparser


def parse_config():
    config = configparser.ConfigParser()
    config.read("config.ini")
    raw_transformation = config["PERFORM RAW TRANSFORMATION"].getboolean('transform')
    files = config["INPUT FILES"]
    output = config["OUTPUT FILES"]
    file_specifications = config['FILE SPECIFICATIONS']

    return raw_transformation, files, output, file_specifications
