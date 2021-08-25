import configparser
import sys


def parse_config():
    config = configparser.ConfigParser()
    config.read("config.ini")
    files = config["INPUT FILES"]
    file_specifications = config['FILE SPECIFICATIONS']
    params = config['PARAMETERS']


    return(files, file_specifications, params)
