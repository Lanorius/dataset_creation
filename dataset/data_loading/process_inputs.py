import configparser
import sys


def parse_config():
    config = configparser.ConfigParser()
    config.read("config.ini")
    files = config["INPUT FILES"]

    return(files)
