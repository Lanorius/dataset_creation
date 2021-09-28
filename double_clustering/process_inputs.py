import configparser


def parse_config():
    config = configparser.ConfigParser()
    config.read("config.ini")
    tasks_to_perform = [config["PERFORM SUB-STEPS"].getboolean('create_drug_matrix'),
                        config["PERFORM SUB-STEPS"].getboolean('create_target_cluster'),
                        config["PERFORM SUB-STEPS"].getboolean('create_drug_cluster'),
                        config["PERFORM SUB-STEPS"].getboolean('update_drug_target_interactions')]

    files = config["INPUT FILES"]
    output = config["OUTPUT FILES"]
    params = config['PARAMETERS']

    return tasks_to_perform, files, output, params
