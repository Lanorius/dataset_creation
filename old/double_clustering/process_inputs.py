import configparser


def parse_config():
    config = configparser.ConfigParser()
    config.read("config.ini")
    tasks_to_perform = [config["PERFORM SUB-STEPS"].getboolean('create_drug_cluster'),
                        config["PERFORM SUB-STEPS"].getboolean('create_target_cluster'),
                        config["PERFORM SUB-STEPS"].getboolean('update_drug_target_interactions'),
                        config["PERFORM SUB-STEPS"].getboolean('prepare_data_for_chemVAE'),
                        config["PERFORM SUB-STEPS"].getboolean('save_key_errors')]

    file = config["INPUT FILES"]
    output = config["OUTPUT FILES"]
    params = config['PARAMETERS']

    return tasks_to_perform, file, output, params
