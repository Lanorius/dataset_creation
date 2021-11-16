import configparser


def parse_config():
    config = configparser.ConfigParser()
    config.read("config.ini")
    tasks_to_perform = [config["PERFORM RAW TRANSFORMATION"].getboolean('transform'),
                        config["PERFORM SUB-STEPS"].getboolean('create_drug_cluster'),
                        config["PERFORM SUB-STEPS"].getboolean('create_target_cluster'),
                        config["PERFORM SUB-STEPS"].getboolean('update_drug_target_interactions'),
                        config["PERFORM SUB-STEPS"].getboolean('prepare_data_for_chemVAE'),
                        config["PERFORM SUB-STEPS"].getboolean('save_key_errors')]
    files = config["INPUT FILES"]
    file_specifications = config['FILE SPECIFICATIONS']
    output = config["OUTPUT FILES"]
    params = config['PARAMETERS']

    return tasks_to_perform, files, file_specifications, output, params
