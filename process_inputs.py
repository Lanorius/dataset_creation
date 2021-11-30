import configparser


def parse_config():
    config = configparser.ConfigParser()
    config.read("config.ini")
    tasks_to_perform = [config["PERFORM STEPS"].getboolean('transform'),
                        config["PERFORM STEPS"].getboolean('create_drug_cluster'),
                        config["PERFORM STEPS"].getboolean('create_target_cluster'),
                        config["PERFORM STEPS"].getboolean('update_drug_target_interactions')]

    sub_tasks_to_perform = [config["PERFORM SUB-STEPS"].getboolean('transform_Kd_to_pKd'),
                            config["PERFORM SUB-STEPS"].getboolean('save_frequency_data')]

    files = config["INPUT FILES"]
    file_specifications = config['FILE SPECIFICATIONS']
    output = config["OUTPUT FILES"]
    params = config['PARAMETERS']

    return tasks_to_perform, sub_tasks_to_perform, files, file_specifications, output, params
