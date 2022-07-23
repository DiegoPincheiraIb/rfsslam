"""
csv_list_generator.py
    Generates csv of timestamp's csv.
"""
import os
import csv
import yaml


def get_csv_names():
    """
    Reads csv files contained in folder, and stores them onto list.
    """
    with open(
            "src/python_script/config/cfg_parameters.yaml",
            "r", encoding='utf8') as yaml_obj:
        chosen_id = yaml.load(yaml_obj)["frames"]["chosen_id"]
    main_path = "data/rgbd/" + chosen_id + "/"
    list_csv = os.listdir(main_path + "csv_files/")
    list_csv = [[item_obj] for item_obj in list_csv]
    return list_csv, main_path, chosen_id


def save_list_onto_csv(list_obj, name_obj):
    """
    Saves list onto csv file.

    Parameters
    ---------
    list_obj: list
        A
    name_obj: str
        B
    """
    file = open(name_obj + ".csv", 'w+', newline='', encoding="utf-8")
    with file:
        write = csv.writer(file)
        write.writerows(list_obj)


def main():
    """
    Executes main code.
    """
    list_obj, main_path, chosen_id = get_csv_names()
    list_obj.sort()
    save_list_onto_csv(
        list_obj,
        main_path
        + "timestamp_list_" + chosen_id)
    print("Cantidad total: ", len(list_obj))
    return 0


if __name__ == "__main__":
    main()
