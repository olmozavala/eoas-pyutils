import os
from os import walk, listdir
from os.path import join
def create_folder(output_folder):
    """
    It creates a folder only if it doesn't exist
    Args:
        output_folder:
    """
    if not(os.path.exists(output_folder)):
        os.makedirs(output_folder)

def all_files_in_folder(input_folder, file_ext=None):
    """
    Gets all the filenames in a folder
    :param input_folder:
    :param file_ext:
    :return:
    """
    paths = []
    file_names = []
    for root,d_names,f_names in os.walk(input_folder):
        for f in f_names:
            if file_ext is None or f.find(file_ext) != -1:
                paths.append(os.path.join(root, f))
                file_names.append(f)

    return file_names, paths


def str2bool(cstr: str) -> bool:
    """
    It compares a string with anything like true, and it returns True or False
    Args:
        cstr:

    Returns:
        Boolean value of the string
    """
    return cstr in ['True', 'true', 't', True]

if __name__ == '__main__':
    print(all_files_in_folder("/home/olmozavala/Dropbox/MyProjects/OZ_LIB/eoas-ai-template/test_data/GOMb0.04"))