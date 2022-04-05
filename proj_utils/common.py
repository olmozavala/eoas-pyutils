import os


def create_folder(output_folder):
    """ It only creates a folder if it doesn't exist"""
    if not(os.path.exists(output_folder)):
        os.makedirs(output_folder)
