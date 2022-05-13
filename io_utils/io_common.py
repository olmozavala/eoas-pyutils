import os
from os import walk, listdir
from os.path import join

def create_folder(output_folder):
    """ It simply verifies if a folder already exists, if not it creates it"""
    if not(os.path.exists(output_folder)):
        os.makedirs(output_folder)


