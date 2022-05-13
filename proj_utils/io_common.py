import os

def create_folder(output_folder):
    """
    It creates a folder only if it doesn't exist
    Args:
        output_folder:
    """
    if not(os.path.exists(output_folder)):
        os.makedirs(output_folder)

def str2bool(cstr: str) -> bool:
    """
    It compares a string with anything like true, and it returns True or False
    Args:
        cstr:

    Returns:
        Boolean value of the string
    """
    return cstr in ['True', 'true', 't', True]

