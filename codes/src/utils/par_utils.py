import os
class SimpleNamespace(object):
    '''
    Useful when defining a set of parameters
    '''
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def getListOfFiles(dirName):
    # create a list of file and sub directories 
    # names in the given directory 
    listOfFile = os.listdir(dirName)
    allFiles = list()
    # Iterate over all the entries
    for entry in listOfFile:
        # Create full path
        fullPath = os.path.join(dirName, entry)
        # If entry is a directory then get the list of files in this directory 
        if os.path.isdir(fullPath):
            allFiles = allFiles + getListOfFiles(fullPath)
        else:
            allFiles.append(fullPath)
                
    return allFiles


def read_last_line(file_name):
    with open(file_name, 'r') as f:
        last_line = f.readlines()[-1]
        
    last_line_list = last_line.split(',')

    return last_line_list