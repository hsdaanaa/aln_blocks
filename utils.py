#! usr/bin/env Python3
# Functions to extract file paths, list file names and to create new folders
#------------------------------------------------------
# Dependencies
import sys, os, datetime
#------------------------------------------------------
def get_f_pathlist_from_folder(folder_path, ext = None):
    """takes an input path to a folder and returns a list containing
    paths to files in the folder. if <ext> is None, all paths are returned
    regardless of the file extension. if ext is e.g. ['.txt'] the function
    will return paths for all file names with .txt as the extension.
    
    note that only paths in the top-level directory are returned.
    
    parameters
    ----------
    folder_path: str
       path to folder to fetch paths to files
    ext: None or list
       specifies whether to fetch files with a particular extension
  
    returns
    -------
    list"""
    
    list_of_files = []
        
    try:
        assert os.path.isdir(folder_path) == True
    except:
        print('\tError in get_f_pathlist_from_folder')
        print('\t\tpath to folder was invalid')
        return -1 
    
    try: 
        assert ext == None or isinstance(ext, list)
    except:
        print('\tError in get_f_pathlist_from_folder')
        print('\t\tinput for <ext> was invalid. can only be None or list type')
        return -1 
    try:          
        folder_path  = os.path.abspath(folder_path)
        folder_items = os.listdir(folder_path)

        # loop over folder items
        for item in folder_items:
            # get item path
            path_to_item = os.path.join(folder_path, item) 
            # check path validity
            if os.path.isfile(path_to_item) == True: 
                # check ext setting and append path to item as necessary 
                if ext != None:   
                    if type(ext) == list: 
                        for extension in ext:
                            if item.endswith(extension):
                                list_of_files.append(path_to_item)
                    else:
                        item.endswith(ext)
                        list_of_files.append(path_to_item)
                else:
                    list_of_files.append(path_to_item)
    except Exception as error:
        print('\tError in get_f_pathlist_from_folder')
        print('\t\tgot: {}'.format(error))
        return -1 
    return list_of_files
#-------------------
def make_new_folder(path_to_folder, folder_name):
    """creates a a new folder
    
    parameters
    ----------
    path_to_folder: str
        path to the folder where a new folder should be created
    folder_name: str
        name of the new folder
    
    returns 
    -------
    str"""
    try:
        assert isinstance(path_to_folder, str)
    except:
        print('\tError in make_new_folder')
        print('\t\t<path_to_folder> was not a string. got: {}'.format(type(path_to_folder)))
        return -1
    try:
        assert os.path.isdir(path_to_folder) == True
    except:
        print('\tError in make_new_folder')
        print('\t\tinput folder path was an invalid. got: {}'.format(path_to_folder))
        return -1
    
    # new folder name
    new_folder_path = os.path.join(path_to_folder, folder_name)
    # check if new_folder_path already exists
    try:
         assert os.path.isdir(new_folder_path) == False
    except:
        print('\tError in make_new_folder')
        print('\t\tCannot create new folder. A folder with the same name exists: {}'.format(folder_name))
        return -1

    # make folder
    try:
        os.makedirs(new_folder_path, exist_ok = False)
    except Exception as error: 
        print('\tError in make_new_folder')
        print('\t\tGot: {}'.format(error))
        return -1
        
    return new_folder_path
#---------------------------------------------------------------
def write_flist_for_folder(path_to_folder): 
    """creates a file containing the number of files in a folder
    and the names of the files. the output file is named 00_filelist.txt.
    This file name is excluded from file counts.
    
    parameters
    ----------
    path_to_folder: str
        path to the folder to fetch file names
    
    returns
    -------
    None"""
    
    try:
        assert isinstance(path_to_folder, str)
    except:
        print('\tError in write_flist_for_folder')
        print('\t\t<path_to_folder> was not a string. got: {}'.format(type(path_to_folder)))
        return -1
                    
    # extract file paths
    try:
        file_path_list = get_f_pathlist_from_folder(path_to_folder)
    except Exception as error: 
        print('\tError in write_flist_for_folder')
        print('\t\tgot: {}'.format(error))
    
    try:
        # retain file names
        file_name_list = [os.path.basename(i) for i in file_path_list]
        file_name_list = [i for i in file_name_list if not i == '00_filelist.txt']

        # create a file called 00_filelist.txt and write file names
        with open(path_to_folder + '/00_filelist.txt', 'w') as file:
            file.write('# files: {}\n'.format(len(file_name_list)))
            for file_name in file_name_list:
                file.write(file_name  + '\n')
    except Exception as error:
        print('\tError in write_flist_for_folder')
        print('\t\tgot: {}'.format(error))
#------------------------------------------------------
def get_current_date():
    """returns current date YMD (e.g. 220101 means 1st Jan 2022)"""
    current_date  = datetime.datetime.now()
    current_year  = current_date.year
    current_month = '0{}'.format(current_date.month) if current_date.month < 10 else current_date.month
    current_day   = '0{}'.format(current_date.day) if current_date.day < 10 else current_date.day

    return '{}{}{}'.format(current_year, current_month,current_day)[2:]
#------------------------------------------------------
def map_dict_vals(input_dict, func_name): 
    new_dict = {}
    
    for key in input_dict:
        value         = input_dict[key]
        new_value     = func_name(value)
        new_dict[key] = new_value
        
    return new_dict

import os, shutil, sys
#------------------------------------------------------
def get_all_paths_in_dir_and_subdirs(dir_path, verbose = 0):
    """takes an input path to a directory, and returns a list
    of all pathsin the directory
    
    strategy
    --------
    -extracts file/folder names from input dir into a list
    -loops list and checks if subdirs exist. If true,
    does recursion, else, appends the aboslute file path to 
    a global variable.the latter is the base case of the recursion.
    *the recursion call allows extraction of files from subdirs,
    subsubdirs...etc
        
    notes
    -----
    -uses:
       os.listdir to extract the names of all files and folders in dir.
       os.path.isdir to check if file/folder is a dir
       os.path.abspath to extract full path of file.
    -handles invalid input type and also when dir is not found.
    
    parameters
    ----------
    dir_path: str
       path to dir
    verbose: 0,1
       set to 1 for debugging 
       
    return
    ------
     list"""
    
    # check_inputs(str, 'your input was not a string', dir_path)
    all_paths_in_dirs_and_subdirs = []
    try:
        # extract names of files and folders in dir into a list
        list_of_items_in_folder = os.listdir(dir_path)

        # loops list of names of files and folders, and appends filesnames to a global variable. 
        # does recursion for dirs to extract its file names
        for item in list_of_items_in_folder:
            
            full_item_path = os.path.join(dir_path,item)
            
            if (os.path.isdir(full_item_path) == True): # checks if item is a dir
                # subdir_path = os.path.join(os.path.abspath(dir_path),item) # extracts full path of item(subdir) to use for recursion
                all_paths_in_dirs_and_subdirs += get_all_paths_in_dir_and_subdirs(full_item_path)
                all_paths_in_dirs_and_subdirs.append(full_item_path)

            else:
                # this is the base case of the recursion
                # file = os.path.join(dir_path,item)
                all_paths_in_dirs_and_subdirs.append(full_item_path) # appends file name to global variable

        return all_paths_in_dirs_and_subdirs
    
    except FileNotFoundError:    
        print('your input path does not exist: {}'.format(dir_path))
#---------------------------------------------------------------------
def conv_nested_list_to_list(nlist_var, verbose = 0, none_if_error = True): 
    """takes an input list of lists (nested list). and 
    extracts all 'inner' list elements into the main list,
    making a 1 dimensional list.
    
    strategy
    --------
    -loops list and appends list items to a new list var.
    when the item is a list(inner list), the inner list is
    also looped and its items are appended to the new list. 
    
    parameters
    ----------
    nlist_var: list of list
        nested list to process
    verbose: 0,1
        optional parameter for debugging

    returns
    -------
    list"""

    # prints error message for invalid input
    try:
        assert isinstance(nlist_var, list)
    except:
        if none_if_error == True:
            return None
        else:
            print('your input was not a list')
    #check_input_type(list,'your input was not a list', nlist_var)
    
    new_list = [] # this allows appending items of list and inner list 
    
    
    # loops list and checks if each item is a list
    # if true, loops inner list and appends its elements to the new list
    # else, appends the items to the new list
    for item in nlist_var: 
        if isinstance(item, list) == True:
            items = conv_nested_list_to_list(item)
            new_list.extend(items)
        else:
            new_list.append(item)
            
        if verbose == 1: 
            print("current list item: {0}".format(item))
            
    if verbose == 1: 
        print('created list: {}'.format(new_list))

    return new_list