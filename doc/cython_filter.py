#!/usr/bin/python

import sys

default_val_key = "="
func_key = " def "

def parse_types(line):
    """
    Removes types from cython function arguments.
    """
    # if this is present, check the line for a
    # function definition
    
    if func_key in line:
        return get_args(line)
    else:
        return line

def get_args(line):
    """
    Returns the comma-separated arguments 
    (includes types, and default values)
    """

    # make sure that open/close parenthesis appear in this line
    assert( "(" in line and ")" in line)
    
    # get everything bewteen the first "("
    # and last ")"
    open_idx = line.find("(")+1
    # search line backward and use negative of index for close_idx
    close_idx = -(line[::-1].find(")")+1)
    # slice out the arguments and split on comma
    prefix = line[:open_idx]
    args = line[open_idx:close_idx].split(',')
    suffix = line[close_idx:]

    for i in range(len(args)):
        if arg_has_type(args[i]):
            args[i] = remove_type(args[i])
            
    return prefix+", ".join(args)+suffix

def arg_has_type(arg):
    """
    Returns True if the argument type is specified.
    Returns False if not.
    """

    # break into elements of the argument by spaces
    elements = arg.split()

    # if there is only the varaible name, there's no type
    # early exit
    if len(elements) == 1:
        return False

    # if there are two elements, then there must be a type
    # (or this isn't a proplerly formulated argument
    if len(elements) == 2:
        return True

    # both of the following assume a space between the equals
    # sign and other arguments
    
    # if any two arguments in the list of arguments isn't
    # separated by an equals sign, then a type exists
    
    if default_val_key in elements and len(elements) > 3:
        return True

    if default_val_key not in elements and len(elements) > 2:
        return True

    # return false otherwise
    return False

def remove_type(arg):
    """
    Extracts the variable name from an argument
    known to contain a type.
    Also returns the variable default value if present.
    """
    # split on spaces
    elements = arg.split()
    
    # if there is an equals sign,
    # we want to preserve the default value
    if default_val_key in elements:
        idx = elements.index(default_val_key)
        # return the variable, separator, and value
        return " ".join(elements[idx-1:idx+2])
    else:
        #otherwise the last element should be the variable name
        return elements[-1]

def main():
    """
    Main function
    """
    # open the file (first argument)
    f = open(sys.argv[1],'r')
    
    # if not a cython file, no filter
    if sys.argv[1].split(".")[-1] != "pyx":
        for line in f.readlines():
            print line
    # if this is a cython file, parse
    # types from function arguments
    else:
        for line in f.readlines():
            print parse_types(line)

    #close the file
    f.close()

# run main function
main()
