#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##########################################
#     handle the file related  operation #
#          2017.12.06                    #
##########################################
__author__ = "K.R.Chow"
__version__ = "v1.1"

import sys, os
from functools import reduce


# accepted parameters for functions
def AcceptArgs(arg, *args):
    if arg not in args:
        args
        sys.stderr.write('Unkown parameters for "' + arg + '"\n')
        sys.exit()

# show argparse help when no aruguments
def ArgparseHelp(parser):
    if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

# return file status
def FileExist(file, error=None):
    if file is None:
        if error is not None:
            sys.stderr.write('The file is NoneType\n')
        sys.exit()
    else:
        if os.path.exists(file)==False:
            if error is not None:
                sys.stderr.write('The file', file, 'is not existed!\n')
            sys.exit()

# map list to dict structure
def List2Dict(listT, separate=2):
    return dict(zip(listT[::separate], listT[separate - 1::separate]))

# map list to str list
def List2Str(listT):
    return map(str, listT)

# map list to int list
def List2Int(listT):
    return map(int, listT)

# retrive the file name and remove sufix
def DeSufix(file):
    return os.path.splitext(os.path.basename(file))[0]

# retrive the intersected keys in list consisted of dicts
# listDict = [dict1, dict2, dict3 ...]
# return a set
def DictKeysSet(*listDict):
    newList = []
    for i in listDict:
        newList.append(i.keys())
    if len(newList) == 1:
        return set(newList[0])
    else:
        return reduce(lambda x,y: set(x) & set(y), newList)

# retrive the intersected elements in lists
# nestList = [list1, list2, list3 ...]
# return a set
def ListInter(*nestList):
    if len(nestList) == 1:
        return set(nestList[0])
    else:
        return reduce(lambda x,y: set(x) & set(y), nestList)

