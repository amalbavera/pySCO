#############################################################################################################################

##       #### ########  ########     ###    ########  #### ########  ######
##        ##  ##     ## ##     ##   ## ##   ##     ##  ##  ##       ##    ##
##        ##  ##     ## ##     ##  ##   ##  ##     ##  ##  ##       ##
##        ##  ########  ########  ##     ## ########   ##  ######    ######
##        ##  ##     ## ##   ##   ######### ##   ##    ##  ##             ##
##        ##  ##     ## ##    ##  ##     ## ##    ##   ##  ##       ##    ##
######## #### ########  ##     ## ##     ## ##     ## #### ########  ######

#############################################################################################################################

import os
import itertools

import multiprocessing as mp

from pathlib import Path as path
#
### local libraries
#
from ..read import read

#############################################################################################################################

##     ##  #######  ########  ##     ## ##       ########  ######
###   ### ##     ## ##     ## ##     ## ##       ##       ##    ##
#### #### ##     ## ##     ## ##     ## ##       ##       ##
## ### ## ##     ## ##     ## ##     ## ##       ######    ######
##     ## ##     ## ##     ## ##     ## ##       ##             ##
##     ## ##     ## ##     ## ##     ## ##       ##       ##    ##
##     ##  #######  ########   #######  ######## ########  ######

#############################################################################################################################
### Perform a read/write operation on a spin-state

def spin_state_set_io(state:str=None, parent_path:path=None, io_operation:any=None, option:str="read") -> any:

    if option == "read":
        state_data = io_operation(filesdir=parent_path/state)
    else:
        raise RuntimeError("unrecognized option '{option}'")

    return state_data

#############################################################################################################################
### Permutations between reference states

def permute_spin_staes(ls:str=None, hs:str=None, length:int=0, existsin:path=None, skip_reference=True) -> list:

    states = []

    start = 1 if skip_reference else 0
    stop  = 0 if skip_reference else 1

    for idx in range(start,length+stop):
        permutations = itertools.permutations(idx*ls + (length - idx)*hs)

        for state in permutations:
            spinstate = "".join(state)
            states.append(spinstate)

    states = list(set(states))

    if existsin is not None:
        remove_list = [state for state in states if not os.path.exists(existsin/state)]

        for state in remove_list:
            states.remove(state)
            print(f"{RuntimeWarning.__name__}: unable to locate {state}, skipping configuration")

    return states

#############################################################################################################################
### Calculate permutations between two spin-states

def generate_mixed_spin_states(ls:any=None, hs:any=None, ms:str=None) -> list:

    if not ms.lower() == "all": raise RuntimeError("unrecognized option '{ms}'")

    class_low_spin_state   = ls.__class__.__name__
    class_high_spin_state  = hs.__class__.__name__

    read_class_attribute   = getattr(read, class_low_spin_state)

    parent_low_spin_state  = path(ls.filesdir.parent)
    parent_high_spin_state = path(hs.filesdir.parent)

    name_low_spin_state    = ls.filesdir.name
    name_high_spin_state   = hs.filesdir.name

    list_low_spin_states   = list(set(name_low_spin_state))
    list_high_spin_states  = list(set(name_high_spin_state))

    low_spin_state         = list_low_spin_states[0]
    high_spin_state        = list_high_spin_states[0]

    length_low_spin_state  = len(name_low_spin_state)

    if not parent_low_spin_state == parent_high_spin_state: raise RuntimeError(f"parent directory for ls and hs differs")
    if not class_low_spin_state  == class_high_spin_state:  raise RuntimeError(f"cannot mix data from {class_low_spin_state} and {class_high_spin_state}")
    if not len(list_low_spin_states) == len(list_high_spin_states): raise RuntimeError(f"more than two spin-states detected: ls = {list_low_spin_states}, hs = {list_high_spin_states}")

    list_mixed_states = permute_spin_staes(ls=low_spin_state, hs=high_spin_state, length=length_low_spin_state, existsin=parent_low_spin_state)

    total_states      = len(list_mixed_states)

    cpu_avail         = mp.cpu_count()

    cpus  = total_states if cpu_avail > total_states else cpu_avail
    pool  = mp.Pool(cpus)

    mparg = [(i, parent_low_spin_state, read_class_attribute, "read") for i in list_mixed_states]

    mixed_states = pool.starmap_async(spin_state_set_io, mparg).get()

    return mixed_states

#############################################################################################################################
