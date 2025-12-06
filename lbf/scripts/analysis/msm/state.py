import sys
import os
import fnmatch
import inspect
import pickle

from collections import defaultdict

class MacrostateMap:

    '''
    Mapping from a state description, i.e. 3 for trimer, the canonical graph representation, 
    or a tuple of coordinates (3,2), into a unique index to describe the state. The indices are assigned
    in the order the states are seen in trajectory data analysis. 

    Uses dicts to convert between states and indices. Unknown states can be added to the map as they are found
    '''

    def __init__(self, verbose=False):

        self.__verbose = verbose

        #init dicts for the mappings
        self.__toIndex = defaultdict(lambda: None)
        self.__toState = defaultdict(lambda: None)

        return

    def update_maps(self, state, verbose = False):
        #check if the supplied state is in the map. If not add it

        if state not in self.__toIndex:

            L = self.get_num_states()

            self.__toIndex[state] = L
            self.__toState[L] = state

            if self.__verbose or verbose:
                print("{} added to map with index {}".format(state, L))

            self.__been_updated = True
        
        return

    def state_to_index(self, state):

        if state in self.__toIndex:
            return self.__toIndex[state]
        
        return -1

    def index_to_state(self, index):

        if index < self.get_num_states():
            return self.__toState[index]
        
        return None

    def get_num_states(self):

        return len(self.__toIndex)
    
    def remove_entry(self, index):
        '''
        Remove the state in the specified index. All greater indices must be shifted
        down by 1 to account for the removal. 

        DANGER: if removing multiple states by index at once, you MUST go from largest
        index to smallest index. Otherwise, the indexing will not be consistent and 
        unintended states will be removed. 
        '''

        # Get the corresponding state and total number of states
        state = self.index_to_state(index)
        num_states = self.get_num_states()

        # Remove the state and index from each dict
        removed_state = self.__toState.pop(index)
        removed_index = self.__toIndex.pop(state, None)  # Use pop with a default value

        # Shift the indices for all states with greater index
        for shift_index in range(index + 1, num_states):

            # Get the corresponding state
            shift_state = self.index_to_state(shift_index)

            # Remove the original key-value pairs
            removed_shift_state = self.__toState.pop(shift_index)
            removed_shift_index = self.__toIndex.pop(shift_state, None)  # Use pop with a default value

            # Update the dictionaries with shifted values
            self.__toState[shift_index - 1] = removed_shift_state
            if removed_shift_index is not None:
                self.__toIndex[removed_shift_state] = shift_index - 1

        return

    def filter_by_size(self, sizes, verbose = False):
        #return all indices for states with size in sizes

        #if only an int is given, put in list so code is general
        if isinstance(sizes, int):
            sizes = [sizes]

        #init output list for indices
        output_indices = []

        #perform dictionary comprehension for each size in sizes
        for size in sizes:

            indices = [self.__toIndex[state] for state in self.__toIndex.keys() if state.get_size() == size ]
            output_indices += indices
            
            if verbose:
                print("Listing all states of size {}...".format(size))
                for index in indices:
                    print("Index {}, {}".format(index, self.__toState[index]))
                print()

        return output_indices






