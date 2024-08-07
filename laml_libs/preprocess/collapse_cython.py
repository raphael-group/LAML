#cimport cython
import numpy as np
OFFSET = 33

def hamming_distance(first,second):
    #cdef int i
    d = 0
    length = len(first)
    
    for i in range(length):
        if first[i] != second[i]:
            d += 1
            
    return d

#@cython.boundscheck(False)
def hamming_distance_matrix(seqs):
    #cdef int i, j, k, d, n, seq_length

    ints = np.array([list(s.encode()) for s in seqs])
    ints_view = ints
    n, seq_length = ints.shape

    ds = np.zeros((n, n), int)
    ds_view = ds

    for i in range(n):
        for j in range(i + 1, n):
            d = 0
            for k in range(seq_length):
                if ints_view[i, k] != ints_view[j, k]:
                    d += 1

            ds_view[i, j] = d

    return ds

#@cython.boundscheck(False)
def register_corrections(ds, max_UMI_distance, UMIs):
    n = len(ds)
    corrections = {}

    # Moving from least common to most common, register a correction
    # from a UMI to the most common UMI that is within Hamming distance
    # max_UMI_distance of it.
    for j in range(n - 1, -1, -1):
        for i in range(j - 1, -1, -1):
            if ds[i, j] <= max_UMI_distance:
                corrections[UMIs[j]] = UMIs[i]
    
    # If a correction points to a UMI that is itself going to be corrected,
    # propogate this correction through.  
    for from_, to in list(corrections.items()):
        while to in corrections:
            to = corrections[to]

        corrections[from_] = to
    
    return corrections

def hq_hamming_distance(first_seq, second_seq, first_qual, second_qual, min_q):
    d = 0
    length = len(first_seq)
    floor = min_q + OFFSET
    
    for i in range(length):
        if (first_seq[i] != second_seq[i]) and (first_qual[i] >= floor) and (second_qual[i] >= floor):
            d += 1
            
    return d

def hq_mismatches_from_seed(seed, seq, qual, min_q):
    d = 0
    length = len(seq)
    floor = min_q
    
    for i in range(length):
        if (seq[i] != seed[i]) and (qual[i] >= floor):
            d += 1
            
    return d
