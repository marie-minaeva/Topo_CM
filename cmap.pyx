#from cython.parallel import *

cpdef list decomposition(list data):
    cdef list vector = []
    cdef list out_genes = []
    cdef list out_vector = []
    cdef result = []
    cdef int i, j
    for j in range(len(data)):
        out_genes = []
        out_vector_1 = []
        out_vector_2 = []
        for gene, val in data[j].items():
            print(val)
            print("Here 1")
            if len(val) == 2:
                print("Here 2")
                out_vector_1.append(1.0)
                out_vector_2.append(1.0)
            if len(val) == 1 and val[0] == 0:
                print("Here 3")
                out_vector_1.append(1.0)
                out_vector_2.append(0.0)
            if len(val) == 1 and val[0] == 1:
                print("Here 4")
                out_vector_1.append(0.0)
                out_vector_2.append(1.0)

        vector.append(out_vector_1)
        vector.append(out_vector_2)
    return vector