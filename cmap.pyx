cpdef list decomposition(list data):
    cdef list vector = []
    cdef list out_genes = []
    cdef list out_vector = []
    cdef result = []
    cdef int i, j
    for j in range(len(data)):
        out_genes = []
        out_vector = []
        for i in range(len(data[j][1])):

            if data[j][1][i] in data[j][0]:
                out_genes.append(data[j][1][i])
                out_vector.append(1.0)
            else:
                out_genes.append(data[j][1][i])
                out_vector.append(0.0)
        vector.append(out_genes)
        vector.append(out_vector)
    #print(vector)
    return vector