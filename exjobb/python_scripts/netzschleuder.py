import graph_tool.all as gt
import scipy.io as sio
import numpy as np
import os

output_path = '../'

def ToObjectArray(my_array):
    object_array = np.zeros((len(my_array),), dtype=object)
    for i in range(len(my_array)):
        object_array[i] = my_array[i]
    return object_array

g = gt.collection.ns["moreno_sociograms/grade_8"]

g_data = g.get_edges([g.edge_index])


# for s, t, i in g.iter_edges([g.edge_index]):
#    print(s, t, i)


# reformat
data = []
formated_data = []
for graph in data:
    n = int(max(graph[:, 0].max(),graph[:, 1].max()))
    graph_matrix = np.zeros((n, n))
    for edge in graph:
        graph_matrix[int(edge[0])-1, int(edge[1])-1] = edge[2]
    formated_data.append(graph_matrix)

obj_arr = ToObjectArray(formated_data)

sio.savemat(os.path.join(output_path, 'netzschleuder.mat'), mdict={'data': obj_arr})