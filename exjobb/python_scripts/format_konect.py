import scipy.io as sio
import numpy as np
import os
import tarfile

output_path = '../'

def extract_files():
    for files in dir_list:
        with tarfile.open(os.path.join(path, files)) as file:
            for f in file:
                if 'out.' in f.name:
                    file.extract(f.name, './Data', filter = 'data')

def ToObjectArray(my_array):
    object_array = np.zeros((len(my_array),), dtype=object)
    for i in range(len(my_array)):
        object_array[i] = my_array[i]
    return object_array

path = './konect_RawData'
dir_list = os.listdir(path)

extract_files()

path = './Data'
dir_list = os.listdir(path)

tags_dict = {
    "Animal network": ["dolphins", "moreno_mac"],
    "Computer communication network": ["p2p-Gnutella08", "p2p-Gnutella09", "p2p-Gnutella06", "p2p-Gnutella05"],
    "Interaction network": ["dimacs10-football", "moreno_crime"],
    "Online communitie": ["librec-filmtrust-trust", "ego-facebook", "petster-hamster-household", "petster-friendships-hamster", "petster-hamster-friend", "petster-hamster", "librec-ciaodvd-trust", "soc-sign-bitcoinalpha", "soc-sign-bitcoinotc"],
    "Electrical network": ["opsahl-powergrid"],
    "Traffic network": ["subelj_euroroad", "opsahl-openflights", "opsahl-usairport"],
    "Social network": ["iceland", "moreno_train", "moreno_highschool", "moreno_innovation", "arenas-jazz", "moreno_oz"],
             }

data = []
data_info = []

for files in dir_list:
    this_path = os.path.join(path, files)
    this_dir = os.listdir(this_path)[0]
    info = []
    with open(os.path.join(this_path, this_dir), "r") as file: # out.
        this_data = []
        for line in file:
            row = line.split()
            if row[0] == '%': # matlab comments
                info.append(line[2:])
            else:
                if len(row) < 3: # if unweighted, set all weights to 1
                    row.append(1)
                this_data.append(row)
        relevant_tag = ""
        file_name = (file.name).split("out.")[1]
        for tags in tags_dict:
            if file_name in tags_dict[tags]:
                relevant_tag = tags
        try:
            data.append(np.array(this_data, dtype=float))
            data_info.append({'name': file_name, 'description': info, 'tag': relevant_tag})
        except:
            print("Error in {}, could not form array".format(this_dir))

# reformat
formated_data = []
for graph in data:
    n = int(max(graph[:, 0].max(),graph[:, 1].max()))
    graph_matrix = np.zeros((n, n))
    for edge in graph:
        graph_matrix[int(edge[0])-1, int(edge[1])-1] = edge[2]
    formated_data.append(graph_matrix)

obj_arr = ToObjectArray(formated_data)
obj_arr_info = ToObjectArray(data_info)

# Note: graph saved with columns: in node, out node, weight
sio.savemat(os.path.join(output_path, 'konect.mat'), mdict={'data': obj_arr, 'info': obj_arr_info})