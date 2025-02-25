import graph_tool.all as gt
import graph_tool
import scipy.io as sio
import numpy as np
import os

output_path = '../'

def GetData(graph_name):
    g = gt.collection.ns[graph_name]
    return [graph_tool.spectral.adjacency(g), g.is_directed()]

def ToObjectArray(my_array):
    object_array = np.zeros((len(my_array),), dtype=object)
    for i in range(len(my_array)):
        object_array[i] = my_array[i]
    return object_array

tags_dict = {
    "Animal network": ["dolphins", "macaques", "fresh_webs/AkatoreA", "fresh_webs/AkatoreB", "fresh_webs/Berwick", "fresh_webs/Blackrock", "fresh_webs/Broad", "fresh_webs/Canton", "fresh_webs/Catlins", "fresh_webs/Coweeta1", "fresh_webs/Coweeta17", "fresh_webs/DempstersAu", "fresh_webs/DempstersSp", "fresh_webs/DempstersSu", "fresh_webs/German", "fresh_webs/Healy", "fresh_webs/Kyeburn", "fresh_webs/LilKyeburn", "fresh_webs/Martins", "fresh_webs/Narrowdale", "fresh_webs/NorthCol", "fresh_webs/Powder", "fresh_webs/Stony", "fresh_webs/SuttonAu", "fresh_webs/SuttonSp", "fresh_webs/SuttonSu", "fresh_webs/Troy", "fresh_webs/Venlaw"],
    "Biological network": ["psi", "interactome_pdz", "celegans_metabolic", "plant_pol_kato", "yeast_transcription"],
    "Brain network": ["cintestinalis", "celegansneural", "celegans_2019/male_chemical","celegans_2019/male_gap_junction","celegans_2019/hermaphrodite_chemical","celegans_2019/hermaphrodite_gap_junction","celegans_2019/male_chemical_synapse","celegans_2019/hermaphrodite_chemical_synapse","celegans_2019/male_chemical_corrected","celegans_2019/male_gap_junction_corrected","celegans_2019/hermaphrodite_chemical_corrected","celegans_2019/hermaphrodite_gap_junction_corrected"],
    "Computer communication network": ["email_company"],
    "Economic netowrk": ["fao_trade", "product_space/HS", "product_space/SITC"],
    "Ecological network": ["foodweb_baywet", "foodweb_little_rock", "messal_shale"],
    "Interaction network": ["reality_mining", "football_tsevans", "football", "dom/Strauss_2019d", "dom/Franz_2015a", "dom/Franz_2015e", "dom/Franz_2015d", "dom/Shimoji_2014c", "contact", "kidnappings", "blumenau_drug", "copenhagen/sms", "copenhagen/fb_friends", "copenhagen/calls", "copenhagen/bt"],
    "Online communitie": ["revolution"],
    "Electrical network": [], # power grid, etc.
    "Scientific computing": [],
    "Social network": ["moreno_sociograms/grade_8", "cs_department", "terrorists_911", "train_terrorists", "highschool", "law_firm", "baseball/player-player", "swingers", "student_cooperation", "jazz_collab", "residence_hall", "facebook_friends", "spanish_highschools/1", "spanish_highschools/2", "spanish_highschools/6", "spanish_highschools/11_1", "spanish_highschools/11_2", "spanish_highschools/11_4", "spanish_highschools/11_5", "spanish_highschools/11_6", "spanish_highschools/11_7", "spanish_highschools/11_9", "spanish_highschools/11_10"],
    "Traffic network": ["london_transport", "eu_airlines"], # roads, airports etc.
    "Technological network": ["wiki_science"],
}

data_info = []

# create list of all graphs to include everything in a single .mat file for matlab
data = []
count = 0
for tags in tags_dict:
    for graph_name in tags_dict[tags]:
        try:
            [d, i] = GetData(graph_name)
            if i: info = "directed"
            else: info = "undirected"
            data_info.append({'name': graph_name, 'description': info, 'tag': tags})
            data.append(d)
            count += 1
        except(KeyError): 
            print("KeyError at {}, graph was not added".format(graph_name))

# Note: graphs saved as adjacency matrixes (need to reformat in matlab)
sio.savemat(os.path.join(output_path, 'netzschleuder.mat'), mdict={'data': data, 'info': data_info})
print("Matlab file was created with {} graphs.".format(count))