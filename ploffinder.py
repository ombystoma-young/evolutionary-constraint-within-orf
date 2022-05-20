import csv
from collections import OrderedDict


def observed_ptv(variants_file, coords_list):
    """
    processes data on observed protein truncation variants
    :param variants_file: (str) the name of the file containing the variants
    :param coords_list: (list) with window coordinates
    :return: two (lists): total allele counts per window and mean allele numbers per windows
    (if n = 0, returns total mean for this window)
    """
    coord_dict = OrderedDict()
    if coords_list[0][0] > coords_list[0][1]:
        for coords in coords_list:
            coord_dict[(coords[1], coords[0])] = [0, 0, 0]
    else:
        for coords in coords_list:
            coord_dict[(coords[0], coords[1])] = [0, 0, 0]

    variants_dict = {}
    total_an = 0
    with open(variants_file, 'r') as in_f:
        variants_table = csv.reader(in_f, delimiter='\t')
        for row in variants_table:
            variants_dict[row[1]] = [int(row[2]), int(row[3])]
            total_an += int(row[3])

    for variant in variants_dict:
        for coords in coord_dict:
            if int(coords[0]) <= int(variant) <= int(coords[1]):
                coord_dict[coords][0] += variants_dict[variant][0]
                coord_dict[coords][1] += variants_dict[variant][1]
                coord_dict[coords][2] += 1
    an_list = []
    ac_list = []

    for coords in coord_dict:
        ac_list.append(coord_dict[coords][0])
        if coord_dict[coords][2] != 0:
            an_list.append(coord_dict[coords][1] / coord_dict[coords][2])
        else:
            an_list.append(total_an / len(variants_dict))
    return ac_list, an_list
