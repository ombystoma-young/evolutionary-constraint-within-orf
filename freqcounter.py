import re
from itertools import islice
import csv


def fasta_parsing(file_name: str):
    """
    :param file_name: path to FASTA file
    :return: tuple(gene_name (str), chrm (str):  chromosome, strand (str): "+" or "-",
    list_of_all_things (list):
    [num_of_exon (int): number of exon, coordinates (tuple): 1-based [start, end], sequence (str)]
    """

    fasta_sign = r'>'
    list_of_all_things = []
    with open(file_name) as protein_coding_sequence:
        first_line = protein_coding_sequence.read()
        fasta_name = first_line.strip().split(':')
        chrm, strand, gene_name = parse_name_string(fasta_name)
    with open(file_name) as protein_coding_sequence:
        for line in protein_coding_sequence:
            # > finding:
            if re.search(fasta_sign, line):
                fasta_name = line.strip().split(':')
                num_of_exon, coords = get_gene_coordinates(fasta_name)
                list_of_all_things.append([num_of_exon, coords])
            # processing sequence:
            else:
                list_of_all_things[-1].append(line.strip())

    return gene_name, chrm, strand, list_of_all_things


def parse_name_string(string_name: list):
    """
    takes list from fasta_parsing function and return gene name, chromosome, strand
    :param string_name: list contains parts of >string of FASTA file.
    :return: tuple(chrm (str):  chromosome, strand (str): "+" or "-", gene_name (str))
    """
    chrm = string_name[6]
    strand = string_name[4]
    gene_name = string_name[3].split('=')[1]
    return chrm, strand, gene_name


def get_gene_coordinates(string_name: list):
    """
    takes list from fasta_parsing function and return gene coordinates: #exon, position.
    :param string_name: (list) contains parts of >string of FASTA file.
    :return: (tuple): (num_of_exon (int): number of exon, coordinates (tuple): 1-based [start, end])
    """
    num_of_exon = string_name[2]
    coord_from_str = tuple(string_name[-1].split('-'))
    coords = (int(coord_from_str[0]) + 1, int(coord_from_str[1]))
    return num_of_exon, coords


def get_gene_name(file_name: str):
    """
    parse only gene name from input file
       :param file_name: path to FASTA file
       :return: gene_name (str)
       """
    with open(file_name) as protein_coding_sequence:
        first_line = protein_coding_sequence.read()
        gene_name_str = first_line.strip().split(':')[3]
        gene_name = gene_name_str.split('=')[1]
    return gene_name


def freq_table_to_dict(freq_table_file: str):
    """
    creates a dictionary of mutation frequency in a three-nucleotide context based on a file
    (with dropped methylated trinucleotides)
    :param freq_table_file: path to file containing mutation frequencies
    (1 col - trinucleotide context, 2 col - ref base,
    3 - alt base, 4 - methylation level, 5 - mut_frequency)

    :return: freq_dict: dict where keys are (contexts, alternative bases) (tuple)
    with values are mutation frequencies (float)
    """
    origin_freq_dict = {}
    with open(freq_table_file) as freq_table:
        table = csv.reader(freq_table, delimiter='\t')
        for row in table:
            if row[3] == '0':
                origin_freq_dict[(row[0], row[2])] = float(row[4])
    completed_freq_dict = {}
    for key, value in origin_freq_dict.items():
        completed_freq_dict[reverse_complement(key[0]), reverse_complement(key[1])] = value
        completed_freq_dict[key] = value
    return completed_freq_dict


def reverse_complement(sequence: str):
    """
    :param sequence: (str) DNA string
    :return: (str) reverse-complement sequence DNA string
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'n': 'n'}
    complement_sequence = []
    for nucleotide in sequence:
        complement_sequence.append(complement[nucleotide])
    complement_sequence.reverse()
    reverse_complement_sequence = ''.join(map(str, complement_sequence))
    return reverse_complement_sequence


def sequence_to_subsequences(sequence: list, window_size: int):
    """
    splits the  sequence into subsequences of the specified size
    :param sequence: (list) of values
    :param window_size: (int)
    :return: subsequences (list) of lists with the desired size.
    if the remainder of the division of the number of values in the sequence by
     the window size is not zero, the last subsequence's size increases.
    """
    num_of_windows = len(sequence) // window_size
    length_to_split = [window_size for _ in range(num_of_windows)]
    mod = divmod(len(sequence), window_size)[1]
    if mod != 0:
        length_to_split[-1] += divmod(len(sequence), window_size)[1]
    subsequences = [list(islice(sequence, 0, length_to_split[0]))]
    cum_sum_of_len = length_to_split[0]
    for i in range(1, num_of_windows):
        subsequences.append(list(islice(sequence, cum_sum_of_len, cum_sum_of_len + length_to_split[i])))
        cum_sum_of_len += length_to_split[i]
    return subsequences


def coordinates_per_window(fasta_name: str, window_size: int):
    """
    coordinates of the windows for which the calculation will be carried out
    :param fasta_name:  (str) name of fasta file
    :param window_size: (int) the size of the desired division of the gene into windows
    :return: (list) of (tuples) 1-based coordinates of scanning windows
    """
    gene_name, chrm, strand, list_of_all_things = fasta_parsing(fasta_name)
    coords = [element[1] for element in list_of_all_things]
    list_of_coords = []
    for coord in coords:
        list_of_coords.extend(i for i in range(coord[0] + 1, coord[1]))

    coords = sequence_to_subsequences(list_of_coords, window_size)
    coords_of_windows = [(subsequence[0], subsequence[-1]) for subsequence in coords]

    return coords_of_windows


def frequencies_count(fasta_name: str, freq_table_name: str):
    """
    counting of mutation frequencies per codons
    :param fasta_name: (str) name of FASTA file
    :param freq_table_name: (str) name of mutation rate table (tsv)
    :return: syn_freqs_per_codon, miss_freqs_per_codon, lof_freqs_per_codon (lists)
    of frequencies per codons
    """
    gencode = {
        'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
        'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
        'TTA': 'L', 'TCA': 'S', 'TAA': 'Stop', 'TGA': 'Stop',
        'TTG': 'L', 'TCG': 'S', 'TAG': 'Stop', 'TGG': 'W',

        'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
        'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
        'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
        'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',

        'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
        'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
        'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
        'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',

        'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
        'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
        'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
        'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G',
    }

    other_variants = {
        'A': ['T', 'G', 'C'],
        'T': ['A', 'G', 'C'],
        'G': ['A', 'T', 'C'],
        'C': ['A', 'T', 'G']
    }

    gene_name, _, strand, list_of_all_things = fasta_parsing(fasta_name)
    freq_dict = freq_table_to_dict(freq_table_name)
    sequence = [element[2] for element in list_of_all_things]

    if strand == '-':
        sequence = [reverse_complement(subsequence) for subsequence in sequence]

    syn_freqs_per_codon = []
    miss_freqs_per_codon = []
    lof_freqs_per_codon = []

    """
    Start of counting:
    """
    sequence_index = 1
    for j in range(len(sequence)):  # needed for boundary conditions processing
        subsequence = sequence[j]
        len_of_subs = len(subsequence)
        for index in range(1, len_of_subs - 1):
            syn_freqs_per_codon.append(0)
            miss_freqs_per_codon.append(0)
            lof_freqs_per_codon.append(0)

            cur_nucl = subsequence[index]
            context = f'{subsequence[index - 1]}{cur_nucl}{subsequence[index + 1]}'

            for var_nucl in other_variants[cur_nucl]:
                cur_codon = None
                var_codon = None
                reminder = divmod(sequence_index, 3)[1]

                # some left boundary conditions:
                if index == 1 and reminder == 1:  # (NNN) | (XNN)
                    cur_codon = cur_nucl + subsequence[index + 1] + subsequence[index + 2]
                    var_codon = var_nucl + subsequence[index + 1] + subsequence[index + 2]
                elif index == 1 and reminder == 2:  # (N | XN)(NNN)
                    cur_codon = sequence[j-1][-2] + cur_nucl + subsequence[index + 1]
                    var_codon = sequence[j-1][-2] + var_nucl + subsequence[index + 1]
                elif index == 1 and reminder == 0:  # (NN | X)(NNN)
                    cur_codon = sequence[j-1][-3] + sequence[j-1][-2] + cur_nucl
                    var_codon = sequence[j-1][-3] + sequence[j-1][-2] + var_nucl

                elif index == 2 and reminder == 0:  # (N | NX)(NNN)
                    cur_codon = sequence[j - 1][-2] + subsequence[index - 1] + cur_nucl
                    var_codon = sequence[j - 1][-2] + subsequence[index - 1] + var_nucl

                # some right boundary conditions:
                elif index == (len_of_subs - 3) and reminder == 1:  # (NNN)(XN | N)
                    cur_codon = cur_nucl + subsequence[index + 1] + sequence[j + 1][1]
                    var_codon = var_nucl + subsequence[index + 1] + sequence[j + 1][1]

                elif index == (len_of_subs - 2) and reminder == 2:  # (NNN)(NX | N)
                    cur_codon = subsequence[index - 1] + cur_nucl + sequence[j + 1][1]
                    var_codon = subsequence[index - 1] + var_nucl + sequence[j + 1][1]

                elif index == (len_of_subs - 2) and reminder == 1:  # (NNN)(X | NN)
                    if j != len(sequence) - 1:  # temp
                        cur_codon = cur_nucl + sequence[j + 1][1] + sequence[j + 1][2]
                        var_codon = var_nucl + sequence[j + 1][1] + sequence[j + 1][2]
                    else:
                        raise IndexError('Incomplete coding sequence: '
                                         'the coding sequence does not end with a last codon end')

                # standard cases:
                elif reminder == 1:
                    cur_codon = cur_nucl + subsequence[index + 1] + subsequence[index + 2]
                    var_codon = var_nucl + subsequence[index + 1] + subsequence[index + 2]

                elif reminder == 2:
                    cur_codon = subsequence[index - 1] + cur_nucl + subsequence[index + 1]
                    var_codon = subsequence[index - 1] + var_nucl + subsequence[index + 1]

                elif reminder == 0:
                    cur_codon = subsequence[index - 2] + subsequence[index - 1] + cur_nucl
                    var_codon = subsequence[index - 2] + subsequence[index - 1] + var_nucl

                # comparing and filling in frequency lists
                if gencode[cur_codon] == gencode[var_codon]:
                    syn_freqs_per_codon[-1] += freq_dict[(context, var_nucl)]

                elif gencode[cur_codon] != gencode[var_codon] and gencode[var_codon] == 'Stop':
                    lof_freqs_per_codon[-1] += freq_dict[(context, var_nucl)]

                elif gencode[cur_codon] != gencode[var_codon]:
                    miss_freqs_per_codon[-1] += freq_dict[(context, var_nucl)]
            sequence_index += 1

    return syn_freqs_per_codon, miss_freqs_per_codon, lof_freqs_per_codon


def frequencies_per_window(fasta_name: str, window_size: int):
    """
    returns mutation rates (U) for windows
    :param fasta_name: (str) name of FASTA file containing CDS of gene of interest
    :param window_size: (int) a count of nucleotides in one window
    :return: (tuple) of (lists): mutation rates for synonymous, missense, lost-of-function per window
    """

    dict_name = './data/mutation_rates.tsv'
    syn_freqs_per_codon, miss_freqs_per_codon, lof_freqs_per_codon = frequencies_count(fasta_name, dict_name)

    syn_to_window = sequence_to_subsequences(syn_freqs_per_codon, window_size)
    miss_to_window = sequence_to_subsequences(miss_freqs_per_codon, window_size)
    lof_freqs_per_codon = sequence_to_subsequences(lof_freqs_per_codon, window_size)

    syn_mut_per_window = [sum(subsequence) for subsequence in syn_to_window]
    miss_mut_per_window = [sum(subsequence) for subsequence in miss_to_window]
    lof_mut_per_window = [sum(subsequence) for subsequence in lof_freqs_per_codon]

    return syn_mut_per_window, miss_mut_per_window, lof_mut_per_window


def generate_output(fasta_name: str, window_size: int, out_suffix: str):
    """
    generate .tsv file contains 1- chr, 2- 3- [start, end] coordinates of window (1-based),
    4- strand, 5- Synonymous, 6- Missense, 7- Lost of function mutation rates, 8- Gene name
    :param fasta_name: (str) name of FASTA file containing CDS of gene of interest
    :param window_size: (int) a count of nucleotides in one window
    :param out_suffix: (str) suffix in writing file name
    :return: Message of successfully ended writing to file.
    """
    gene_name, chrm, strand, list_of_all_things = fasta_parsing(fasta_name)
    syn_mut, missense_mut, lof_mut = frequencies_per_window(fasta_name, window_size)
    coords = coordinates_per_window(fasta_name, window_size)
    out_file_name = f'{gene_name}_{out_suffix}.tsv'
    with open(out_file_name, 'w') as csvfile:
        tablewriter = csv.writer(csvfile, delimiter='\t')
        tablewriter.writerow(['Chr', 'Start', 'End', 'Strand', 'FreqSyn', 'FreqMiss', 'FreqLoF', 'Gene'])
        for i in range(len(lof_mut)):
            tablewriter.writerow([chrm, coords[i][0], coords[i][1], strand, syn_mut[i],
                                  missense_mut[i], lof_mut[i], gene_name])
    return 'File is ready!'
