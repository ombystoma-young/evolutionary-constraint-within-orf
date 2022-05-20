from emissionprobgetter import emission_probabilities
from freqcounter import frequencies_per_window, coordinates_per_window, get_gene_name
from ploffinder import observed_ptv
from hmmalg import decoder


def conservativeness_estimator(genefilename: str, variantsfilename: str, window_size: int, transitions: tuple):
    """
    takes .fasta file, .tsv file with variants for one gene, parameters of hmm model
    and returns most likelihood states for this gene
    :param genefilename: (str) name of FASTA file containing CDS of gene of interest
    :param variantsfilename: (str) name of TSV file containing protein-truncating variants of gene of interest
    :param window_size: (int) desired window size
    :param transitions: (tuple): desired transitions probabilities (Cons->Cons, Not->Not)
    :return: (tuple) result: (list) path of states, prob: (float) likelihood of this states path.
    """
    u_list = frequencies_per_window(genefilename, window_size)[2]
    coords_list = coordinates_per_window(genefilename, window_size)
    n_list, an_list = observed_ptv(variantsfilename, coords_list)
    gene_name = get_gene_name(genefilename)

    emissions_cons, emissions_not = emission_probabilities(u_list=u_list,
                                                           allele_number_list=an_list,
                                                           n_list=n_list,
                                                           gene_name=gene_name)
    plof_per_window_len = n_list[0] / window_size
    path, _ = decoder(lof_first_window=plof_per_window_len, trans_p=transitions,
                      observations=n_list, emit_p_cons=emissions_cons,
                      emit_p_not=emissions_not)

    return gene_name, coords_list, u_list, n_list, an_list, path
