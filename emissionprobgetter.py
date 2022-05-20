from scipy import integrate
from scipy.stats import poisson
from scipy.stats import invgauss
import csv


def pois_function(s, n: int, nu: float):
    """
    P(n | S_het, nu)
    :param s: integrable parameter, selective effect against heterozygous PTV
    :param n: (int) observed genic PTV counts
    :param nu: (float) expected genic PTV counts
    :return: Poisson probability mass function for integration
    """
    return poisson.pmf(n, nu / s)


def ig_function(s, alpha: float, beta: float):
    """
    P(S_het | alpha, beta)
    :param s: integrable parameter, selective effect against heterozygous PTV
    :param alpha: (float) mean-parameter of inverse Gaussian distribution
    :param beta: (float) shape-parameter of inverse Gaussian distribution
    :return: inverse Gaussian probability density function
    """
    mu = alpha / beta
    scale = beta
    return invgauss.pdf(s, mu=mu, scale=scale)


def integrated_function(s, n: int, nu: float, alpha: float, beta: float):
    """
    function that is integrated
    :param s: integrable parameter, selective effect against heterozygous PTV
    :param n: (int) observed genic PTV counts
    :param nu: (float) expected genic PTV counts
    :param alpha: mean-parameter of inverse Gaussian distribution
    :param beta: shape-parameter of inverse Gaussian distribution
    :return: callable function that will be integrated
    """
    return pois_function(s, n, nu) * ig_function(s, alpha, beta)


def get_emission_prob(n: int, nu: float, alpha: float, beta: float, s_het: str):
    """
    numerically counts of integral which is emission probability
    :param n: (int) observed genic PTV counts
    :param nu: (float) expected genic PTV counts
    :param alpha: (float) mean-parameter of inverse Gaussian distribution
    :param beta: (float) shape-parameter of inverse Gaussian distribution
    :param s_het: (str)
    :return: emission_prob: (float) numerically counted integral
    """
    if s_het == 's_not':
        a = 0
        b = 0.01
    elif s_het == 's_cons':
        a = 0.01
        b = 1
    else:
        raise NameError('S_het should be "s_cons" or "s_not"')

    emission_prob = integrate.quad(func=integrated_function, a=a, b=b, args=(n, nu, alpha, beta))[0]
    return emission_prob


def get_alpha_beta(gene_name: str):
    """
    takes alpha and beta parameters from file alpha_beta_for_genes.tsv (part of
    comprehensive_gene_annotation.tsv, provided by Yury Barbitoff)
    :param gene_name: (str) gene short name, like 'A1BG' or 'SHH'
    :return: tuple of alpha (float) and beta (float) parameters of inverse Gaussian distribution
    """
    alpha = None
    beta = None
    with open('data/alpha_beta_for_genes.tsv', 'r') as sourse_file:
        table = csv.reader(sourse_file, delimiter='\t')
        for row in table:
            if row[0] == gene_name:
                alpha = float(row[1])
                beta = float(row[2])
    if alpha is None:
        raise NameError('Can not find alpha, beta for this gene')
    return alpha, beta


def emission_probabilities(u_list: list, allele_number_list: list, n_list: list, gene_name: str):
    """
    calculates the emission probabilities required for the Viterbi algorithm to work.
    :param u_list: (list): mutation rates for lost-of-function mutation rates per window
    :param allele_number_list: (list) of allele counts per window
    :param n_list: (list) of observed genic PTV counts
    :param gene_name:  (str) gene short name, like 'A1BG' or 'SHH'
    :return: (tuple): e_cons_s, e_not_s are (lists) of emission probabilities
    for Conservative and Not Conservative states in HMM model
    """

    e_cons_s = []
    e_not_s = []
    alpha, beta = get_alpha_beta(gene_name)
    for i in range(len(u_list)):
        u_rate = u_list[i]
        allele_number = allele_number_list[i]
        n = n_list[i]
        nu = allele_number * u_rate

        e_cons = get_emission_prob(n=n, nu=nu, alpha=alpha, beta=beta, s_het='s_cons')
        e_not = get_emission_prob(n=n, nu=nu, alpha=alpha, beta=beta, s_het='s_not')

        # normalize probabilities by one, e_cons_i + e_not_i = 1:
        if (e_cons + e_not) == 0:
            e_cons_s.append(0.5)
            e_not_s.append(0.5)
        else:  # CELSR1 gene!
            e_cons_s.append(e_cons / (e_cons + e_not))
            e_not_s.append(e_not / (e_cons + e_not))

    return e_cons_s, e_not_s
