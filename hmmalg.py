import numpy as np


def forward(v, a, b, initial_distribution):
    """

    :param v:  (np.array) of observations (allele counts per window)
    :param a: (np.array) of transition probabilities at the i-th step
    :param b: (np.array) of emission probabilities
    :param initial_distribution: (np.array) of transition probabilities of first state
    :return: alpha (np.array) of auxiliary quantities for the Baum-Welch algorithm
    """
    alpha = np.zeros((v.shape[0], a.shape[0]))
    alpha[0, :] = initial_distribution * b[:, 0]

    for t in range(1, v.shape[0]):
        for j in range(a.shape[0]):
            alpha[t, j] = alpha[t - 1] @ a[:, j] * b[j, t]

    return alpha


def backward(v, a, b):
    """

    :param v:  (np.array) of observations (allele counts per window)
    :param a: (np.array) of transition probabilities at the i-th step
    :param b: (np.array) of emission probabilities
    :return: beta (np.array) of auxiliary quantities for the Baum-Welch algorithm
    """
    beta = np.zeros((v.shape[0], a.shape[0]))

    beta[v.shape[0] - 1] = np.ones((a.shape[0]))
    for t in range(v.shape[0] - 2, -1, -1):
        for j in range(a.shape[0]):
            beta[t, j] = (beta[t + 1] * b[:, t + 1]) @ a[j, :]

    return beta


def baum_welch_algo(observations, tran_vec, em_vec, start_transitions, n_iter=100):
    """
    returns transition probabilities corresponding to the maximum likelihood (searches for a local maximum)
    :param observations: (list) of observations (allele counts per window)
    :param tran_vec: (list) of zero approximation for transition probabilities
    :param em_vec: (list) of emission probabilities
    :param start_transitions: (list) of transition probabilities of first state
    :param n_iter: (int) number of iterations of the procedure
    :return: (list) of maximum-likelihood transition probabilities
    """
    observations = np.array(observations)
    tran_vec = np.array(tran_vec)
    em_vec = np.array(em_vec)
    start_transitions = np.array(start_transitions)
    m = tran_vec.shape[0]
    t_max = len(observations)
    for n in range(n_iter):
        # estimation step
        alpha = forward(observations, tran_vec, em_vec, start_transitions)
        beta = backward(observations, tran_vec, em_vec)
        xi = np.zeros((m, m, t_max - 1))
        for t in range(t_max - 1):
            denominator = (alpha[t, :].T @ tran_vec * em_vec[:, t + 1].T) @ beta[t + 1, :]
            for i in range(m):
                numerator = alpha[t, i] * tran_vec[i, :] * em_vec[:, t + 1].T * beta[t + 1, :].T
                xi[i, :, t] = numerator / denominator
        gamma = np.sum(xi, axis=1)
        # maximization step
        tran_vec = np.sum(xi, 2) / np.sum(gamma, axis=1).reshape((-1, 1))
    return tran_vec


def viterbi_algo(observations: list, states: tuple, start_transitions: dict,
                 emissions: list, transitions: dict):
    """

    :param observations: (list) of observations (allele counts per window)
    :param states: (tuple) states of hidden Markov Model
    :param start_transitions: (dict) of transition probabilities of first state
    :param emissions: (list) of emission probabilities
    :param transitions: (dict) of transition probabilities
    :return: (list) of maximum-likelihood optimal states' pathway and (float) maximum-likelihood probability
    """
    # all existing paths:
    v = [{}]
    # boundary condition (start):
    for i in range(len(states)):
        st = states[i]
        v[0][st] = {'prob': start_transitions[st] * emissions[i][0], 'prev': None}

    for k in range(1, len(observations)):
        v.append({})
        for i in range(len(states)):
            st = states[i]
            max_prob = 0
            prev_st_selected = None
            for j in range(len(states)):
                prev_st = states[j]
                tr_prob = v[k - 1][prev_st]['prob'] * transitions[prev_st][st] * emissions[i][k]
                if tr_prob > max_prob:
                    max_prob = tr_prob
                    prev_st_selected = prev_st

            v[k][st] = {'prob': max_prob, 'prev': prev_st_selected}

    opt = []
    max_prob = 0.0
    best_st = None
    # Get most probable state and its backtrack
    for st, data in v[-1].items():
        if data['prob'] > max_prob:
            max_prob = data['prob']
            best_st = st
    opt.append(best_st)
    previous = best_st

    # Follow the backtrack till the first observation
    for k in range(len(v) - 2, -1, -1):
        opt.insert(0, v[k + 1][previous]['prev'])
        previous = v[k + 1][previous]['prev']

    return opt, max_prob


def decoder(lof_first_window: float, trans_p: tuple,
            observations: list, emit_p_cons: list, emit_p_not: list):
    """

    :param lof_first_window: allele count per window size if first window
    :param trans_p: (tuple) of zero approximation for transition probabilities
    :param observations: (list) of observations (allele counts per window)
    :param emit_p_cons: (list) of emission probabilities for Conservative state in HMM model
    :param emit_p_not: (list) of emission probabilities for Not Conservative states in HMM model
    :return: (list) of maximum-likelihood optimal states' pathway  and (float) maximum-likelihood probability
    """
    states = ('Cons', 'Not')

    # the formation of this lists depends on the order of states:
    transitions = [[trans_p[0], 1 - trans_p[0]],
                   [1 - trans_p[1], trans_p[1]]]
    emissions = [emit_p_cons, emit_p_not]
    start_p = [1 - lof_first_window, lof_first_window]

    trans_p = baum_welch_algo(observations, transitions, emissions, start_p)
    start_p = {'Cons': 1 - lof_first_window, 'Not': lof_first_window}
    transitions = {'Cons': {'Cons': trans_p[0][0], 'Not': trans_p[0][1]},
                   'Not': {'Cons': trans_p[1][0], 'Not': trans_p[1][1]}}
    path, prob = viterbi_algo(observations=observations, states=states, start_transitions=start_p,
                              emissions=emissions, transitions=transitions)
    return path, prob
