# single network fwd simulation

import sys
import numpy as np
import casadi as cas
from .network import KineticModel
import time


def _sample(dE, transi_matrix, sample_method):
    if sample_method == 'elementwise':
        n = len(dE)
        i = int(n * np.random.uniform())
        scale = np.random.normal(0,1)
        deltaE = scale * transi_matrix[:, i]
    if sample_method == 'augment':
        n = len(dE)
        scale = np.random.normal(0, 1, n)
        deltaE = transi_matrix.dot(scale)
    newE = np.array(dE) + np.array(deltaE)
    return list(newE)


def evalZ(dE, reactorlist, condi_reactor_list, prior_list, evidence_list):
    log_prior_list, log_likeli_list = [], []
    for reactor, conditionlist, prior_info, evidence_info in zip(reactorlist,
                                                                 condi_reactor_list,
                                                                 prior_list,
                                                                 evidence_list):
        log_likeli, log_prior = reactor.eval_prob(dE, conditionlist, evidence_info, prior_info)
        log_prior_list.append(log_prior)
        log_likeli_list.append(log_likeli)
    assert all(prior == log_prior_list[0] for prior in log_prior_list)
    return sum(log_likeli_list) + log_prior_list[0]

def sn_hessian_cal(dE, reactorlist, condi_reactor_list, prior_list, evidence_list, est=False, delta=0.1):

    E0 = np.array(dE)
    Z0 = evalZ(dE, reactorlist, condi_reactor_list, prior_list, evidence_list)

    hess = np.zeros((len(dE), len(dE)))

    Zlist = np.zeros_like(dE)
    for i in range(len(E0)):
        Ei = np.copy(E0)
        Ei[i] += delta
        Zlist[i] = evalZ(Ei, reactorlist, condi_reactor_list, prior_list, evidence_list)

    if est:
        grad = (Zlist - Z0) / delta
    else:
        for i in range(len(E0)):
            print('------ %2dth degree of freedom ------' %i)
            for j in range(len(E0)):
                # construct delta energy
                Eij = np.copy(E0)
                Eij[i] += delta
                Eij[j] += delta

                Zij = evalZ(Eij, reactorlist, condi_reactor_list, prior_list, evidence_list)

                hess[i, j] = (Zij - Zlist[i] - Zlist[j] + Z0) / (delta**2)
        hess = (hess.T + hess) / 2.
    return hess

def evaluate_prob(dE, reactorlist, condi_reactor_list, prior_list, evidence_list):
    log_prior_list, log_likeli_list = [], []
    for reactor, conditionlist, prior_info, evidence_info in zip(reactorlist,
                                                                 condi_reactor_list,
                                                                 prior_list,
                                                                 evidence_list):
        log_likeli, log_prior = reactor.eval_prob(dE, conditionlist, evidence_info, prior_info)
        log_prior_list.append(log_prior)
        log_likeli_list.append(log_likeli)
    assert all(prior == log_prior_list[0] for prior in log_prior_list)
    return log_prior_list, log_likeli_list


def sn_bayesian_infer(reactorlist, condi_reactor_list, prior_list, evidence_list,
                      ntot, nbuf, dE_start, transi_matrix,
                      sample_method='elementwise', constraint=True,
                      step_write=100, save_result=True, save_step=10000, report=None):
    '''
    Bayesian Analysis on Sinlge reaction network chemistry with different type of reactors
    '''
    reactor = reactorlist[0]
    prior_info = prior_list[0]

    out = ''
    # Summary of Evidence and Prior
    out += '\n\n'
    out += 'Transition Matrix:  \n'
    out += str(transi_matrix) + '\n'
    out += '==' * 20 + '\n'
    out += 'Total Step = ' + str(ntot) + '\n'
    out += 'Buffer Step = ' + str(nbuf) + '\n'
    out += 'Save Result: ' + str(save_result) + '\n'
    out += 'Starting Point = ' + str(dE_start) + '\n'
    out += 'Thermo Constraint (True: check, False: not check) = ' + str(constraint) + '\n'
    out += '==' * 20 + '\n'

    tic = time.clock()

    # prior and posterior calculation
    log_prior_list, log_likeli_list = evaluate_prob(dE_start, reactorlist, condi_reactor_list, prior_list, evidence_list)
    log_prior_prev = log_prior_list[0]
    log_likeli_prev = sum(log_likeli_list)
    log_posterior_prev = log_prior_prev + log_likeli_prev

    header = '{0:^10s}  {1:^15s}  {2:^15s}  {3:^15s}  {4:^15s}  {5:^15s}'. \
               format('step', 'log(prior)', 'log(likelihood)', 'log(posterior)', 'accept%', 'infeasi%') + '\n'
    header += '==' * 30
    out += header + '\n'
    print(header)


    Edis = []
    Eprev = np.copy(dE_start)

    jump, infeasi = 0, 0
    for i in range(ntot):
        if i % step_write == 0:
            output = '{0:^10d}  {1:^15.2f}  {2:^15.2f}  {3:^15.2f}  {4:^15.2f}  {5:^15.2f}'. \
                format(i, log_prior_prev, log_likeli_prev, log_posterior_prev,
                       jump / float(i + 1) * 100, infeasi / float(i + 1) * 100)
            out += output + '\n'
            print(output)

        Estar = _sample(Eprev, transi_matrix, sample_method)

        # Define Constarint: True: Satisfy // False: not Satisfy
        thermoConstraint = reactor.CheckThermoConsis(Estar) or (not constraint)
        boundConstraint = len(np.argwhere(np.array(Estar) - np.array(prior_info['lbound']) < 0)) == 0 and \
                          len(np.argwhere(np.array(Estar) - np.array(prior_info['ubound']) > 0)) == 0

        if (not thermoConstraint) or (not boundConstraint):
            # REJECT
            pass
        else:
            try:
                log_prior_list, log_likeli_list = evaluate_prob(Estar, reactorlist, condi_reactor_list, prior_list, evidence_list)
                log_prior_star = log_prior_list[0]
                log_likeli_star = sum(log_likeli_list)
                log_posterior_star = log_prior_star + log_likeli_star
                UU = np.random.uniform()
                AA = min(1, np.exp(log_posterior_star - log_posterior_prev))
                # AA = min(1, posterior_star / posterior_prev)
                if UU < AA:
                    Eprev = np.copy(Estar)
                    log_prior_prev = log_prior_star
                    log_likeli_prev = log_likeli_star
                    log_posterior_prev = log_posterior_star

                    jump += 1
                else:
                    # REJECT
                    pass
            except:
                infeasi += 1
                pass
        if i >= nbuf:
            Edis.append(Eprev)
            if i % save_step == 0:
                np.save('edis_breakpoint_%d.npy' % i, Edis)

    out += '==' * 30 + '\n'
    toc = time.clock()
    out += 'Total time = %.2f hr\n' % ((toc - tic) / 3600.)

    if report is not None:
        with open(report, 'a') as fp:
            fp.write(out)
            fp.close()
    return Edis


