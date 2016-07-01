#!/usr/bin/env python3
import copy
import sys


import numpy as np

###
# Load data
###
with open('fake_data.txt') as f:
    header = f.readline()
    cols = header.rstrip().split('\t')
    samples = cols[1:]
    rows = []
    counts = []
    field_gen = (row.rstrip().split('\t') for row in f)
    rows =[(r[0], [float(e) for e in r[1:]]) for r in field_gen]
    otus, counts = zip(*rows)
    counts_array = np.array(counts).transpose()

###
# Perform SparCC correlation estimation
###
cor_list = []
cov_list = []
for i in range(20):
    print('Running iteration', i)
    #Step 1 (estimate component fractions):
    np.random.seed(18724)
    #counts = np.load('input_table.npy')
    n, m = np.shape(counts_array)
    fractions = np.matrix(np.ones((n, m)))
    fractions /= fractions.sum(axis = 1)
    for i in range(n):
            count_row = counts_array[i, ]
            dirichlet_param = count_row + 1
            fractions[i, ] = np.random.mtrand.dirichlet(dirichlet_param)

    #Step 2 (calculate log-ration variances):
    variance_array = np.zeros((m ,m))
    for i in range(m-1):
        for j in range(i+1, m):
            log_ratio = np.array(np.log(np.array(fractions)[: ,i]/np.array(fractions)[: ,j]))
            variance_array[i,j] = variance_array[j,i] = np.var(log_ratio, ddof=1)

    variance_matrix = np.matrix(variance_array)


    #Step 3 (calculate component variations):
    k = len(variance_matrix)
    variance_vector = variance_matrix.sum(axis = 1)
    cov_array = np.matrix(np.zeros((k, k)))
    cov_vector = (cov_array - np.diag(np.diag(cov_array))).sum(axis=1)
    mod_array = np.matrix(np.ones((k, k)) + np.diag([m-2]*m))
    basis_variance = mod_array.I * (variance_vector + 2 * cov_vector)
    # Minimum variance as argument?
    basis_variance[basis_variance < 0] = 1e-4

    #Step 4 (calculate cor and cov from component and fraction variation):
    k = len(basis_variance)
    basis_covariance = np.matrix(np.zeros((k, k)))
    basis_correlation = np.matrix(np.zeros((k, k)))

    # Combinations of all basis_variance rows plus variance_matrix[i,j]
    for i in range(k-1):
        basis_covariance[i, i] = basis_variance[i]
        basis_correlation[i, i] = 1.0
        for j in range(i+1, k):
            basis_cov_elem = 0.5*(basis_variance[i] + basis_variance[j] - variance_matrix[i,j])
            basis_cor_elem = basis_cov_elem / np.sqrt(basis_variance[i]) / np.sqrt(basis_variance[j])
            if np.abs(basis_cor_elem) > 1: #check if got valid correlation value (-1 < cor < 1)
                basis_cor_elem = np.sign(basis_cor_elem)
                basis_cov_elem = basis_cov_elem * np.sqrt(basis_variance[i]) * np.sqrt(basis_variance[j])
            basis_covariance[i,j] = basis_cov_elem
            basis_covariance[j,i] = basis_covariance[i,j]
            basis_correlation[i,j] = basis_cor_elem
            basis_correlation[j,i] = basis_correlation[i,j]
            if np.isnan(basis_correlation[i,j]):
                print(basis_variance[i], basis_variance[j], variance_matrix[i,j])

    basis_covariance[k-1,k-1] = basis_variance[k-1]
    basis_correlation[k-1,k-1] = 1.0

    #Step 5 (exclude highly correlated pairs):
    # Argument for xiter and theshold (10, 0.1 respectively)
    excluded = None
    variance_matrix_mun = copy.deepcopy(variance_matrix)
    for i in range(10):
        # Exclude pairs function
        basis_cor_temp = abs(basis_correlation - np.diag(np.diag(basis_correlation)))
        if excluded:
            basis_cor_temp[excluded] = 0
        exclude_idx_arrays = [np.atleast_2d(r) for r in np.where(basis_cor_temp==basis_cor_temp.max())]
        i = exclude_idx_arrays[0][0,1]
        j = exclude_idx_arrays[1][0,1]
        exclude = (np.matrix([i,j]), np.matrix([j,i]))
        if basis_cor_temp.max() > 0.1:
            mod_array[exclude] -= 1
            for i in exclude[0]:
                mod_array[i,i] -= 1
            if excluded:
                excluded = (np.r_[excluded[0],exclude[0]], np.r_[excluded[1],exclude[1]])
            else:
                excluded = exclude
        else:
            # If we're out pairs to exclude, break (?)
            break
        ###
        # Recalculate component variations after excluding ----
        k = len(variance_matrix_mun)
        if excluded:
            variance_matrix_mun[excluded] = 0
        variance_vector = variance_matrix_mun.sum(axis = 1)
        cov_array = np.matrix(np.zeros((k, k)))
        cov_vector = (cov_array - np.diag(np.diag(cov_array))).sum(axis=1)
        basis_variance = mod_array.I * (variance_vector + 2 * cov_vector)
        # Minimum variance as argument?
        basis_variance[basis_variance < 0] = 1e-4
        ###
        # Recalculate correlation and covariance ----
        k = len(basis_variance)
        basis_covariance = np.matrix(np.zeros((k, k)))
        basis_correlation = np.matrix(np.zeros((k, k)))
        # Combinations of all basis_variance rows plus variance_matrix[i,j]
        for i in range(k-1):
            basis_covariance[i, i] = basis_variance[i]
            basis_correlation[i, i] = 1.0
            for j in range(i+1, k):
                basis_cov_elem = 0.5*(basis_variance[i] + basis_variance[j] - variance_matrix[i,j])
                basis_cor_elem = basis_cov_elem / np.sqrt(basis_variance[i]) / np.sqrt(basis_variance[j])
                if np.abs(basis_cor_elem) > 1: #check if got valid correlation value (-1 < cor < 1)
                    basis_cor_elem = np.sign(basis_cor_elem)
                    basis_cov_elem = basis_cov_elem * np.sqrt(basis_variance[i]) * np.sqrt(basis_variance[j])
                basis_covariance[i,j] = basis_cov_elem
                basis_covariance[j,i] = basis_covariance[i,j]
                basis_correlation[i,j] = basis_cor_elem
                basis_correlation[j,i] = basis_correlation[i,j]
                if np.isnan(basis_correlation[i,j]):
                    print(basis_variance[i], basis_variance[j], variance_matrix[i,j])
        basis_covariance[k-1,k-1] = basis_variance[k-1]
        basis_correlation[k-1,k-1] = 1.0
        # If conditional is meet, calculate as center-log ratio (this doesn't seem
        # to happen often)
        if np.max(np.abs(basis_correlation)) > 1.0:
            n,k  = fractions.shape
            temp = np.log(fractions)
            m = np.mean(temp, axis = 1)
            z = temp - np.tile(m,(1,k))
            basis_covariance = np.cov(z, rowvar = 0)
            basis_correlation = np.corrcoef(z, rowvar = 0)
            basis_variance = np.diag(basis_covariance)
    # Finally add basis_variance, basis_correlation, and basis_covariance to
    # list
    cor_list.append(basis_correlation)
    cov_list.append(np.diag(basis_covariance))

###
# Calculate the final tables
###
print('Finishing up')
cor_array = np.array(cor_list)
var_med = np.median(cov_list,axis = 0) #median covariance
cor_med = np.median(cor_array,axis = 0) #median covariance
x,y     = np.meshgrid(var_med,var_med)
cov_med = cor_med * x**0.5 * y**0.5
