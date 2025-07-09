import numpy as np
import scipy.sparse as sp
import pandas as pd
from scipy.sparse.linalg import eigs
import torch



def select_edge_info(mut_filtered):
    co_mut = np.dot(mut_filtered.transpose(), mut_filtered)
    sum_co_mut = np.sum(co_mut, axis=1)
    mut_diag = np.diagonal(co_mut)
    overlap_time = sum_co_mut / mut_diag
    sample_idx = np.where(overlap_time >= 0.05 * mut_filtered.shape[1])[0]
    return mut_filtered[:, sample_idx]

def filter_data(gene_feature, gene_name, exp_data, met_data, sample_id, ppi, mut_matrix, exp_normal, met_normal, sample_cluster):
    com_gene = np.intersect1d(ppi[:, 0:1].flatten(), gene_name.flatten())
    com_idx = [gene in com_gene for gene in gene_name]
    node_name = gene_name[com_idx, :]
    mut_filtered =  mut_matrix[com_idx, :]
    col_sum = np.sum(mut_filtered, axis=0)
    col_idx = [(10 <= cs <= gene_name.shape[0] * 0.15) for cs in col_sum]
    mut_filtered = mut_filtered[:, col_idx]

    feature_filtered = np.nan_to_num(gene_feature[com_idx, :])
    exp_filtered = exp_data[np.ix_(com_idx, col_idx)]
    met_filtered = met_data[np.ix_(com_idx, col_idx)]
    sample_filtered = sample_id[col_idx, :]
    exp_normal_filtered = exp_normal[com_idx, :]
    met_normal_filtered = met_normal[com_idx, :]
    # co_mut_matrix = mut_filtered
    co_mut_matrix = select_edge_info((mut_filtered > 0) + 0)
    cluster = sample_cluster[col_idx, :]
    return feature_filtered, node_name, exp_filtered, met_filtered, sample_filtered, co_mut_matrix, exp_normal_filtered, met_normal_filtered, cluster



def sigmoid(x):
    s = 1 / (1 + np.exp(-x))
    return s


def ppi_limitation(ppi, node_name, weight_type):
    if weight_type != 0:
        node_name = node_name[:, 0]
        num_node = node_name.shape[0]
        ppi_row = np.isin(ppi[:, 0], node_name) & np.isin(ppi[:, 1], node_name)
        ppi = ppi[ppi_row, :]
        gene_list = node_name.tolist()
        idx_0 = [gene_list.index(gene) for gene in ppi[:, 0]]
        idx_1 = [gene_list.index(gene) for gene in ppi[:, 1]]
        edge = np.asarray(list(map(int, ppi[:, 2])))
        print("\n\tWeighted PPI network...")
        if weight_type == 1:
            edge = (edge > 0) + 0
        elif weight_type == 2:
            edge = sigmoid(edge / np.max(edge))
        else:
            edge = edge / np.max(edge)
        ppi_network = sp.coo_matrix((edge, (idx_0, idx_1)), shape=(num_node, num_node), dtype=np.float32).toarray()
        print("\tnode number: {}\tedge number: {}".format(num_node, len(idx_0)))
        ppi_network = np.maximum(ppi_network, ppi_network.transpose()).astype(np.float32)
        return ppi_network
    return 0



def normalize_adj(adj, sparse=True):
    # adj = sp.coo_matrix(adj)
    row_sum = adj.sum(axis=1)
    idx_0 = np.where(row_sum == 0.)[0]
    row_sum[idx_0] = 1.
    d_inv_sqrt = np.power(row_sum, -0.5).flatten()
    # d_inv_sqrt[np.isinf(d_inv_sqrt)] = 0
    d_inv_sqrt[idx_0] = 0.
    d_mat_inv_sqrt = np.diagflat(d_inv_sqrt)
    res = adj.dot(d_mat_inv_sqrt).transpose().dot(d_mat_inv_sqrt)
    if sparse:
        return sp.coo_matrix(res)
    else:
        return res



def subtract_lower_support(polys):
    for i in range(1, len(polys)):
        for j in range(0, i):
            if j == 0:
                polys[i][(np.abs(polys[j].todense()) > 0.0001)] = 0
            else:
                polys[i][(np.abs(polys[j]) > 0.0001)] = 0
    return polys



def chebyshev_polynomials(adj, k, sparse=True, subtract_support=True):
    """Calculate Chebyshev polynomials up to order k. Return a list of sparse matrices (tuple representation)."""
    adj_normalized = normalize_adj(adj, sparse=sparse)
    laplacian = sp.eye(adj.shape[0]) - adj_normalized
    largest_eigval, _ = eigs(laplacian)
    scaled_laplacian = (2. / largest_eigval[0].real) * laplacian - sp.eye(adj.shape[0])

    t_k = list()
    if sparse:
        t_k.append(sp.eye(adj.shape[0]))
    else:
        t_k.append(np.eye(adj.shape[0]))
    t_k.append(scaled_laplacian)

    def chebyshev_recurrence(t_k_minus_one, t_k_minus_two, scaled_lap):
        s_lap = sp.csr_matrix(scaled_lap, copy=True)
        return 2 * s_lap.dot(t_k_minus_one) - t_k_minus_two

    for _ in range(2, k + 1):
        t_k.append(chebyshev_recurrence(t_k[-1], t_k[-2], scaled_laplacian))

    if subtract_support:
        t_k = subtract_lower_support(t_k)
    return t_k



def get_support_matrices(adj, poly_support):
    if poly_support > 0:
        support = chebyshev_polynomials(adj, poly_support)
        num_supports = 1 + poly_support
    else:
        support = [np.eye(adj.shape[0])]
        num_supports = 1
    ppi_graph = support[0]
    for i in range(1, num_supports):
        ppi_graph += support[i]
    return ppi_graph



def get_hypergraph_matrix(gene_name, co_data_net):
    slice_num = len(co_data_net)
    hp_graph = list()
    for s in range(slice_num):
        slice_s = co_data_net[s]
        hp_graph.append(get_support_matrices(slice_s, poly_support=1))

    slice_num = len(hp_graph)
    hp_graph_arr = hp_graph[0].toarray()[:, :, np.newaxis]

    for s in range(1, slice_num):
        hp_graph_arr = np.concatenate((hp_graph_arr, hp_graph[s].toarray()[:, :, np.newaxis]), axis=2)

    print("\t\tThe number of genes in Hypergraph: {}.".format(gene_name.shape[0]))
    print("\t\tThe number of slices in Hypergraph: {}.".format(slice_num))
    return hp_graph_arr


def calculate_co_association(tumor_data, normal_data, gene_name, sample_cluster, net_split):
    gene_num = gene_name.shape[0]
    cluster_num = sample_cluster[:, 0].astype(np.int32)
    cor_nor = np.corrcoef(normal_data + 1e-4 * np.random.normal(0, 1, [gene_num, normal_data.shape[1]]))
    diff_co_data = list()
    for cluster in range(1, np.max(cluster_num) + 1):
        sample_idx = np.where(cluster_num == cluster)[0]
        tumor_cluster = tumor_data[:, sample_idx]
        cor_cluster = np.corrcoef(tumor_cluster + 1e-4 * np.random.normal(0, 1, [gene_num, sample_idx.shape[0]]))
        if net_split <= 0:
            diff_cor = ((np.abs(cor_cluster) > 0.6) ^ (np.abs(cor_nor) > 0.6)) + 0
        else:
            diff_idx = (np.abs(cor_cluster-cor_nor) > net_split) + 0
            diff_cor = np.abs(cor_cluster-cor_nor) * diff_idx
        diff_co_data.append(diff_cor.astype(np.float32))
    return diff_co_data


def construct_association_network(diff_co_data, ppi_network):
    cluster_num = len(diff_co_data)
    co_data_net = list()
    print("\n\tEdge number in differential co-association network...")
    for cluster in range(cluster_num):
        diff_cor = diff_co_data[cluster]
        co_data_net.append(diff_cor * ppi_network)
        edge_num = np.sum(co_data_net[cluster])
        print("\t\tFor sample cluster {}: {}".format(cluster, edge_num))
    return co_data_net
