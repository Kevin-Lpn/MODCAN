import os
import pickle
import math
from collections import Counter
from multiprocessing import Process

import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import linregress, median_abs_deviation
import networkx as nx
from networkx.algorithms.efficiency_measures import global_efficiency
from sklearn.cluster import AgglomerativeClustering

def check_input_style(input_list):

    comm_file = input_list[0]
    weight_file = input_list[1]
    control_file = input_list[2]
    case_file = input_list[3]

    return comm_file, weight_file, control_file, case_file


def make_community_dict(comm_file, cancer):

    df_comm = pd.read_csv(comm_file, sep='\t')

    # convert gene-community table to dictionary
    comm_dict = {}
    for idx, row in df_comm.iterrows():
        if row["community"] in comm_dict.keys():
            comm_dict[row["community"]].append(row["gene"])
        else:
            comm_dict[row["community"]] = [row["gene"]]

    # save module: gene-list dictionary
    file_path = "./TCGA-{}/comm_dict.pickle".format(cancer)
    with open(file_path, "wb") as f:
        pickle.dump(comm_dict, f)

    print("The number of community: ", len(comm_dict))


def make_edge_list_per_community(weight_file, cancer):

    # load module: gene-list dictionary
    file_path = "./TCGA-{}/comm_dict.pickle".format(cancer)
    with open(file_path, "rb") as f:
        comm_dict = pickle.load(f)

    # load weight matrix
    df_weight = pd.read_csv(weight_file, header=None, skiprows=1)
    # df_weight.columns = df_weight.columns.astype(str)

    with open("./TCGA-{}/gene.txt".format(cancer), "r") as f:
        labels = [line.strip() for line in f.readlines()]
    df_weight.index = labels
    df_weight.columns = labels

    os.mkdir("./TCGA-{}/community_edges".format(cancer))
    for comm in comm_dict:
        path = "./TCGA-{}/community_edges/".format(cancer) + comm
        os.mkdir(path)

    for community in sorted(list(comm_dict.keys())):
        community_genes = sorted(comm_dict[community])
        df_weight_module = df_weight.loc[community_genes, community_genes]

        # make edge tables (gene1, gene2, weight)
        gene1s = []
        gene2s = []
        weights = []
        for i in range(0, len(community_genes)):
            for j in range(i+1, len(community_genes)):
                gene1s.append(community_genes[i])
                gene2s.append(community_genes[j])
                weights.append(df_weight_module.loc[community_genes[i], community_genes[j]])


        df_community_edges = pd.DataFrame(columns=["gene1", "gene2", "weight"])
        df_community_edges["gene1"] = gene1s
        df_community_edges["gene2"] = gene2s
        df_community_edges["weight"] = weights

        top10 = np.percentile(weights, 90)
        # include only top10% (weight) egdes in each module
        df_community_edges_filtered = df_community_edges[df_community_edges["weight"] >= top10]
        # save edge list in each module
        df_community_edges_filtered.to_csv("./TCGA-{}/community_edges/".format(cancer) + community + "/df_edges_filtered.csv", index=False)


def convert_to_median_rank_scores(control_file, case_file, cancer):

    # load gene expression profiles (normalized) of control samples
    exp_control = pd.read_csv(control_file, header=None,  sep='\t')
    with open("./TCGA-{}/gene.txt".format(cancer), "r") as f:
        columns_control = [line.strip() for line in f.readlines()]
    with open("./TCGA-{}/sample_normal.txt".format(cancer), "r") as f:
        index_control = [line.strip() for line in f.readlines()]
    exp_control.index = index_control
    exp_control.columns = columns_control
    exp_control = exp_control.transpose()
    samples_control = list(exp_control.columns)
    genes_control= list(exp_control.index)

    # convert expression values to rank
    rank_control = exp_control.rank(method="max", ascending=True)
    rank_values = np.array(rank_control.values)
    rank_median = np.median(rank_values, axis=1)

    # make gene: median rank dictionary
    gene_rank = {}
    for i in range(0, len(genes_control)):
        gene_rank[genes_control[i]] = rank_median[i]

    sorted_values = []
    for sample in samples_control:
        values = sorted(list(exp_control[sample]))
        sorted_values.append(values)

    # meausure median value per rank
    sorted_values = np.array(sorted_values)
    sorted_values = sorted_values.astype("float32")
    medians = np.median(sorted_values, axis=0)

    # normalize gene expression profiles of control samples using median rank scores
    mrs_reference = pd.DataFrame(index=list(exp_control.index), columns=list(exp_control.columns))
    for i in tqdm(range(0, len(samples_control))):
        ranks = list(rank_control[samples_control[i]])
        counters = Counter(ranks)
        rank_dict = {}
        for rank in counters.keys():
            if int(rank) == 1 and counters[rank] != 1:
                rank_dict[rank] = medians[0]
            elif counters[rank] == 1:
                rank_dict[rank] = medians[int(rank - 1)]
            else:
                medians_filtered = medians[int(rank - counters[rank]):int(rank)]
                medians_mean = np.mean(medians_filtered)
                rank_dict[rank] = medians_mean

        qns = [rank_dict[rank] for rank in ranks]
        mrs_reference[samples_control[i]] = qns

    mrs_reference = mrs_reference.reset_index()
    mrs_reference = mrs_reference.rename(columns={"index": "Symbol"})

    # save normalized gene expression profiles of control samples
    mrs_reference.to_csv("./TCGA-{}/control_exp_normalized.csv".format(cancer), index=False)
    print(mrs_reference.head(5))

    exp_case = pd.read_csv(case_file, header=None, sep='\t')
    with open("./TCGA-{}/gene.txt".format(cancer), "r") as f:
        columns_case = [line.strip() for line in f.readlines()]
    with open("./TCGA-{}/sample_data.txt".format(cancer), "r") as f:
        index_case = [line.strip() for line in f.readlines()]
    exp_case.index = index_case
    exp_case.columns = columns_case
    exp_case = exp_case.transpose()
    samples_case = list(exp_case.columns)
    genes_case= list(exp_case.index)
    # genes_intersections = set(genes_control).intersection(set(genes_case))
    # exp_case = exp_case[sorted(list(genes_intersections))]

    rank_case = exp_case.rank(method="max", ascending=True)

    gene_missed = set(gene_rank.keys()) - set(genes_case)

    # infer ranks of missing genes from control samples
    for gene in gene_missed:
        rank_case.loc[gene, :] = [gene_rank[gene]] * len(samples_case)

    rank_case = rank_case.sort_index(ascending=True)
    rank_case = rank_case.rank(method="max", ascending=True)

    # normalize gene expression profiles of case samples using median rank scores
    mrs_case = pd.DataFrame(index=list(rank_case.index), columns=list(rank_case.columns))
    for i in tqdm(range(0, len(samples_case))):
        ranks = list(rank_case[samples_case[i]])
        counters = Counter(ranks)
        rank_dict = {}
        for rank in counters.keys():
            if counters[rank] == 1:
                rank_dict[rank] = medians[int(rank - 1)]
            else:
                # medians_filtered = medians[int(rank - 1):int(rank - 1 + counters[rank])]
                medians_filtered = medians[int(rank - counters[rank]):int(rank)]
                medians_mean = np.mean(medians_filtered)
                rank_dict[rank] = medians_mean

        qns = [rank_dict[rank] for rank in ranks]
        mrs_case[samples_case[i]] = qns

    mrs_case = mrs_case.reset_index()
    mrs_case = mrs_case.rename(columns={"index": "Symbol"})

    # save normalized gene expression profiles of case samples
    mrs_case.to_csv("./TCGA-{}/case_exp_normalized.csv".format(cancer), index=False)
    print(mrs_case.head(5))


def make_linear_equation(mrs_control, genes, idx_list, cancer):

    for i in tqdm(range(0, len(idx_list))):

        gene1s = []
        gene2s = []
        slopes = []
        intercepts = []
        rs = []
        meds = []
        mads = []

        x = np.array(list(mrs_control[genes[idx_list[i]]]))
        x = x.astype("float")

        for j in range(idx_list[i]+1, len(genes)):

            y = np.array(list(mrs_control[genes[j]]))
            y = y.astype("float")

            slope, intercept, r, p, se = linregress(x, y)

            distances = (slope * x) - y + intercept
            distances = np.abs(distances) / math.sqrt(slope * slope + 1)

            mad = median_abs_deviation(distances)
            med = np.median(distances)

            gene1s.append(genes[idx_list[i]])
            gene2s.append(genes[j])

            meds.append(med)
            mads.append(mad)
            slopes.append(slope)
            intercepts.append(intercept)
            rs.append(abs(r))

        df_linregress = pd.DataFrame(columns=["gene1", "gene2", "slope", "intercept", "r", "med", "mad"])
        df_linregress["gene1"] = gene1s
        df_linregress["gene2"] = gene2s

        df_linregress["slope"] = slopes
        df_linregress["intercept"] = intercepts
        df_linregress["r"] = rs
        df_linregress["med"] = meds
        df_linregress["mad"] = mads

        df_linregress.to_csv("./TCGA-{}/linear_regression_modeling/".format(cancer) + genes[idx_list[i]] + "_linear.csv", index=False)


def multiprocess_measure_linear_regression(cancer):

    os.mkdir("./TCGA-{}/linear_regression_modeling".format(cancer))

    # load normalized gene expression profiles of the case samples
    mrs_control = pd.read_csv("./TCGA-{}/control_exp_normalized.csv".format(cancer))
    mrs_control = mrs_control.set_index("Symbol")
    mrs_control = mrs_control.transpose()

    # load community: gene-list dictionary
    with open("./TCGA-{}/comm_dict.pickle".format(cancer), "rb") as f:
        comm_dict = pickle.load(f)

    assigned_genes = []
    for module in comm_dict.keys():
        assigned_genes += comm_dict[module]

    mrs_control = mrs_control[assigned_genes]

    idxs = [i for i in range(0, len(assigned_genes))]

    split_num = len(assigned_genes) // 1000
    idx_splited = []
    for num in range(0, split_num):
        idx_splited.append(idxs[1000 * (num) : 1000 * (num + 1)])

    idx_splited.append(idxs[1000 * split_num :])

    procs = []
    for index, idx_list in enumerate(idx_splited):
        proc = Process(target=make_linear_equation,
                       args=(mrs_control, assigned_genes, idx_list, cancer))
        procs.append(proc)
        proc.start()

    for proc in procs:
        proc.join()


def estimate_individualized_perturbed_edges(mrs_case, case_samples, genes, cancer):

    for i in tqdm(range(0, len(genes))):
        df_linregress_gene = pd.read_csv("./TCGA-{}/linear_regression_modeling/".format(cancer) + genes[i] + "_linear.csv")
        df_linregress_gene1 = df_linregress_gene[df_linregress_gene["r"] >= 0.3]
        num, _ = df_linregress_gene.shape
        if num != 0:
            gene2s = list(df_linregress_gene1["gene2"])
            slopes = list(df_linregress_gene1["slope"])
            slopes = np.array(slopes).reshape(-1, 1)
            intercepts = list(df_linregress_gene1["intercept"])
            intercepts = np.array(intercepts).reshape(-1, 1)
            mads = list(df_linregress_gene1["mad"])
            mads = np.array(mads).reshape(-1, 1)
            mads = np.where(mads == 0, 1e-9, mads)
            meds = list(df_linregress_gene1["med"])
            meds = np.array(meds).reshape(-1, 1)

            x = np.array(mrs_case.loc[genes[i], :].values)
            x = np.tile(x, (len(gene2s), 1))
            y = np.array(mrs_case.loc[gene2s, :].values)

            distances = (slopes * x) - y + intercepts
            distances = np.abs(distances) / np.sqrt(slopes * slopes + 1)


            zscores = (distances - meds) / (1.486 * mads)

            df_ppn = pd.DataFrame(columns=["gene1", "gene2"] + case_samples)
            df_ppn["gene1"] = [genes[i]] * len(gene2s)
            df_ppn["gene2"] = gene2s
            df_ppn.iloc[:, 2:] = zscores

            df_ppn.to_csv("./TCGA-{}/individualized_perturbed_edges/".format(cancer) + genes[i] + "_perturbed_edges.csv", index=False)


def multiprocess_estimate_individualized_perturbed_edges(cancer):

    os.mkdir("./TCGA-{}/individualized_perturbed_edges".format(cancer))

    # load normalized gene expression profiles of the cases samples
    mrs_case = pd.read_csv("./TCGA-{}/case_exp_normalized.csv".format(cancer))
    mrs_case = mrs_case.set_index("Symbol")
    mrs_case = mrs_case.transpose()
    print(mrs_case.head(5))

    # load community: gene-list dictionary
    with open("./TCGA-{}/comm_dict.pickle".format(cancer), "rb") as f:
        comm_dict = pickle.load(f)

    assigned_genes = []
    for module in comm_dict.keys():
        assigned_genes += comm_dict[module]

    mrs_case = mrs_case[assigned_genes]

    split_num = len(assigned_genes) // 1000

    genes_splited = []
    case_samples = sorted(list(mrs_case.index))
    for num in range(0, split_num):
        genes_splited.append(assigned_genes[1000 * num : 1000 * (num + 1)])

    genes_splited.append(assigned_genes[1000 * split_num :])

    # df_corrected = df_corrected.sort_index()
    mrs_case = mrs_case.transpose()
    # print(df_corrected.shape)

    procs = []
    for index, gene_list in enumerate(genes_splited):
        proc = Process(target=estimate_individualized_perturbed_edges, args=(mrs_case, case_samples, gene_list, cancer))
        procs.append(proc)
        proc.start()

    for proc in procs:
        proc.join()


def combine_perturbed_edges_per_community(genes, community_list, comm_dict, columns, cancer):
    for community in community_list:
        comm_genes = comm_dict[community]
        comm_genes = sorted(comm_genes)

        gene1s = []
        gene2s = []
        values = []

        for i in range(0, len(comm_genes)):
            if comm_genes[i] + "_perturbed_edges.csv" in genes:
                df_ppn = pd.read_csv("./TCGA-{}/individualized_perturbed_edges/".format(cancer) + comm_genes[i] + "_perturbed_edges.csv")

                df_ppn_filtered = df_ppn[df_ppn["gene2"].isin(comm_genes[i+1:])]
                row, _ = df_ppn_filtered.shape

                gene1 = []
                gene2 = []
                for symbol in list(df_ppn_filtered["gene2"]):
                    [temp1, temp2] = sorted([comm_genes[i], symbol])
                    gene1.append(temp1)
                    gene2.append(temp2)

                value = list(df_ppn_filtered.iloc[:, 2:].values)

                gene1s = gene1s +gene1
                gene2s = gene2s +gene2
                values = values +value

            else:
                pass

        df_comm_perturbed_edges = pd.DataFrame(columns=columns)
        df_comm_perturbed_edges["gene1"] = gene1s
        df_comm_perturbed_edges["gene2"] = gene2s
        df_comm_perturbed_edges.iloc[:, 2:] = values
        # print(df_module_ppn.head(5))
        df_comm_perturbed_edges.to_csv("./TCGA-{}/individualized_perturbed_edges_per_community/".format(cancer) + community + "/df_comm_perturbed_edges.csv", index=False)


def multiprocess_combine_perturbed_edges_per_module(cancer):

    os.mkdir("./TCGA-{}/individualized_perturbed_edges_per_community".format(cancer))
    genes = os.listdir("./TCGA-{}/individualized_perturbed_edges".format(cancer))

    with open("./TCGA-{}/comm_dict.pickle".format(cancer), "rb") as f:
        comm_dict = pickle.load(f)

    communities = sorted(list(comm_dict.keys()))

    for community in communities:
        path = "./TCGA-{}/individualized_perturbed_edges_per_community/".format(cancer) + community
        os.mkdir(path)

    comm_splited = []
    for num in range(0, len(communities)):
        comm_splited.append(communities[1 * num : 1 * (num + 1)])

    mrs_case = pd.read_csv("./TCGA-{}/case_exp_normalized.csv".format(cancer))
    mrs_case = mrs_case.set_index("Symbol")

    columns = ["gene1", "gene2"] + list(mrs_case.columns)

    procs = []
    for index, community_list in enumerate(comm_splited):
        proc = Process(target=combine_perturbed_edges_per_community,
                       args=(genes, community_list, comm_dict, columns, cancer))
        procs.append(proc)
        proc.start()

    for proc in procs:
        proc.join()


def quantify_community_cohesion_scores(communities, threshold, cancer):
    for community in communities:
        df_comm_perturbed_edges = pd.read_csv("./TCGA-{}/individualized_perturbed_edges_per_community/".format(cancer) + community + "/df_comm_perturbed_edges.csv")
        df_community_edges = pd.read_csv("./TCGA-{}/community_edges/".format(cancer) + community + "/df_edges_filtered.csv")

        num, _ = df_community_edges.shape
        if num != 0:
            subgraph = nx.Graph()
            for idx, row in df_community_edges.iterrows():
                subgraph.add_edge(row["gene1"], row["gene2"], weight=row["weight"])
            num, _ = df_community_edges.shape

            ne = global_efficiency(subgraph)

            df_merged = pd.merge(df_comm_perturbed_edges, df_community_edges,
                                 left_on=["gene1", "gene2"],
                                 right_on=["gene1", "gene2"],
                                 how="inner")
            samples = list(df_comm_perturbed_edges.columns[2:])

            scores = []

            for i in tqdm(range(0, len(samples))):
                df_comm_sample = df_merged[["gene1", "gene2", "weight", samples[i]]]
                df_comm_sample_filtered = df_comm_sample[abs(df_comm_sample[samples[i]]) > threshold]

                sample_edges_removed = [(gene1, gene2) for gene1, gene2 in zip(list(df_comm_sample_filtered["gene1"]), list(df_comm_sample_filtered["gene2"]))]
                num_sample, _ = df_comm_sample_filtered.shape

                if num_sample != 0:

                    # ign = individualized genetic network
                    ign = subgraph.copy()
                    ign.remove_edges_from(sample_edges_removed)
                    ign_ne = global_efficiency(ign)

                    for gene, original_weight in total_original_weights.items():
                        removed_weight = total_remove_weights.get(gene, 0)
                        connectivity_loss[gene] = removed_weight / original_weight if original_weight != 0 else 0

                    connectivity_loss_dict[samples[i]] = connectivity_loss

                else:
                    ign_ne = ne

                ccs = ign_ne / ne
                scores.append(ccs)

            ccs_dict = {}
            for i in range(0, len(samples)):
                ccs_dict[samples[i]] = scores[i]

            with open("./TCGA-{}/community_cohesion_scores/ccs_".format(cancer) + community + ".pickle", "wb") as f:
                pickle.dump(ccs_dict, f)

            with open("./TCGA-{}/connectivity_loss_scores/connectivity_loss_".format(cancer) + community + ".pickle", "wb") as f:
                pickle.dump(connectivity_loss_dict, f)

        else:
            samples = list(df_comm_perturbed_edges[2:])
            ccs_dict = {}
            connectivity_loss_dict = {}

            for i in range(0, len(samples)):
                ccs_dict[samples[i]] = 0
                connectivity_loss_dict[samples[i]] = {}

            with open("./TCGA-{}/community_cohesion_scores/ccs_".format(cancer) + community + ".pickle", "wb") as f:
                pickle.dump(ccs_dict, f)

            with open("./TCGA-{}/connectivity_loss_scores/connectivity_loss_".format(cancer) + community + ".pickle", "wb") as f:
                pickle.dump(connectivity_loss_dict, f)


def multiprocess_quantify_community_cohesion_scores_and_connectivity_loss_scores(cancer):

    os.mkdir("./TCGA-{}/community_cohesion_scores".format(cancer))
    os.mkdir("./TCGA-{}/connectivity_loss_scores".format(cancer))

    with open("./TCGA-{}/comm_dict.pickle".format(cancer), "rb") as f:
        comm_dict = pickle.load(f)

    communities = list(comm_dict.keys())

    community_splited = []
    for num in range(0, len(communities)):
        community_splited.append(communities[1 * (num) : 1* (num + 1)])

    procs = []
    for index, community_list in enumerate(community_splited):
        proc = Process(target=quantify_community_cohesion_scores_and_connectivity_loss_scores,
                       args=(community_list, 3.00, cancer))
        procs.append(proc)
        proc.start()

    for proc in procs:
        proc.join()


def make_df_ccs_and_df_cls(cancer):

    mrs_case = pd.read_csv("./TCGA-{}/case_exp_normalized.csv".format(cancer))
    mrs_case = mrs_case.set_index("Symbol")

    samples = list(mrs_case.columns)
    genes = list(mrs_case.index)

    with open("./TCGA-{}/comm_dict.pickle".format(cancer), "rb") as f:
        comm_dict = pickle.load(f)

    community_list = sorted(list(comm_dict.keys()))

    df_ccs = pd.DataFrame(index=community_list, columns=samples)
    df_cls = pd.DataFrame(index=genes, columns=samples)

    for i in tqdm(range(0, len(community_list))):
        if community_list[i] not in []:

            with open("./TCGA-{}/community_cohesion_scores/ccs_".format(cancer) + community_list[i] + ".pickle", "rb") as f:
                module_efficiency = pickle.load(f)

            with open("./TCGA-{}/connectivity_loss_scores/connectivity_loss_".format(cancer) + community_list[i] + ".pickle", "rb") as f:
                connectivity_loss = pickle.load(f)

            for sample in samples:
                df_ccs.loc[community_list[i], sample] = module_efficiency[sample]

                for gene in genes:
                    df_cls.loc[gene, sample] = connectivity_loss[sample].get(gene)

        df_ccs = df_ccs.astype('float')
        df_cls = df_cls.astype('float')
    df_ccs.to_csv("./TCGA-{}/community_cohesion_scores/df_ccs.csv".format(cancer), index=True)
    df_cls.to_csv("./TCGA-{}/connectivity_loss_scores/df_cls.csv".format(cancer), index=True)


if __name__ == '__main__':
    # -----------------------------------------------------
    # Checking for input from the command line:
    # -----------------------------------------------------
    #
    # [1] file providing the genes and the corresponding genetic community
    #     (gene, community, .csv format)
    #
    # [2] file providing the network edges and their weights
    #     (index: gene, column: gene)
    #
    # [3] file providing the gene expression profiles of control samples
    #     (index: sample, column: gene)
    #
    # [4] file providing the gene expression profiles of case samples (at least one sample)
    #     (index: sample, column: gene)
    #
    #

    os.chdir("")
    input_list = ["./df_modules.txt", "./df_tom_similarity_none.csv", "./exp_normal_normalized.txt", "./exp_data_normalized.txt"]
    comm_file, weight_file, control_file, case_file = check_input_style(input_list)
    cancer = ""
    make_community_dict(comm_file, cancer)
    make_edge_list_per_community(weight_file, cancer)
    convert_to_median_rank_scores(control_file, case_file, cancer)
    multiprocess_measure_linear_regression(cancer)
    multiprocess_estimate_individualized_perturbed_edges(cancer)
    multiprocess_combine_perturbed_edges_per_module(cancer)
    multiprocess_quantify_community_cohesion_scores_and_connectivity_loss_scores(cancer)
    make_df_ccs_and_df_cls(cancer)