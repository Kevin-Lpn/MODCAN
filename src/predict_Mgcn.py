import os
import MgcnPre
import numpy as np
import MgcnIO
import torch
import utils
import gc
import argparse
import Mgcn

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'


def parse_args():
    parser = argparse.ArgumentParser(description='Load model and predict.')
    parser.add_argument('-ct', '--cancer', help='cancer type',
                        dest='cancer',
                        default='',
                        type=str
                        )
    parser.add_argument('-md', '--model_dir', help='model dir',
                        dest='model_dir',
                        default='',
                        type=str
                        )
    parser.add_argument('-cv', '--cv_idx', help='model with cross validation: i',
                        dest='cv_idx',
                        default=,
                        type=int
                        )
    args = parser.parse_args()
    return args


def run_Mgcn_model(model, gene_feature, hp_graph, gene_name, non_driver, output_dir, cgc):
    print("\n\tPredict...")
    gene_label = MgcnPre.get_gene_label(gene_name, cgc, non_driver)
    all_mask = np.ones(gene_label.shape[0], dtype=bool)[:, np.newaxis]
    model.eval()
    output, logits, loss, performance = Mgcn.predict_and_performance(model=model,
                                                                     feature=gene_feature,
                                                                     hp_graph=hp_graph,
                                                                     label=gene_label[:, [0]],
                                                                     mask=all_mask)
    print("\n\t\t************************************",
          "\tPred_acc= {:.6f}".format(performance[0]),
          "\tPred_auc= ", "{:.6f}".format(performance[1]),
          "\tPred_aupr= ", "{:.6f}".format(performance[2]),
          "\t************************************"
          )
    MgcnIO.save_predictions(output_dir, gene_name, logits)
    return performance


if __name__ == '__main__':
    args = parse_args()
    cancer_type = args.cancer
    model_dir = args.model_dir
    cv_idx = args.cv_idx

    print("\n")
    print("\t\t************************************************************************")
    print("\t\t***                                                                  ***")
    print("\t\t************************************************************************")
    params_dir = "../data/MGCN/training_{}/{}".format(cancer_type, model_dir)
    args_dict, input_file = MgcnIO.load_hyper_params(params_dir)
    if cv_idx != -1:
        args_dict['cv_idx'] = cv_idx
    print("\nHyper parameters and datasets:")
    print(args_dict)

    input_data_path = args_dict['data'][9:]
    print("\nLoad data from: {}...".format(input_data_path))
    data = MgcnIO.load_hdf_data(input_data_path, network_name='Network')
    (gene_feature, feature_name, gene_name, exp_data, met_data, sample_id, ppi, non_driver,
     mut_matrix, exp_normal, met_normal, sample_cluster, cgc) = data

    output_dir = MgcnIO.create_model_dir(cancer_type, is_training=False)
    print("\nCreate output dir: {}".format(output_dir))

    print("\nData preprocessing...")
    (gene_feature, gene_name, exp_data, met_data, sample_id,
     mut_matrix, exp_normal, met_normal, sample_cluster) = utils.filter_data(gene_feature=gene_feature,
                                                                 gene_name=gene_name,
                                                                 exp_data=exp_data,
                                                                 met_data=met_data,
                                                                 sample_id=sample_id,
                                                                 ppi=ppi,
                                                                 mut_matrix=mut_matrix,
                                                                 exp_normal=exp_normal,
                                                                 met_normal=met_normal,
                                                                 sample_cluster=sample_cluster
                                                                 )

    print("\n\tIntegrate and scale gene feature...")
    gene_feature = MgcnPre.scale_gene_feature(gene_feature, scale_type=0)

    del (data, feature_name, sample_id, mut_matrix)
    gc.collect()

    print("\nCreate difference co-association network...")
    diff_co_data_exp = utils.calculate_co_association(exp_data, exp_normal, gene_name, sample_cluster, args_dict['net_split'])
    diff_co_data_met = utils.calculate_co_association(met_data, met_normal, gene_name, sample_cluster, args_dict['net_split'])

    diff_co_data = list()
    for exp_diff, met_diff in zip(diff_co_data_exp, diff_co_data_met):
        combined_diff = 0.7 * exp_diff + 0.3 * met_diff
        diff_co_data.append(combined_diff)

    print("\nCreate adjacent matrix of PPI network...")
    ppi_network = utils.ppi_limitation(ppi, gene_name, weight_type=args_dict['weight_type'])

    print("\nCreate differential co-association network...")
    co_data = utils.construct_association_network(diff_co_data, ppi_network)

    print("\n\tCreate Weighted Hypergraph...")
    hp_graph = utils.get_hypergraph_matrix(gene_name=gene_name, co_data_net=co_data)

    del (diff_co_data_exp, diff_co_data_met, diff_co_data, ppi_network, co_data)
    gc.collect()

    print("\nLoad model...")
    model = Mgcn.MGCN(n_input=gene_feature.shape[1],
                      n_edge=args_dict['n_edge'],
                      hidden_dims=args_dict['hidden_dims'],
                      learning_rate=args_dict['lr'],
                      weight_decay=args_dict['decay'],
                      dropout=args_dict['dropout'],
                      pos_loss_multiplier=args_dict['loss_mul'],
                      logging=True
                      )

    model_path = "{}/cv_{}/saved_model_cv_{}.pkl".format(params_dir, args_dict['cv_idx'], args_dict['cv_idx'])
    model.load_state_dict(torch.load(model_path, ))
    print("\n\tLoad model from '{}'".format(model_path, weights_only=True))

    print("\nRunning...")
    performance = run_Mgcn_model(model, gene_feature, hp_graph, gene_name, non_driver, output_dir, cgc)
    print("\nFinished!")
