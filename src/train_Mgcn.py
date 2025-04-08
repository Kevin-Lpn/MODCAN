import argparse
import os
import sys
import numpy as np

import Mgcn
import MgcnIO
import utils
import MgcnPre
import gc






os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

def parse_args():
    parser = argparse.ArgumentParser(description='Train IHGC with cross_validation and save model to file')
    parser.add_argument('-e', '--epochs', help='number of epochs',
                        dest='epochs',
                        default=1000,
                        type=int
                        )
    parser.add_argument('-lr', '--learning_rate', help='Learning Rate',
                        dest='lr',
                        default=.001,
                        type=float
                        )
    parser.add_argument('-hd', '--hidden_dims',
                        help='Hidden Dimensions (number of filters per layer).',
                        nargs='+',
                        dest='hidden_dims',
                        default=[64, 128])
    parser.add_argument('-lm', '--loss_mul',
                        help='Number of times, false negatives are weighted higher than false positives',
                        dest='loss_mul',
                        default=2,
                        type=float
                        )
    parser.add_argument('-wd', '--weight_decay', help='Weight Decay',
                        dest='decay',
                        default=5e-4,
                        type=float
                        )
    parser.add_argument('-do', '--dropout', help='Dropout Percentage',
                        dest='dropout',
                        default=.5,
                        type=float
                        )
    parser.add_argument('-wt', '--weight_type', help='Weight type of Hypergraph',
                        dest='weight_type',
                        default=2,
                        type=int
                        )
    parser.add_argument('-ns', '--network_split', help='[0, 1)',
                        dest='net_split',
                        default=0.2,
                        type=float
                        )
    parser.add_argument('-d', '--data', help='Path to HDF5 container with data',
                        dest='data',
                        default='',
                        type=str,
                        # required=True
                        )
    parser.add_argument('-cv', '--cv_runs', help='Number of cross validation runs',
                        dest='cv_runs',
                        default=10,
                        type=int
                        )
    parser.add_argument('-dr', '--driver', help='cgc',
                        dest='driver',
                        default='cgc',
                        type=str
                        )
    parser.add_argument('-xr', '--xr', help='xiao rong',
                        dest='xr',
                        default='',
                        type=str
                        )
    args = parser.parse_args()
    return args


def single_cv_run(gene_feature, hp_graph, y_train, train_mask, y_val, val_mask, y_test, test_mask,
                  args_dict, model_dir, cv_run):
    print("\n\tConstruct model...")
    hidden_dims = [int(x) for x in args_dict['hidden_dims']]
    model = Mgcn.MGCN(n_input=gene_feature.shape[1],
                      n_edge=hp_graph.shape[2],
                      hidden_dims=hidden_dims,
                      learning_rate=args_dict['lr'],
                      weight_decay=args_dict['decay'],
                      dropout=args_dict['dropout'],
                      pos_loss_multiplier=args_dict['loss_mul'],
                      logging=True
                      )

    # fit the model
    print('\n\tFit model for cv_{}...\n'.format(cv_run))
    model = Mgcn.fit_model(model=model,
                           feature=gene_feature,
                           hp_graph=hp_graph,
                           y_train=y_train,
                           train_mask=train_mask,
                           y_val=y_val,
                           val_mask=val_mask,
                           epochs=args_dict['epochs'],
                           model_dir=model_dir
                           )
    _, _, test_loss, test_performance = Mgcn.predict_and_performance(model=model,
                                                                     feature=gene_feature,
                                                                     hp_graph=hp_graph,
                                                                     label=y_test,
                                                                     mask=test_mask)
    print("\n\t\t************************************",
          "\ttest_acc= {:.6f}".format(test_performance[0]),
          "\ttest_auc= {:.6f}".format(test_performance[1]),
          "\ttest_aupr= {:.6f}".format(test_performance[2]),
          "\t************************************"
          )
    MgcnIO.write_train_test_sets(model_dir, y_train, train_mask, y_val, val_mask, y_test, test_mask)
    return model, test_performance


def run_all_cvs(gene_feature, exp_data, met_data, exp_normal, met_normal, gene_name, sample_cluster, ppi, non_driver, args_dict, output_dir, cgc):
    print("\nCreate difference co-expression network...")
    diff_co_data_exp = utils.calculate_co_association(exp_data, exp_normal, gene_name, sample_cluster, args_dict['net_split'])
    diff_co_data_met = utils.calculate_co_association(met_data, met_normal, gene_name, sample_cluster, args_dict['net_split'])

    print("\nCreate adjacent matrix of PPI network...")
    ppi_network = utils.ppi_limitation(ppi, gene_name, weight_type=args_dict['weight_type'])

    print("\nCreate differential co-expression network...")
    co_data_exp_net = utils.construct_association_network(diff_co_data_exp, ppi_network)
    co_data_met_net = utils.construct_association_network(diff_co_data_met, ppi_network)

    print("\n\tCreate Weighted Hypergraph...")
    hp_graph_exp = utils.get_hypergraph_matrix(gene_name=gene_name, co_data_net=co_data_exp_net)
    hp_graph_met = utils.get_hypergraph_matrix(gene_name=gene_name, co_data_net=co_data_met_net)
    hp_graph = np.concatenate((hp_graph_exp, hp_graph_met), axis=2)

    del (diff_co_data_exp, diff_co_data_met, ppi_network, co_data_exp_net, co_data_met_net, hp_graph_exp, hp_graph_met)
    gc.collect()

    print("\nGet gene label and label mask...")
    gene_label = MgcnPre.get_gene_label(gene_name, cgc, non_driver)
    k_sets, y_test, test_mask = MgcnPre.cross_validation_sets(gene_label, folds=args_dict['cv_runs'])

    print("\nModel construction with cross validation...")
    performance_measures = []
    for cv_run in range(args_dict['cv_runs']):
        model_dir = os.path.join(output_dir, 'cv_{}'.format(cv_run))
        os.mkdir(model_dir)
        y_train, y_val, train_mask, val_mask = k_sets[cv_run]
        model, test_performance = single_cv_run(gene_feature, hp_graph, y_train, train_mask, y_val, val_mask,
                                                y_test, test_mask, args_dict, model_dir, cv_run)

        print("\n\tPrediction of cv_{}".format(cv_run))
        all_mask = np.ones_like(test_mask, dtype=bool)
        output, logits, loss, performance = Mgcn.predict_and_performance(model=model,
                                                                         feature=gene_feature,
                                                                         hp_graph=hp_graph,
                                                                         label=gene_label[:, [0]],
                                                                         mask=all_mask)

        print("\n\t************************************",
              "\tPred_acc= {:.6f}".format(performance[0]),
              "\tPred_auc= ", "{:.6f}".format(performance[1]),
              "\tPred_aupr= ", "{:.6f}".format(performance[2]),
              "\t************************************"
              )
        performance_measures.append((test_performance, performance, loss))
        MgcnIO.save_predictions(model_dir, gene_name, logits)

        del (model, output, logits, test_performance, performance, loss)
        gc.collect()

    print("\nSave performance_measures...")
    cv_idx = MgcnIO.save_performance_measures(output_dir, performance_measures)

    print("\nSave hyper parameters...")
    data_rel_to_model = os.path.relpath(args_dict['data'], output_dir)
    args_dict['data'] = data_rel_to_model
    args_dict['n_edge'] = hp_graph.shape[2]
    args_dict['cv_idx'] = cv_idx
    MgcnIO.write_hyper_params(args_dict, args_dict['data'], os.path.join(output_dir, 'hyper_params.txt'))

    return performance_measures


if __name__ == '__main__':
    args = parse_args()
    if not args.data.endswith('.h5'):
        print("\nData is not a hdf5 container. Exit now.")
        sys.exit(-1)
    args_dict = vars(args)

    print("\n")
    print("\t\t************************************************************************")
    print("\t\t***                                                                  ***")
    print("\t\t************************************************************************")
    print("\nHyper parameters and datasets:")
    print(args_dict)
    input_data_path = args.data

    print("\nLoad data from: {}...".format(input_data_path))
    data = MgcnIO.load_hdf_data(input_data_path, network_name='Network')
    (gene_feature, feature_name, gene_name, exp_data, met_data, sample_id, ppi, non_driver,
     mut_matrix, exp_normal, met_normal, sample_cluster, cgc) = data
    cancer_type = input_data_path[11:20]

    output_dir = MgcnIO.create_model_dir(cancer_type, is_training=True)
    print("\nCreate output dir: {}".format(output_dir))

    print("\nPre-processing...")
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

    print("\nModel construction and prediction performance...")
    performance_measures = run_all_cvs(gene_feature, exp_data, met_data, exp_normal, met_normal, gene_name, sample_cluster, ppi,
                                       non_driver, args_dict, output_dir, cgc)
    print('\nFinished!')




