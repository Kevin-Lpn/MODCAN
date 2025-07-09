import numpy as np
import os
import h5py
from datetime import datetime



def btype_to_str(data):
    n_row = data.shape[0]
    n_col = data.shape[1]
    str_data = []
    for i in range(n_row):
        for j in range(n_col):
            str_data.append(data[i, j].decode())
    str_data = np.array(str_data).reshape(n_row, n_col)
    return str_data

def load_hdf_data(path, network_name='Network'):
    with h5py.File(path, 'r') as f:
        gene_feature = f['Gene_Feature'][:]
        feature_name = f['Feature_Name'][:]
        gene_name = f['Gene_Name'][:]
        exp_data = f['Exp_Data'][:]
        met_data = f['Met_Data'][:]
        sample_id = f['Sample_ID'][:]
        ppi = f[network_name][:]
        non_driver = f['Non_Driver_Gene'][:]
        mut_matrix = f['Mutation_Matrix'][:]
        exp_normal = f['Exp_Normal'][:]
        met_normal = f['Met_Normal'][:]
        sample_cluster = f['Sample_Cluster'][:]
        cgc = f['CGC'][:]

    feature_name = btype_to_str(feature_name)
    gene_name = btype_to_str(gene_name)
    sample_id = btype_to_str(sample_id)
    ppi = btype_to_str(ppi)
    non_driver = btype_to_str(non_driver)
    non_driver = non_driver[:, np.newaxis]
    sample_cluster = btype_to_str(sample_cluster)
    cgc = btype_to_str(cgc)
    cgc = cgc[:, np.newaxis]

    return (gene_feature, feature_name, gene_name, exp_data, met_data, sample_id, ppi, non_driver,
            mut_matrix, exp_normal, met_normal, sample_cluster, cgc)


def write_hyper_params(args_dict, input_file, file_name):
    with open(file_name, 'w') as f:
        for arg in args_dict:
            f.write('{}\t{}\n'.format(arg, args_dict[arg]))
        f.write('{}\n'.format(input_file))
    print("Hyper-Parameters saved to {}".format(file_name))


def str_to_num(s):
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s


def load_hyper_params(model_dir):
    file_name = os.path.join(model_dir, 'hyper_params.txt')
    input_file = None
    with open(file_name, 'r') as f:
        args = {}
        for line in f:
            if '\t' in line:
                key, value = line.split('\t')
                if value.startswith('['):
                    f = lambda x: "".join(c for c in x if c not in ['\"', '\'', ' ', '\n', '[', ']'])
                    l = [int(f(i)) for i in value.split(',')]
                    args[key.strip()] = l
                else:
                    args[key.strip()] = str_to_num(value.strip())
            else:
                input_file = line.strip()
    print("Hyper-Parameters read from {}".format(file_name))
    return args, input_file


def create_model_dir(cancer_type, is_training=True):
    if is_training:
        root_dir = '../data/MGCN/training_{}'.format(cancer_type)
    else:
        root_dir = '../prediction/{}'.format(cancer_type)
    if not os.path.isdir(root_dir):
        os.mkdir(root_dir)
        print("Created Output Subdir")
    date_string = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    os.mkdir(os.path.join(root_dir, date_string))
    model_path = os.path.join(root_dir, date_string)
    return model_path


def save_predictions(output_dir, node_names, predictions):
    with open(os.path.join(output_dir, 'predictions.tsv'), 'w') as f:
        f.write('Name\tProb_neg\tProb_pos\n')
        for pred_idx in range(predictions.shape[0]):
            f.write('{}\t{}\t{}\n'.format(node_names[pred_idx, 0],
                                          predictions[pred_idx, 0],
                                          predictions[pred_idx, 1])
                    )


def write_train_test_sets(out_dir, y_train, train_mask, y_val, val_mask, y_test, test_mask):
    np.save(os.path.join(out_dir, 'y_train.npy'), y_train)
    np.save(os.path.join(out_dir, 'train_mask.npy'), train_mask)
    np.save(os.path.join(out_dir, 'y_val.npy'), y_val)
    np.save(os.path.join(out_dir, 'val_mask.npy'), val_mask)
    np.save(os.path.join(out_dir, 'y_test.npy'), y_test)
    np.save(os.path.join(out_dir, 'test_mask.npy'), test_mask)


def read_train_test_sets(in_dir):
    y_train = np.load(os.path.join(in_dir, 'y_train.npy'))
    train_mask = np.load(os.path.join(in_dir, 'train_mask.npy'))
    y_val = np.load(os.path.join(in_dir, 'y_val.npy'))
    val_mask = np.load(os.path.join(in_dir, 'val_mask.npy'))
    y_test = np.load(os.path.join(in_dir, 'y_test.npy'))
    test_mask = np.load(os.path.join(in_dir, 'test_mask.npy'))
    return y_train, train_mask, y_val, val_mask, y_test, test_mask


def save_performance_measures(output_dir, performance_measures):
    cv_aupr = list()
    with open(os.path.join(output_dir, 'performance_measures.tsv'), 'w') as f:
        for idx in range(len(performance_measures)):
            performance = performance_measures[idx]
            cv_aupr.append(performance[0][2])
            f.write('cv_{}\tloss= {}\n'.format(idx, performance[2]))
            f.write('Test_Acc\tTest_AUC\tTest_AUPR\tAcc\tAUC\tAUPR\n')
            f.write('{}\t{}\t{}\t{}\t{}\t{}\n\n'.format(performance[0][0],
                                                        performance[0][1],
                                                        performance[0][2],
                                                        performance[1][0],
                                                        performance[1][1],
                                                        performance[1][2])
                    )
        return cv_aupr.index(max(cv_aupr))