import numpy as np
from sklearn import preprocessing
import sklearn.model_selection

def scale_minmax(data, axis=0):
    mi = np.expand_dims(data.min(axis=axis), axis=axis)
    ma = np.expand_dims(data.max(axis=axis), axis=axis)
    data_scaled = (data - mi) / (ma - mi)
    return data_scaled


def scale_gene_feature(data, scale_type=0):
    if scale_type == 0:
        data = scale_minmax(data, axis=0)
    elif scale_type == 1:
        data = scale_minmax(data, axis=1)
    elif scale_type == 2:
        min_max_scaler = preprocessing.MinMaxScaler()
        data = min_max_scaler.fit_transform(data)
    elif scale_type == 3:
        data = preprocessing.StandardScaler().fit_transform(data)
    if scale_type < 3:
        data = np.clip(a=data, a_min=1e-4, a_max=1 - 1e-4)
    return data


def get_gene_label(gene_name, driver_gene, non_driver, label_type=True):
    gene_name = gene_name[:, 0]
    driver_gene = driver_gene[:, 0]
    driver_label = np.asarray([g in driver_gene for g in gene_name])
    if label_type:
        non_driver = non_driver[:, 0]
        non_driver_label = np.asarray([g in non_driver for g in gene_name])
        all_labeled = driver_label + non_driver_label
        gene_label = np.concatenate([driver_label[:, np.newaxis], all_labeled[:, np.newaxis]], axis=1)
    else:
        gene_label = np.concatenate([driver_label[:, np.newaxis], np.ones_like(driver_label)[:, np.newaxis]], axis=1)
    return gene_label

def get_y_from_indices(y_label, y_mask, indices):
    assert (y_label.shape[0] == y_mask.shape[0])
    # construct y
    y_sub = np.zeros_like(y_label)
    y_sub[indices] = y_label[indices]
    # construct the mask
    m_sub = np.zeros_like(y_mask)
    m_sub[indices] = 1
    return np.expand_dims(y_sub, -1), np.expand_dims(m_sub, -1)


def cross_validation_sets(y, folds):
    y_label = y[:, 0]
    y_mask = y[:, 1]
    label_idx = np.where(y_mask)[0]
    # train set & test set
    kf_0 = sklearn.model_selection.StratifiedKFold(n_splits=4, shuffle=True)
    k_sets_0 = []
    splits_0 = kf_0.split(label_idx, y_label[label_idx])
    for train, test in splits_0:
        train_label_idx = label_idx[train]
        test_label_idx = label_idx[test]
        y_test, test_mask = get_y_from_indices(y_label, y_mask, test_label_idx)
        k_sets_0.append((train_label_idx, y_test, test_mask))
        break
    train_label_idx = k_sets_0[0][0]
    # cv_sets
    kf = sklearn.model_selection.StratifiedKFold(n_splits=folds, shuffle=True)
    k_sets = []
    splits = kf.split(train_label_idx, y_label[train_label_idx])
    for train, val in splits:
        train_idx = label_idx[train]
        val_idx = label_idx[val]
        # construct y and mask for train_set and test_set
        y_train, train_mask = get_y_from_indices(y_label, y_mask, train_idx)
        y_val, val_mask = get_y_from_indices(y_label, y_mask, val_idx)
        k_sets.append((y_train, y_val, train_mask, val_mask))
    return k_sets, k_sets_0[0][1], k_sets_0[0][2]
