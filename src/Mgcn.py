import torch
import torch.nn as nn
from torch.nn import Linear
import torch.nn.functional as F
from torch.autograd import Variable
from torch.optim import Adam
import torchmetrics
import os
import copy



class MGCN(nn.Module):
    def __init__(self, n_input, n_edge, hidden_dims, learning_rate=0.1, weight_decay=5e-4, dropout=0.5,
                 pos_loss_multiplier=1, logging=True):
        super(MGCN, self).__init__()
        self.n_input = n_input
        self.n_edge = n_edge
        self.hidden_dims = hidden_dims
        self.layer_num = len(hidden_dims)
        self.learning_rate = learning_rate
        self.weight_decay = weight_decay
        self.pos_loss_multiplier = pos_loss_multiplier
        self.logging = logging

        # construct convolution layers
        self.layer_0 = Linear(self.n_input * self.n_edge, self.hidden_dims[0], bias=True)
        self.layer_1 = Linear(self.hidden_dims[0] * self.n_edge, self.hidden_dims[1], bias=True)
        self.layer_2 = Linear(self.hidden_dims[1], 2, bias=True)
        self.dropout = nn.Dropout(p=dropout)


    def calculate_slice_tensor(self, data):
        slice_num = data.shape[2]
        mgcn_x = data[:, :, 0]
        for s in range(1, slice_num):
            mgcn_x = torch.cat((mgcn_x, data[:, :, s]), dim=1)
        return mgcn_x

    def forward(self, feature, hp_graph):
        mgcn_x0 = self.calculate_slice_tensor(torch.einsum("nij,ik->nkj", hp_graph, feature))
        mgcn_x0 = self.dropout(F.relu(self.layer_0(mgcn_x0)))
        mgcn_x1 = self.calculate_slice_tensor(torch.einsum("nij,ik->nkj", hp_graph, mgcn_x0))
        mgcn_x1 = F.relu(self.layer_1(mgcn_x1))
        output = self.layer_2(mgcn_x1)
        return output, F.softmax(output, dim=1)


def masked_cross_entropy_loss(score, label, mask, weight):
    """Cross-entropy loss with masks."""
    label = torch.Tensor(label[:, 0]).long()
    pos_weight = torch.ones(score.shape[0]) + (weight - 1) * label
    loss_func_none = nn.CrossEntropyLoss(reduction="none")
    loss = pos_weight * loss_func_none(score, label)
    mask = torch.Tensor(mask[:, 0])
    mask /= torch.mean(mask)
    loss = torch.matmul(loss, mask) / score.shape[0]
    return loss

def get_performance_metrics(logits, label, mask):
    logits = logits[mask[:, 0], :]
    label = torch.Tensor(label[mask[:, 0], 0]).long()
    # acc auc aupr
    accuracy = torchmetrics.Accuracy(task="multiclass", num_classes=2, top_k=1, average='micro')
    acc = accuracy(logits, label)
    _AUC = torchmetrics.AUROC(task="binary", num_classes=2)
    auc = _AUC(logits[:, 1], label)
    average_precision = torchmetrics.AveragePrecision(task="binary", num_classes=2)
    aupr = average_precision(logits[:, 1], label)
    return acc, auc, aupr


def fit_model(model, feature, hp_graph, y_train, train_mask,
              y_val, val_mask, epochs, model_dir, save_model=True):
    hp_graph = torch.Tensor(hp_graph)
    log_interval = 10
    alpha = 0.
    model_alpha = copy.deepcopy(model)

    optimizer = Adam(model.parameters(), lr=model.learning_rate, weight_decay=model.weight_decay)
    for epoch in range(1, epochs + 1):
        model.train()
        x_tensor = Variable(torch.Tensor(feature))
        output, logits = model(x_tensor, hp_graph)
        loss = masked_cross_entropy_loss(score=output,
                                         label=y_train,
                                         mask=train_mask,
                                         weight=model.pos_loss_multiplier)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        val_performance = get_performance_metrics(logits, y_val, val_mask)
        if val_performance[2] > alpha and epoch >= 500:
            alpha = val_performance[2]
            model_alpha = copy.deepcopy(model)
            print("\t***************************",
                  "\tTrain Epoch: {}\ttrain_loss: {:.6f}".format(epoch, loss.item()),
                  "\tval_acc= {:.6f}".format(val_performance[0]),
                  "\tval_auc= ", "{:.6f}  ".format(val_performance[1]),
                  "\tval_aupr= ", "{:.6f}  ".format(val_performance[2]),
                  "\t***************************",
                  )
        if model.logging:
            if epoch % log_interval == 0 or epoch == 1:
                print("\t\tTrain Epoch: {}\ttrain_loss: {:.6f}".format(epoch, loss.item()),
                      "\tval_acc= {:.6f}".format(val_performance[0]),
                      "\tval_auc= ", "{:.6f}  ".format(val_performance[1]),
                      "\tval_aupr= ", "{:.6f}  ".format(val_performance[2]),
                      )
            if epoch == epochs:
                print('\n\t\tTraining completed')

        if save_model:
            print("\n\t\tSave model to {}".format(model_dir))
            model_save_path = os.path.join(model_dir, 'saved_model_cv_{}.pkl'.format(model_dir[-1]))
            torch.save(model_alpha.state_dict(), model_save_path)

    return model_alpha


def predict_and_performance(model, feature, hp_graph, label, mask):
    x_tensor = Variable(torch.Tensor(feature))
    hp_graph = torch.Tensor(hp_graph)
    model.eval()
    output, logits = model(x_tensor, hp_graph)
    loss = masked_cross_entropy_loss(score=output,
                                     label=label,
                                     mask=mask,
                                     weight=model.pos_loss_multiplier)
    performance_metrics = get_performance_metrics(logits,label,mask)
    return output, logits, loss, performance_metrics
