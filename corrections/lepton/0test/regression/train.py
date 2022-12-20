import json
import matplotlib.pyplot as plt
import torch
import torch.nn as nn
import torch.nn.functional as F
import uproot
import pandas as pd
import ROOT
import numpy as np

class regressor(nn.Module):
    def __init__(self, n_inputs):
        super(regressor, self).__init__()
        self.output = nn.Sequential(
            nn.Linear(n_inputs, 20),
            nn.ELU(),
            nn.Linear(20, 50),
            nn.ELU(),
            nn.Linear(50,25),
            nn.Softmax(),
            nn.Linear(25,10),
            nn.ELU(),
            nn.Linear(10,2),
            nn.LeakyReLU()
        )

    def forward(self, x):
        out = self.output(x)

        return out

# TODO neue Idee
# Schritt 1: Histogrammierte Massen-Verteilung eingeben und zwei Werte finden so, dass die Verteilung damit skaliert/verschmiert mit der anderen uebereinstimmt -> finde mu, sigma
# Schritt 2: Herausfinden, wie sich die Veraenderung von m durch die Veraenderung von p1, p2 ergibt
#               ggf aus p1, p2 usw. mit mu, sigma jeweils 0, 0 die Verteilung von m lernen, dann mu, sigma dazugeben und die Verteilung von m_data lernen?

class kl_loss(nn.Module):
    def __init__(self):
        super(kl_loss, self).__init__()

    def forward(self, corr, mc, dt):
        pt_1 = torch.mul(corr[:,0], mc[:,0]) #TODO elementwise multiplication
        pt_2 = torch.mul(corr[:,1], mc[:,1])
        
        px_1 = torch.mul(pt_1, torch.cos(mc[:,4]))
        px_2 = torch.mul(pt_2, torch.cos(mc[:,5]))

        py_1 = torch.mul(pt_1, torch.sin(mc[:,4]))
        py_2 = torch.mul(pt_2, torch.sin(mc[:,5]))

        pz_1 = torch.mul(pt_1, torch.sinh(mc[:,6]))
        pz_2 = torch.mul(pt_2, torch.sinh(mc[:,7]))

        p1 = torch.sqrt(torch.mul(pt_1, pt_1) + torch.mul(pz_1, pz_1))
        p2 = torch.sqrt(torch.mul(pt_2, pt_2) + torch.mul(pz_2, pz_2))

        p1p2 = torch.mul(p1, p2)

        d_ang = torch.acos(torch.div(torch.mul(px_1, px_2)+torch.mul(py_1, py_2)+torch.mul(pz_1, pz_2), p1p2))

        mz_mc = torch.sqrt(torch.mul(2 * p1p2, 1-torch.cos(d_ang)))

        mz_dt = dt[:,8]

        hist_mc = torch.histc(mz_mc, bins=100, min=60, max=120)
        hist_dt = torch.histc(mz_dt, bins=100, min=60, max=120)

        return F.kl_div(hist_mc, hist_dt, reduction='mean', log_target=False)


class naive_loss(nn.Module):
    def __init__(self):
        super(naive_loss, self).__init__()

    def forward(self, corr, mc, dt):
        pt_1 = torch.mul(corr[:,0], mc[:,0]) #TODO elementwise multiplication
        pt_2 = torch.mul(corr[:,1], mc[:,1])

        #print(pt_1, mc[:,0])
        
        px_1 = torch.mul(pt_1, torch.cos(mc[:,4]))
        px_2 = torch.mul(pt_2, torch.cos(mc[:,5]))

        py_1 = torch.mul(pt_1, torch.sin(mc[:,4]))
        py_2 = torch.mul(pt_2, torch.sin(mc[:,5]))

        pz_1 = torch.mul(pt_1, torch.sinh(mc[:,6]))
        pz_2 = torch.mul(pt_2, torch.sinh(mc[:,7]))

        p1 = torch.sqrt(torch.mul(pt_1, pt_1) + torch.mul(pz_1, pz_1))
        p2 = torch.sqrt(torch.mul(pt_2, pt_2) + torch.mul(pz_2, pz_2))

        p1p2 = torch.mul(p1, p2)

        d_ang = torch.acos(torch.div(torch.mul(px_1, px_2)+torch.mul(py_1, py_2)+torch.mul(pz_1, pz_2), p1p2))

        mz_mc = torch.sqrt(torch.mul(2 * p1p2, 1-torch.cos(d_ang)))

        #print(corr[:,0])

        mz_true = 91.1876*torch.ones(corr[:,0].size())

        mz_diff = mz_mc - mz_true
        #print(mc[:,8],mz_mc)

        mz_diff_sq = torch.mul(mz_diff, mz_diff)

        MSE2 = torch.sum(mz_diff_sq)

        return MSE2


def plot(x='m_vis'):
    plt.hist(mc[x], bins=50, color='g', histtype='step', density=True)
    plt.hist(dt[x], bins=50, color='b', histtype='step', density=True)
    plt.show()
    plt.close()
    return

if __name__ == '__main__':
    device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
    quants = ["pt_1", "pt_2", "q_1", "q_2", "phi_1", "phi_2", "eta_1", "eta_2", "m_vis"]

    """
    ev_mc = uproot.open("/ceph/moh/CROWN_samples/Run3V02/ntuples/2022/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/mm//DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X_*0.root:ntuple")
    ev_dt = uproot.open('/ceph/moh/CROWN_samples/Run3V02/ntuples/2022/SingleMuon_Run2022C-PromptReco-v1/mm/SingleMuon_Run2022C-PromptReco-v1_*0.root:ntuple')
    print(ev_mc)
    mc = ev_mc.arrays(quants, library="pd")
    dt = ev_dt.arrays(quants, library="pd")
    """
    
    df_mc = ROOT.RDataFrame("ntuple", "/ceph/moh/CROWN_samples/Run3V02/ntuples/2022/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/mm//DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X_*0.root")
    df_dt = ROOT.RDataFrame("ntuple", '/ceph/moh/CROWN_samples/Run3V02/ntuples/2022/SingleMuon_Run2022C-PromptReco-v1/mm/SingleMuon_Run2022C-PromptReco-v1_*.root')
    
    mc = pd.DataFrame(df_mc.AsNumpy(columns = quants))
    dt = pd.DataFrame(df_dt.AsNumpy(columns = quants))

    #print(len(mc['pt_1']), len(dt['pt_1']))

    #plot()

    # train
    classifier = regressor(len(quants)).to(device)
    optimizer = torch.optim.Adam(classifier.parameters(), lr=0.01)
    criterion = naive_loss()

    n_epochs = 10
    batch_size = 5000

    n_batches = int(len(mc['pt_1'])/batch_size)

    # Keep track of the losses 
    train_losses = []
    val_losses = []

    for ep in range(n_epochs):

        for i in range(n_batches):
            
            # Reset gradient
            optimizer.zero_grad()
            
            i_start = i*batch_size
            i_stop  = (i+1)*batch_size
            
            # Convert x and y to proper objects for PyTorch
            x = torch.tensor(mc.to_numpy()[i_start: i_stop][:], dtype=torch.float).to(device)
            y = torch.tensor(dt.to_numpy()[i_start: i_stop][:], dtype=torch.float).to(device)

            #print(x)

            # Apply the network 
            net_out = classifier(x)
            #print(net_out)
                    
            # Calculate the loss function
            loss = criterion(net_out, x, y)
                    
            # Calculate the gradients
            loss.backward()
            
            # Update the weights
            optimizer.step()
            
            
        # Calculate predictions for the full training and validation sample
        corr = classifier(torch.tensor(mc.to_numpy(), dtype=torch.float).to(device))
        x_train = torch.tensor(mc.to_numpy(), dtype=torch.float).to(device)
        y_train = torch.tensor(dt.to_numpy(), dtype=torch.float).to(device)
        #y_pred_val = classifier(torch.tensor(dt.to_numpy(), dtype=torch.float).to(device)).detach().cpu().numpy().flatten()


        print(corr, x_train[:,0], x_train[:,8])
        # Calculate aver loss / example over the epoch
        #print(corr.size(), x_train.size(), y_train.size())
        train_loss = criterion(corr, x_train, y_train)
        #val_loss = criterion(torch.tensor(y_pred_val, dtype=torch.float).to(device), torch.tensor(y_val, dtype=torch.float).to(device))
        
        # print some information
        print("Epoch:", ep, "Train Loss:", train_loss.item())
        #print("Epoch:", ep, "Train Loss:", train_loss.item(),  "Test Loss:", val_loss.item())
        
        # and store the losses for later use
        train_losses.append(train_loss.item())
        #val_losses.append(val_loss.item())


