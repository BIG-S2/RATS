import pandas as pd
import torch
import numpy as np
import torch.nn as nn
import torch.optim as optim
import os
import argparse

batch_size = 128
epochs = 3000
lr = 0.00005
use_cuda = torch.cuda.is_available()
device = torch.device("cuda" if use_cuda else "cpu")
print(
    "batch size {}\nepochs {}\nAdam lr {} \n using device {}\n".format(
        batch_size, epochs, lr, device.type
    )
)

parser = argparse.ArgumentParser(description='gRNA prediction model')
parser.add_argument('--seed', help='random seed used to generate cohort')
args = parser.parse_args()
seed = 776
seed1 = args.seed
print("Seed1: "+str(seed1))
data_dir = ""

dat = pd.read_csv(os.path.join(data_dir, "UKB_p12_select_hippo_resid.tsv"), sep="\t")
dat_mat = dat.iloc[:,2:].to_numpy(dtype = np.float32)
#df_reshape = dat_mat.reshape((dat_mat.shape[0],1,100,150), order = 'F')
df_reshape = dat_mat

pmax = np.max(df_reshape)
pmin = np.min(df_reshape)
train_norm = (df_reshape - pmin) / (pmax - pmin)
X = torch.tensor(train_norm, dtype=torch.float32).to(device)
train_loader = torch.utils.data.DataLoader(X, batch_size=batch_size, shuffle=True)


class Encoder(nn.Module):
    def __init__(self, image_size, hidden_size, latent_size):
        super(Encoder, self).__init__()
        self.image_size = image_size
        self.hidden_size = hidden_size
        self.latent_size = latent_size
        self.main = nn.Sequential(
            nn.Linear(self.image_size, 8 * hidden_size, bias=False),
            nn.LeakyReLU(0.2, True),
            nn.Linear(8 * hidden_size, 4 * hidden_size, bias=False),
            nn.LeakyReLU(0.2, True),
            nn.Linear(4 * hidden_size, hidden_size, bias=False),
            nn.LeakyReLU(0.2, True),
            nn.Linear(hidden_size, self.latent_size),
        )
    
    def forward(self, input):
        return self.main(input.view(-1, self.image_size))


class Decoder(nn.Module):
    def __init__(self, image_size, hidden_size, latent_size):
        super(Decoder, self).__init__()
        self.image_size = image_size
        self.hidden_size = hidden_size
        self.latent_size = latent_size
        self.main = nn.Sequential(
            nn.Linear(latent_size, hidden_size, bias=False),
            nn.ReLU(True),
            nn.Linear(hidden_size,  4* hidden_size, bias=False),
            nn.ReLU(True),
            nn.Linear(4 * hidden_size, 8 * hidden_size, bias=False),
            nn.ReLU(True),
            nn.Linear(8 * hidden_size, self.image_size, bias=False),
            nn.Sigmoid()
        )
    
    def forward(self, input):
        return self.main(input)

class Autoencoder(nn.Module):
    def __init__(self, image_size, hidden_size, latent_size, device):
        super(Autoencoder, self).__init__()
        self.image_size = image_size
        self.hidden_size = hidden_size
        self.latent_size = latent_size
        self.device = device
        self.encoder = Encoder(image_size, hidden_size, latent_size)
        self.decoder = Decoder(image_size, hidden_size, latent_size)
    
    def forward(self, input):
        return(self.decoder(self.encoder(input)))


model = Autoencoder(image_size=100*150, latent_size=128, hidden_size=512, device=device).to(device)
optimizer = optim.Adam(model.parameters(), lr=lr)
lossfunc = nn.MSELoss().to(device)

iter = 0
for epoch in range(epochs):
    # Training
    print('Epoch: ' + str(epoch))
    running_loss = 0.0
    for i, batch in enumerate(train_loader, 0):
        optimizer.zero_grad()
        reconstruct = model(batch)
        loss = lossfunc(batch, reconstruct)
        loss.backward()
        optimizer.step()
        iter += 1
        running_loss += loss.item()
        
        if iter % 25 == 24:    # print every 25 mini-batches
            print('[%d, %5d] loss: %.6f' %
                  (epoch + 1, iter + 1, running_loss / 25))
            running_loss = 0.0
            iter = 0

ckptPATH = os.path.join(data_dir, "UKB_p12_select_hippo_AE.pth")
torch.save(model.state_dict(), ckptPATH)

encoded = model.encoder(X)
dat_encoded = pd.DataFrame(encoded.detach().to('cpu').numpy())
dat_encoded.columns = ["Latent" + str(i+1) for i in range(dat_encoded.shape[1])]
dat_encoded['PTID'] = dat['PTID']
dat_encoded['AD_proxy'] = dat['AD_proxy']
cols = list(dat_encoded.columns)
cols = cols[-2:] + cols[:-2]
dat_encoded = dat_encoded[cols]
dat_encoded.to_csv(os.path.join(data_dir,"UKB_p12_select_hippo_encoded.tsv"), sep="\t", index=False)
