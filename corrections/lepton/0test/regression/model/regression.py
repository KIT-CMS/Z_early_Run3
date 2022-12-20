from itertools import zip_longest
from typing import Any, Callable, Tuple, Union
import uproot

import numpy as np
#import numpy.typing as npt
import torch as t
from torch.nn.modules.loss import _Loss

from copy import deepcopy
import os


def custom_histogram(
    ctx: t.autograd.Function,
    input: t.Tensor,
    bins: t.Tensor,
    mask: Union[Any, t.Tensor] = None,
    weights: Union[Any, t.Tensor] = None,
) -> t.Tensor:
    mask = mask.to(t.double) if mask is not None else t.ones_like(input).to(t.double)
    weights = weights if weights is not None else t.ones_like(input).to(t.double)

    left_edge, right_edge = bins[:-1], bins[1:]  # since no t.histogram is available on cuda devices...

    bool_mask = (((input * mask)[..., None] > left_edge) & ((input * mask)[..., None] <= right_edge)).T
    bool_mask = t.einsum("ji...->...ji", bool_mask)

    histogram = t.einsum(
        "ji, i, i -> j",
        bool_mask.to(t.double),
        weights,
        mask,
    )

    ctx.save_for_backward(input, mask, bins, weights, bool_mask)
    return histogram




#def cutted_gaussians_gradient(ctx: t.autograd.Function, output: t.Tensor) -> tuple[t.Tensor, Any, Any, Any]:
def cutted_gaussians_gradient(ctx: t.autograd.Function, output: t.Tensor):
    input, mask, bins, weights, bool_mask = ctx.saved_tensors

    sigma = t.diff(bins) / 2
    mu = bins[1:] - (t.diff(bins) / 2)
    input_shape, hist_shape = input.shape, mu.shape

    input = t.einsum("i, j -> ji", input, t.ones(hist_shape).to(input.device))
    mu = t.einsum("i, j -> ji", t.ones(input_shape).to(input.device), mu)
    sigma = t.einsum("i, j -> ji", t.ones(input_shape).to(input.device), sigma)

    x = (input - mu) * bool_mask
    _f = -((x) / (sigma ** 2)) * (-0.5 * ((x) / (sigma)) ** 2).exp()

    to_pass = t.einsum(
        "ji , j, i, i -> i",
        _f.to(t.double),
        output.to(t.double),
        mask.to(t.double),
        weights.to(t.double),
    )
    return to_pass, None, None, None

#def continous_gaussians_gradient(ctx: t.autograd.Function, output: t.Tensor) -> tuple[t.Tensor, Any, Any, Any]:
def continous_gaussians_gradient(ctx: t.autograd.Function, output: t.Tensor):
    input, mask, bins, weights, bool_mask = ctx.saved_tensors

    sigma = t.diff(bins) / 2
    mu = bins[1:] - (t.diff(bins) / 2)
    input_shape, hist_shape = input.shape, mu.shape

    input = t.einsum("i,j -> ji", input, t.ones(hist_shape).to(input.device))
    mu = t.einsum("i,j -> ji", t.ones(input_shape).to(input.device), mu)
    sigma = t.einsum("i,j -> ji", t.ones(input_shape).to(input.device), sigma)

    x = (input - mu) * t.ones_like(bool_mask)
    _f = -((x) / (sigma ** 2)) * (-0.5 * ((x) / (sigma)) ** 2).exp()

    to_pass = t.einsum(
        "ji , j, i, i -> i",
        _f.to(t.double),
        output.to(t.double),
        mask.to(t.double),
        weights.to(t.double),
    )
    return to_pass, None, None, None





class CustomHistGradient(t.autograd.Function):
    pass

setattr(CustomHistGradient, "forward", staticmethod(custom_histogram))
setattr(CustomHistGradient, "backward", staticmethod(cutted_gaussians_gradient))



class Chi2HistogrammDiffLoss(_Loss):
    def __init__(
        self,
        device: str,
        #bins: Union[int, list, npt.NDArray] = 8,
        bins: Union[int, list] = 8,
        size_average: Union[bool, Any] = None,
        reduce: Union[str, Any] = None,
        reduction: str = "mean",
    ) -> None:
        super().__init__(size_average, reduce, reduction)

        self.device = device

        if isinstance(bins, int):
            self.bins = t.linspace(0, 1, bins + 1).to(self.device).to(t.float32)
        elif isinstance(bins, (list, np.ndarray)) and all(isinstance(it, float) for it in bins):
            self.bins = t.tensor(bins).to(self.device).to(t.float32)

        self.hist_ = CustomHistGradient.apply

    def histograms(
        self,
        input: t.Tensor,
        weights: t.Tensor,
        target: Union[t.Tensor, Any] = None,
        shifted_input: t.Tensor = None,
        ) -> Tuple[t.Tensor, t.Tensor]:

        nominal_data = self.hist_(
            input,
            self.bins,
            None,
            None,
        )

        shifted_mc = self.hist_(
            shifted_input,
            self.bins,
            None,
            weights
        )

        return nominal_data, shifted_mc


    def forward(
        self,
        input: t.Tensor,
        weights: t.Tensor,
        target: Union[t.Tensor, Any] = None,
        shifted_input: t.Tensor = None,
        ) -> t.Tensor:

            n_chunks = 10

            histograms = []
            for (
                input_chunk,
                target_chunk,
                weights_chunk,
                shifted_input_chunk,
            ) in zip_longest(
                t.chunk(input, chunks=n_chunks, dim=0),
                t.chunk(target, chunks=n_chunks, dim=0) if target is not None else [None],
                t.chunk(weights, chunks=n_chunks, dim=0),
                t.chunk(shifted_input, chunks=n_chunks, dim=0),
                fillvalue=None,
            ):
                histograms.append(
                    self.histograms(
                        input=input_chunk,
                        weights=weights_chunk,
                        shifted_input=shifted_input_chunk,
                    )
                )

            nominal, shifted = tuple(map(lambda it: t.stack(it).sum(axis=0), zip(*histograms)))
            nominal = nominal.detach()
            # return ((nominal - shifted) ** 2).sum() / nominal.sum()
            return (((nominal - shifted) / (t.sqrt(nominal) + 1)) ** 2).sum()


class NeuralNetwork(t.nn.Module):
    def __init__(self, input_dim: int, hidden_nodes: tuple = (100,100,)) -> None:
        super(NeuralNetwork, self).__init__()

        self.input_dim = input_dim
        self.hidden_nodes = hidden_nodes

        self.layer1 = t.nn.Linear(input_dim, hidden_nodes[0])
        self.activation1 = t.nn.ReLU()        
        self.layer2 = t.nn.Linear(hidden_nodes[0], hidden_nodes[1])
        self.activation2 = t.nn.Tanh()
        self.layer3 = t.nn.Linear(self.layer1.out_features, 2)
        self.activation3 = t.nn.Tanh()

    def forward(self, x: t.Tensor) -> t.Tensor:
        x = self.activation1(self.layer1(x))
        x = self.activation2(self.layer2(x))
        x = self.layer3(x)
        x = .001*self.activation3(x) + 1
        return x


def get_mz(corr, mc_t):
    pt_1 = t.mul(corr[:,0], mc_t[:,0]) 
    pt_2 = t.mul(corr[:,1], mc_t[:,1])
    
    px_1 = t.mul(pt_1, t.cos(mc_t[:,4]))
    px_2 = t.mul(pt_2, t.cos(mc_t[:,5]))

    py_1 = t.mul(pt_1, t.sin(mc_t[:,4]))
    py_2 = t.mul(pt_2, t.sin(mc_t[:,5]))

    pz_1 = t.mul(pt_1, t.sinh(mc_t[:,6]))
    pz_2 = t.mul(pt_2, t.sinh(mc_t[:,7]))

    p1 = t.sqrt(t.mul(pt_1, pt_1) + t.mul(pz_1, pz_1))
    p2 = t.sqrt(t.mul(pt_2, pt_2) + t.mul(pz_2, pz_2))
    p1p2 = t.mul(p1, p2)

    d_ang = t.acos(t.div(t.mul(px_1, px_2)+t.mul(py_1, py_2)+t.mul(pz_1, pz_2), p1p2))
    mz = t.sqrt(t.mul(2 * p1p2, 1-t.cos(d_ang)))
    return mz



device = "cuda:0"




seed = np.random.randint(100)
#seed = 61
print(seed)
os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":16:8"
# os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
np.random.seed(seed)
t.manual_seed(seed)

t.cuda.manual_seed(seed)
t.cuda.manual_seed_all(seed)
t.backends.cudnn.deterministic = True
t.backends.cudnn.benchmark = False
t.use_deterministic_algorithms(mode=True)


#### new code

quants = ["pt_1", "pt_2", "q_1", "q_2", "phi_1", "phi_2", "eta_1", "eta_2", "m_vis"]

with uproot.open("/ceph/moh/CROWN_samples/Run3V02/ntuples/2022/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/mm//DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X_0.root:ntuple") as f:
    mc = f.arrays(quants, library='pd')
    mc = mc[(mc['m_vis']>60) & (mc['m_vis']<120)]

with uproot.open('/ceph/moh/CROWN_samples/Run3V02/ntuples/2022/SingleMuon_Run2022C-PromptReco-v1/mm/SingleMuon_Run2022C-PromptReco-v1_0.root:ntuple') as f:
    dt = f.arrays(['m_vis'], library='pd')
    dt = dt[(dt['m_vis']>60) & (dt['m_vis']<120)]
n_tot = min(len(mc), len(dt))
print(len(mc), len(dt), n_tot)

mc_t = t.tensor(mc.to_numpy()[: n_tot][:], dtype=t.float).to(device)
m_vis_dt = t.tensor(dt.to_numpy()[: n_tot][:], dtype=t.float).to(device)
m_vis_dt.requires_grad=True

weights = t.ones_like(m_vis_dt)
bins = np.linspace(60,120,100)


model = NeuralNetwork(len(quants)).to(device)
untrained_model = deepcopy(model)
loss_function = Chi2HistogrammDiffLoss(device=device, bins=bins)
optimizer = t.optim.NAdam(model.parameters(), lr=0.001)
scheduler = t.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', patience=10, factor=.15)
_loss = 10e10
_model = None

for epoch in range(10000):
    try:
        corr = model(mc_t)
        #print(corr[:,0], corr[:,1])
        mz = get_mz(corr, mc_t)
        

        loss = loss_function(
            input=m_vis_dt.squeeze(), #prediction.squeeze(),
            target=None,
            weights=weights.squeeze(),
            shifted_input=mz.squeeze(),
        )
        scheduler.step(loss)


        if loss < _loss and loss != 0.0:
            _loss = loss
            _model = deepcopy(model)
            print(f"{epoch: > 4}: {loss.item(): > 10}")
        
        print(f"{epoch: > 4}: {loss.item(): > 10}", f"corr 1: {t.mean(corr[:,0])} +/- {t.std(corr[:,0])}; factor 2: {t.mean(corr[:,1])} +/- {t.std(corr[:,1])}", end="\r")
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    except KeyboardInterrupt:
        break



import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=(20, 5))
ax[0].hist(m_vis_dt.detach().cpu().numpy(), histtype="step", bins=bins, label="Shifted")
ax[0].hist(get_mz(_model(mc_t), mc_t).detach().cpu().numpy(), bins=bins, label="Nominal")
ax[0].legend()
ax[0].set_title("Trained")
ax[0].set(xlabel="f", ylabel="Count", xlim=(60, 120), ylim=(10, None))#, yscale='log')

ax[1].hist(m_vis_dt.detach().cpu().numpy(), histtype="step", bins=bins, label="Shifted")
ax[1].hist(get_mz(t.ones_like(_model(mc_t)), mc_t).detach().cpu().numpy(), bins=bins, label="Nominal")
ax[1].legend()
ax[1].set_title("Pre-trained")
ax[1].set(xlabel="f", ylabel="Count", xlim=(60, 120), ylim=(10, None))#, yscale='log')

#plt.show()
plt.savefig('train_results.png')