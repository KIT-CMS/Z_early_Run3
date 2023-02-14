import numpy as np
from utils import plot2d

def plot_sf():
    # pt on x and eta on y axis
    sf = np.transpose(np.loadtxt('../correction_files/Run3/mm/res_sf_extra.txt'))
    #print(sf)
    bins = {'eta': [-2.4, -1.2, 0., 1.2, 2.4], 'pt': [25, 30, 35, 39, 42, 45, 49, 60, 120, 9999]}
    
    plot2d(
            matrix = sf, 
            outfile = '../plots/sf.pdf', 
            title = 'Resolution scale factors',
            x = r"$p_\mathrm{T}$ (GeV)",
            xbins = bins['pt'],
            y = r"$\eta$",
            ybins = bins['eta'],
            cmin = 0,
            cmax = 2,
    )

def plot_res(corr):
    # pt on x and eta on y axis
    if corr:
        res_mc = np.transpose(np.loadtxt('../correction_files/Run3/mm/corrected_values/resmz_mc_eta_pt_02_corr.txt'))
        res_dt = np.transpose(np.loadtxt('../correction_files/Run3/mm/corrected_values/resmz_dt_eta_pt_02_corr.txt'))
        corr = '_corr'
        title = 'Resolution ratio mc/dt after correction'
    else:
        res_mc = np.transpose(np.loadtxt('../correction_files/Run3/mm/resmz_mc.txt'))
        res_dt = np.transpose(np.loadtxt('../correction_files/Run3/mm/resmz_dt.txt'))
        corr = ''
        title = 'Resolution ratio mc/dt before correction'
    #print(sf)
    bins = {'eta': [-2.4, -1.2, 0., 1.2, 2.4], 'pt': [25, 30, 35, 39, 42, 45, 49, 60, 120, 9999]}

    plot2d(
            matrix = res_mc/res_dt, 
            outfile = f'../plots/ratio_res{corr}.pdf', 
            title = title,
            x = r"$p_\mathrm{T}$ (GeV)",
            xbins = bins['pt'],
            y = r"$\eta$",
            ybins = bins['eta'],
            cmin = 0.8,
            cmax = 1.2,
    )

def plot_scale(corr):
    # pt on x and eta on y axis
    if corr:
        mc = np.transpose(np.loadtxt('../correction_files/Run3/mm/corrected_values/mz_mc_eta_pt_02_corr.txt'))
        dt = np.transpose(np.loadtxt('../correction_files/Run3/mm/corrected_values/mz_dt_eta_pt_02_corr.txt'))
        corr = '_corr'
        title = 'Peak position ratio mc/dt after correction'
    else:
        mc = np.transpose(np.loadtxt('../correction_files/Run3/mm/mz_mc.txt'))
        dt = np.transpose(np.loadtxt('../correction_files/Run3/mm/mz_dt.txt'))
        corr = ''
        title = 'Peak position ratio mc/dt before correction'
    #print(sf)
    bins = {'eta': [-2.4, -1.2, 0., 1.2, 2.4], 'pt': [25, 30, 35, 39, 42, 45, 49, 60, 120, 9999]}

    plot2d(
            matrix = (91.1876+mc)/(91.1876+dt), 
            outfile = f'../plots/ratio_scale{corr}.pdf', 
            title = title,
            x = r"$p_\mathrm{T}$ (GeV)",
            xbins = bins['pt'],
            y = r"$\eta$",
            ybins = bins['eta'],
            cmin = 0.998,
            cmax = 1.002,
    )


if __name__=='__main__':
    plot_scale(True)
