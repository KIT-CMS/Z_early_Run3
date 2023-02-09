"""
collecting some util scripts here
"""

import pandas as pd
from tabulate import tabulate

def FormatOutputForWZ(istring: str, tdrStyle: bool=False):
    labelmaps = {}
    if tdrStyle:
        labelmaps['muplus'] = '\\PWp'
        labelmaps['muminus'] = '\\PWm'
        labelmaps['mumu'] = '\\PZ'
        labelmaps['eplus'] = '\\PWp'
        labelmaps['eminus'] = '\\PWm'
        labelmaps['ee'] = '\\PZ'
        labelmaps['lepplus'] = '\\PWp'
        labelmaps['lepminus'] = '\\PWm'
        labelmaps['leplep'] = '\\PZ'
        labelmaps['Wp'] = '\\PWp'
        labelmaps['Wm'] = '\\PWm'
        labelmaps['Winc'] = '\\PW'
        labelmaps['WOverZ'] = '\\PW / \\PZ'
        labelmaps['WpOverZ'] = '\\PWp / \\PZ'
        labelmaps['WmOverZ'] = '\\PWm / \\PZ'
        labelmaps['WpOverWm'] = '\\PWp / \\PWm'
    else:
        labelmaps['muplus'] = '$\\mathrm{W}^{+}\\rightarrow\\mu^{+}\\nu$'
        labelmaps['muminus'] = '$\\mathrm{W}^{-}\\rightarrow\\mu^{-}\\bar{\\nu}$'
        labelmaps['mumu'] = '$\\mathrm{Z}\\rightarrow\\mu^{+}\\mu^{-}$'
        labelmaps['eplus'] = '$\\mathrm{W}^{+}\\rightarrow e^{+}\\nu$'
        labelmaps['eminus'] = '$\\mathrm{W}^{-}\\rightarrow e^{-}\\bar{\\nu}$'
        labelmaps['ee'] = '$\\mathrm{Z}\\rightarrow e^{+}e^{-}$'
        labelmaps['lepplus'] = '$\\mathrm{W}^{+}\\rightarrow \\ell^{+}\\nu$'
        labelmaps['lepminus'] = '$\\mathrm{W}^{-}\\rightarrow \\ell^{-}\\bar{\\nu}$'
        labelmaps['leplep'] = '$\\mathrm{Z}\\rightarrow \\ell^{+}\\ell^{-}$'
        labelmaps['Wp'] = '$\\mathrm{W}^{+}\\rightarrow \\ell^{+}\\nu$'
        labelmaps['Wm'] = '$\\mathrm{W}^{-}\\rightarrow \\ell^{-}\\bar{\\nu}$'
        labelmaps['Winc'] = '$\\mathrm{W}^{\\pm}\\rightarrow \\ell^{\\pm}\\nu$'
        labelmaps['WOverZ'] = '$\\mathrm{W}^{\pm}/\\mathrm{Z}$'
        labelmaps['WpOverWm'] = '$\\mathrm{W}^{+}/\\mathrm{W}^{-}$'


    procmaps = {}
    procmaps['data'] = 'Data'
    procmaps['sig'] = "Signal"
    procmaps['ewk'] = "EWK"
    procmaps['qcd'] = "QCD"
    procmaps['ttbar'] = "$t\\bar{t}$"

    sysmaps = {}
    sysmaps['lumi         '] = 'Luminosity             '
    sysmaps['recoil       '] = 'Recoil                 '
    sysmaps['QCDbkg       '] = 'QCD background         '
    sysmaps['effstat      '] = 'Efficiency Stat.       '
    sysmaps['prefire      '] = 'Prefire                '
    sysmaps['QCDscale     '] = 'QCD scale              '
    sysmaps['pdfscale     '] = 'PDF + QCD scale        '
    sysmaps['Momentum     '] = 'Muon momentum          '
    sysmaps['stat         '] = 'Stat. uncertainty      '
    sysmaps['binByBinStat '] = 'Bin-by-bin uncertainty '
    sysmaps['effsys       '] = 'Muon efficiency        '
    sysmaps['pdfalphaS    '] = 'PDF + $\\alpha_\\mathrm{S}$'
    sysmaps['mcsec        '] = 'MC normalization       '

    if tdrStyle:
        for key in ['lepplus','lepminus','leplep','WpOverZ','WmOverZ','WpOverWm']:
            istring = istring.replace(key, labelmaps[key])
    else:
        for key in labelmaps.keys():
            istring = istring.replace(key, labelmaps[key])

    for key in procmaps.keys():
        istring = istring.replace(key, procmaps[key])

    for key in sysmaps.keys():
        istring = istring.replace(key, sysmaps[key])

    return istring


def FormatTable(pdict: str, columns: list = None, caption: str = None, label: str = None, precision: int=1, simple: bool=False, tdrStyle: bool=False):
    """
    given a dictionary, print the latex version of the table
    """
    df = pd.DataFrame(pdict, columns=columns)
    df = df.round(precision)
    if simple:
        output = tabulate(df, headers='keys', tablefmt='psql')

    if tdrStyle:
        output = df.to_latex(caption = caption, label = label)
        output = output.replace('\\toprule', '\\hline').replace('\\midrule', '\\hline').replace('\\bottomrule','\\hline')
        output = FormatOutputForWZ(output, tdrStyle=True)

    else:
        output = df.to_latex(caption = caption, label = label)
        output = output.replace('\\toprule', '\\hline').replace('\\midrule', '\\hline').replace('\\bottomrule','\\hline')
        output = FormatOutputForWZ(output)

    return output
