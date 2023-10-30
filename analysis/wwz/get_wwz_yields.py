import argparse
import pickle
import json
import gzip
import os
import numpy as np
import math
import matplotlib.pyplot as plt
import copy

import hist

from topcoffea.scripts.make_html import make_html

from topcoffea.modules import utils
import topcoffea.modules.MakeLatexTable as mlt

# This script opens a pkl file of histograms produced by wwz processor
# Reads the histograms and dumps out the yields for each group of processes
# Example usage: python get_yld_check.py -f histos/tmp_histo.pkl.gz

# Colors in VVV observation
#ZZ    = (240, 155, 205)  #F09B9B
#ttZ   = (0, 208, 145) #00D091
#WZ    = (163, 155, 47) #A39B2F
#tWZ   = (205, 240, 155) #CDF09B
#Other = (205, 205, 205) #CDCDCD
CLR_LST = ["red","blue","#F09B9B","#00D091","#CDF09B","#CDCDCD"]

#f303abb
EWK_REF = {'WWZ': {'sr_4l_sf_A': (2.248077894200833, 0.008356887900028736), 'sr_4l_sf_B': (1.9092815210624394, 0.0076604307766386034), 'sr_4l_sf_C': (0.5507712043909123, 0.004158433973102092), 'sr_4l_of_1': (0.6239818002213724, 0.004367491375600947), 'sr_4l_of_2': (0.7141527369931282, 0.004685804848413748), 'sr_4l_of_3': (1.433722815978399, 0.006633467314852928), 'sr_4l_of_4': (4.9765614088337315, 0.012434105056347348), 'all_events': (143.03255106410325, 0.06621972573333901), '4l_presel': (25.390131740285142, 0.027935006960828946), 'cr_4l_of': (0.45467763889973867, 0.003688619041601364), 'cr_4l_sf': (1.0234111444715381, 0.0056041066032865225)}, 'ZH': {'sr_4l_sf_A': (0.8808398754063091, 0.008183315119750697), 'sr_4l_sf_B': (1.4252621539799293, 0.009591622240915487), 'sr_4l_sf_C': (0.6304575557069256, 0.005248900545492597), 'sr_4l_of_1': (2.8618885583600786, 0.011746159453499655), 'sr_4l_of_2': (1.2643676003663131, 0.007830100610778018), 'sr_4l_of_3': (0.32845933757653256, 0.004211379560470585), 'sr_4l_of_4': (0.14063188295949658, 0.003254882286889511), 'all_events': (137.65612493509434, 0.09089195409138591), '4l_presel': (19.875361077707566, 0.03255127455912648), 'cr_4l_of': (0.3745634370106927, 0.004910825266362168), 'cr_4l_sf': (0.06583264991149917, 0.001721442231232119)}, 'ZZ': {'sr_4l_sf_A': (1.19626292139219, 0.02635255671225019), 'sr_4l_sf_B': (4.06302963554117, 0.0484072064061236), 'sr_4l_sf_C': (2.4055145616730442, 0.036728105407724476), 'sr_4l_of_1': (0.5727219310065266, 0.01833303036813739), 'sr_4l_of_2': (0.562517372865841, 0.01824400150041214), 'sr_4l_of_3': (0.3522159432868648, 0.01446290428662107), 'sr_4l_of_4': (0.4632858518671128, 0.01671940976648722), 'all_events': (20839.935963596385, 3.5284508076238787), '4l_presel': (2996.3448773944438, 1.30265257357745), 'cr_4l_of': (0.8202006802384858, 0.022285230109117902), 'cr_4l_sf': (1683.1491276815104, 0.9747933823102441)}, 'ttZ': {'sr_4l_sf_A': (1.0880661539267749, 0.061487838616742355), 'sr_4l_sf_B': (0.9408668631222099, 0.055421815101919916), 'sr_4l_sf_C': (0.23643477854784578, 0.028831465995237344), 'sr_4l_of_1': (0.26424446457531303, 0.027872369167180076), 'sr_4l_of_2': (0.2673212867230177, 0.0287057168946335), 'sr_4l_of_3': (0.7172879233257845, 0.04581019242996458), 'sr_4l_of_4': (2.0460971421562135, 0.0831325583461516), 'all_events': (6400.391729545197, 5.524013760526012), '4l_presel': (155.00430985842831, 0.735436023352798), 'cr_4l_of': (58.69226784154307, 0.4534277451200942), 'cr_4l_sf': (0.4853159709600732, 0.04021033112162022)}, 'tWZ': {'sr_4l_sf_A': (0.3471746818977408, 0.015240760391803007), 'sr_4l_sf_B': (0.311000743182376, 0.014440510942444753), 'sr_4l_sf_C': (0.08661178010515869, 0.007642970915514243), 'sr_4l_of_1': (0.1106245769187808, 0.008607028292374249), 'sr_4l_of_2': (0.12694021221250296, 0.009157548303152886), 'sr_4l_of_3': (0.21680891863070428, 0.012050015469807407), 'sr_4l_of_4': (0.7259955864865333, 0.022056413718964097), 'all_events': (457.53748542949324, 0.5535594488742449), '4l_presel': (21.40333431190811, 0.11971397281075291), 'cr_4l_of': (7.203838749264833, 0.06944165652499233), 'cr_4l_sf': (0.1605230116401799, 0.010342360895006718)}, 'other': {'sr_4l_sf_A': (0.6415139838354662, 0.20189232553220066), 'sr_4l_sf_B': (1.228122082655318, 0.37931943188323064), 'sr_4l_sf_C': (0.3666284556966275, 0.186784518924625), 'sr_4l_of_1': (0.6906138394260779, 0.25272547622068714), 'sr_4l_of_2': (0.1325544456485659, 0.12057559170557362), 'sr_4l_of_3': (0.08836240391246974, 0.1843592588174136), 'sr_4l_of_4': (1.3523445471655577, 0.2714217795041188), 'all_events': (2676923.753788651, 2177.829454889236), '4l_presel': (52.9465736518614, 2.5140479882052373), 'cr_4l_of': (2.3255568619351834, 0.2304910831012081), 'cr_4l_sf': (11.619387141661718, 0.7341400279219354)}, '$S/\\sqrt{B}$': {'sr_4l_sf_A': [1.7294976252152816, None], 'sr_4l_sf_B': [1.3036088753027868, None], 'sr_4l_sf_C': [0.6714139887147356, None], 'sr_4l_of_1': [2.7234984522148564, None], 'sr_4l_of_2': [1.8956585995614936, None], 'sr_4l_of_3': [1.5029715209596433, None], 'sr_4l_of_4': [2.3890939134622236, None]}, '$S/\\sqrt{S+B}$': {'sr_4l_sf_A': [1.236626367090044, None], 'sr_4l_sf_B': [1.060990525966462, None], 'sr_4l_sf_C': [0.5712075260687947, None], 'sr_4l_of_1': [1.5399388843353916, None], 'sr_4l_of_2': [1.1295961393615663, None], 'sr_4l_of_3': [0.9949549440590747, None], 'sr_4l_of_4': [1.6426155209771396, None]}, 'Sig': {'sr_4l_sf_A': [3.128917769607142, None], 'sr_4l_sf_B': [3.3345436750423687, None], 'sr_4l_sf_C': [1.181228760097838, None], 'sr_4l_of_1': [3.485870358581451, None], 'sr_4l_of_2': [1.9785203373594413, None], 'sr_4l_of_3': [1.7621821535549316, None], 'sr_4l_of_4': [5.117193291793228, None]}, 'Bkg': {'sr_4l_sf_A': [3.273017741052172, None], 'sr_4l_sf_B': [6.543019324501074, None], 'sr_4l_sf_C': [3.095189576022676, None], 'sr_4l_of_1': [1.6382048119266983, None], 'sr_4l_of_2': [1.0893333174499276, None], 'sr_4l_of_3': [1.3746751891558233, None], 'sr_4l_of_4': [4.587723127675417, None]}, 'Zmetric': {'sr_4l_sf_A': [1.5271304041317981, None], 'sr_4l_sf_B': [1.211362434502514, None], 'sr_4l_sf_C': [0.6343417608521935, None], 'sr_4l_of_1': [2.171342366726005, None], 'sr_4l_of_2': [1.5478844248069015, None], 'sr_4l_of_3': [1.285097369234633, None], 'sr_4l_of_4': [2.0756701298826195, None]}}

SOVERROOTB = "$S/\sqrt{B}$"
SOVERROOTSPLUSB = "$S/\sqrt{S+B}$"

# Keegan's yields as of Aug 18, 2023
KEEGAN_YIELDS = {
    "WWZ" : {
        "sr_4l_sf_A": 2.33,
        "sr_4l_sf_B": 1.97,
        "sr_4l_sf_C": 0.572,
        "sr_4l_of_1": 0.645,
        "sr_4l_of_2": 0.738,
        "sr_4l_of_3": 1.48,
        "sr_4l_of_4": 5.14,
    },
    "ZH" : {
        "sr_4l_sf_A": 0.952,
        "sr_4l_sf_B": 1.52,
        "sr_4l_sf_C": 0.677,
        "sr_4l_of_1": 3.05,
        "sr_4l_of_2": 1.35,
        "sr_4l_of_3": 0.348,
        "sr_4l_of_4": 0.152,
    },
    "ttZ" : {
        "sr_4l_sf_A": 1.32*(0.281/0.2529),
        "sr_4l_sf_B": 1.12*(0.281/0.2529),
        "sr_4l_sf_C": 0.269*(0.281/0.2529),
        "sr_4l_of_1": 0.322*(0.281/0.2529),
        "sr_4l_of_2": 0.382*(0.281/0.2529),
        "sr_4l_of_3": 0.833*(0.281/0.2529),
        "sr_4l_of_4": 2.51*(0.281/0.2529),
    },
    "ZZ": {
        "sr_4l_sf_A": 1.32,
        "sr_4l_sf_B": 4.58,
        "sr_4l_sf_C": 2.78,
        "sr_4l_of_1": 0.623,
        "sr_4l_of_2": 0.619,
        "sr_4l_of_3": 0.39,
        "sr_4l_of_4": 0.503,
    },
    "other": {
        "sr_4l_sf_A": 0.607 + 0.00855,
        "sr_4l_sf_B": 0.417 + 0.0571,
        "sr_4l_sf_C": 0.0967 + 0.0106,
        "sr_4l_of_1": 0.27 + 0.107 + 0.0069,
        "sr_4l_of_2": -0.055 + 0.0 + 0.0106,
        "sr_4l_of_3": -0.0555 + 0.176 + 0.0446,
        "sr_4l_of_4": 0.452 + 0.12 + 0.0233,
    },
}


sample_dict_base = {
    "WWZ" : ["WWZJetsTo4L2Nu"],
    "ZH"  : ["GluGluZH","qqToZHToZTo2L"],

    #"qqZZ": ["ZZTo4l"],
    #"ggZZ": ["ggToZZTo2e2mu", "ggToZZTo2e2tau", "ggToZZTo2mu2tau", "ggToZZTo4e", "ggToZZTo4mu", "ggToZZTo4tau"],
    "ZZ"  : ["ZZTo4l", "ggToZZTo2e2mu", "ggToZZTo2e2tau", "ggToZZTo2mu2tau", "ggToZZTo4e", "ggToZZTo4mu", "ggToZZTo4tau"],

    "ttZ" : [
        "TTZToLL_M_1to10",
        "TTZToLLNuNu_M_10",
        "TTZToQQ",
    ],

    "tWZ" : ["tWll"],

    "other" : [

        ##"WWZJetsTo4L2Nu",
        ##"GluGluZH","qqToZHToZTo2L",
        ##"ZZTo4l", "ggToZZTo2e2mu", "ggToZZTo2e2tau", "ggToZZTo2mu2tau", "ggToZZTo4e", "ggToZZTo4mu", "ggToZZTo4tau",
        ##"TTZToLL_M_1to10","TTZToLLNuNu_M_10","TTZToQQ",
        ##"tWll",

        ##"DYJetsToLL_M_10to50_MLM",
        "DYJetsToLL_M_50_MLM",
        "SSWW",
        "ST_antitop_t-channel",
        "ST_top_s-channel",
        "ST_top_t-channel",
        "tbarW_noFullHad",
        "ttHnobb",
        "TTTo2L2Nu",
        "TTWJetsToLNu",
        "TTWJetsToQQ",
        "tW_noFullHad",
        "tZq",
        "VHnobb",
        ##"WJetsToLNu",
        "WWTo2L2Nu",
        "WZTo3LNu",

        "WWW",
        "WZZ",
        "ZZZ",
    ],
}

# Processes indiviudally
sample_dict_base_indiv = {
    "WWZJetsTo4L2Nu":            ["WWZJetsTo4L2Nu"],
    "GluGluZH":                  ["GluGluZH"],
    "qqToZHToZTo2L":             ["qqToZHToZTo2L"],
    "ZZTo4l":                    ["ZZTo4l"],
    "ggToZZTo2e2mu":             ["ggToZZTo2e2mu"],
    "ggToZZTo2e2tau":            ["ggToZZTo2e2tau"],
    "ggToZZTo2mu2tau":           ["ggToZZTo2mu2tau"],
    "ggToZZTo4e":                ["ggToZZTo4e"],
    "ggToZZTo4mu":               ["ggToZZTo4mu"],
    "ggToZZTo4tau":              ["ggToZZTo4tau"],
    "TTZToLL_M_1to10":           ["TTZToLL_M_1to10"],
    "TTZToLLNuNu_M_10":          ["TTZToLLNuNu_M_10"],
    "TTZToQQ":                   ["TTZToQQ"],
    "tWll" :                     ["tWll"],

    ##"DYJetsToLL_M_10to50_MLM": ["DYJetsToLL_M_10to50_MLM"],
    "DYJetsToLL_M_50_MLM":       ["DYJetsToLL_M_50_MLM"],
    "SSWW":                      ["SSWW"],
    "ST_antitop_t-channel":      ["ST_antitop_t-channel"],
    "ST_top_s-channel":          ["ST_top_s-channel"],
    "ST_top_t-channel":          ["ST_top_t-channel"],
    "tbarW_noFullHad":           ["tbarW_noFullHad"],
    "ttHnobb":                   ["ttHnobb"],
    "TTTo2L2Nu":                 ["TTTo2L2Nu"],
    "TTWJetsToLNu":              ["TTWJetsToLNu"],
    "TTWJetsToQQ":               ["TTWJetsToQQ"],
    "tW_noFullHad":              ["tW_noFullHad"],
    "tZq":                       ["tZq"],
    "VHnobb":                    ["VHnobb"],
    ##"WJetsToLNu":              ["WJetsToLNu"],
    "WWTo2L2Nu":                 ["WWTo2L2Nu"],
    "WZTo3LNu":                  ["WZTo3LNu"],

    "WWW" : ["WWW"],
    "WZZ" : ["WZZ"],
    "ZZZ" : ["ZZZ"],
}



# Pass dictionary with the base names for the samples, and return with full list for 4 years
def create_full_sample_dict(in_dict):
    out_dict = {}
    years = ["UL16APV","UL16","UL17","UL18"]
    for proc_group in in_dict.keys():
        out_dict[proc_group] = []
        for proc_base_name in in_dict[proc_group]:
            for year_str in years:
                out_dict[proc_group].append(f"{year_str}_{proc_base_name}")
    return out_dict


sample_dict = create_full_sample_dict(sample_dict_base)




################### Getting and printing yields ###################

# Get the yields in the SR
def get_yields(histos_dict,raw_counts=False,quiet=False):

    yld_dict = {}

    # Look at the yields in one histo (e.g. njets)
    if raw_counts: dense_axis = "njets_counts"
    else: dense_axis = "njets"
    for proc_name in sample_dict.keys():
        yld_dict[proc_name] = {}
        for cat_name in histos_dict[dense_axis].axes["category"]:
            val = sum(sum(histos_dict[dense_axis][{"category":cat_name,"process":sample_dict[proc_name]}].values(flow=True)))
            var = np.sqrt(sum(sum(histos_dict[dense_axis][{"category":cat_name,"process":sample_dict[proc_name]}].variances(flow=True))))
            yld_dict[proc_name][cat_name] = (val,var)

    # Print to screen
    if not quiet:
        for proc in yld_dict.keys():
            print(f"\n{proc}:")
            for cat in yld_dict[proc].keys():
                val = yld_dict[proc][cat]
                print(f"\t{cat}: {val}")

    return yld_dict


# Gets the S/sqrt(B) and puts it into the dict
# Hard coded for the summed values (e.g. looking for "ZH" not "GluGluZH","qqToZHToZTo2L")
def put_s_over_root_b(yld_dict):
    sig_lst = ["WWZ","ZH"]
    bkg_lst = ["ZZ","ttZ","tWZ","other"]
    sig_sum = {"sr_4l_sf_A":0, "sr_4l_sf_B":0, "sr_4l_sf_C":0, "sr_4l_of_1":0, "sr_4l_of_2":0, "sr_4l_of_3":0, "sr_4l_of_4":0}
    bkg_sum = {"sr_4l_sf_A":0, "sr_4l_sf_B":0, "sr_4l_sf_C":0, "sr_4l_of_1":0, "sr_4l_of_2":0, "sr_4l_of_3":0, "sr_4l_of_4":0}
    for proc in yld_dict.keys():
        print(proc)
        for cat in yld_dict[proc].keys():
            if cat not in sig_sum: continue
            val,err = yld_dict[proc][cat]
            print("   ",cat,val)
            if proc in sig_lst:
                sig_sum[cat] += val
            if proc in bkg_lst:
                bkg_sum[cat] += val

    yld_dict[SOVERROOTB] = {}
    yld_dict[SOVERROOTSPLUSB] = {}
    yld_dict["Sig"] = {}
    yld_dict["Bkg"] = {}
    yld_dict["Zmetric"] = {}
    for cat in sig_sum.keys():
        s = sig_sum[cat]
        b = bkg_sum[cat]
        yld_dict[SOVERROOTB][cat]      = [s/math.sqrt(b) , None]
        yld_dict[SOVERROOTSPLUSB][cat] = [s/math.sqrt(s+b) , None]
        yld_dict["Sig"][cat] = [s, None]
        yld_dict["Bkg"][cat] = [b, None]
        yld_dict["Zmetric"][cat] = [math.sqrt(2 * ((s + b) * math.log(1 + s / b) - s)), None] # Eq 18 https://cds.cern.ch/record/2203244/files/1087459_109-114.pdf


# Print yields
def print_yields(yld_dict,print_fom=True):

    # Dump the yields to dict for latex table
    yld_dict_for_printing = {}
    cats_to_print = [ "sr_4l_sf_A", "sr_4l_sf_B", "sr_4l_sf_C", "sr_4l_of_1", "sr_4l_of_2", "sr_4l_of_3", "sr_4l_of_4"]
    for proc in yld_dict.keys():
        yld_dict_for_printing[proc] = {}
        for cat in yld_dict[proc].keys():
            if cat not in cats_to_print: continue
            yld_dict_for_printing[proc][cat] = yld_dict[proc][cat]


    # Print the yields directly
    mlt.print_latex_yield_table(
        yld_dict_for_printing,
        tag="All yields",
        key_order=sample_dict_base.keys(),
        subkey_order=cats_to_print,
        print_begin_info=True,
        print_end_info=True,
        print_errs=True,
        column_variable="subkeys",
    )
    #exit()


    # Compare with other yields, print comparison

    #tag1 = "ewkcoffea"
    #tag2 = "VVVNanoLooper"
    tag1 = "lepRecoSFlepTightSF"
    tag2 = "noSF " # Hard coded ref

    #yld_dict_comp = utils.put_none_errs(KEEGAN_YIELDS)
    yld_dict_comp = EWK_REF

    yld_dict_1 = yld_dict_for_printing
    yld_dict_2 = yld_dict_comp

    pdiff_dict = utils.get_diff_between_nested_dicts(yld_dict_1,yld_dict_2,difftype="percent_diff",inpercent=True)
    diff_dict  = utils.get_diff_between_nested_dicts(yld_dict_1,yld_dict_2,difftype="absolute_diff")

    #utils.print_yld_dicts(yld_dict_1,tag1)
    #utils.print_yld_dicts(yld_dict_2,tag2)
    #utils.print_yld_dicts(pdiff_dict,f"Percent diff between {tag1} and {tag2}")
    #utils.print_yld_dicts(diff_dict,f"Diff between {tag1} and {tag2}")

    procs_to_print = list(sample_dict_base.keys())
    if print_fom:
        procs_to_print.append("Sig")
        procs_to_print.append("Bkg")
        procs_to_print.append(SOVERROOTB)
        procs_to_print.append(SOVERROOTSPLUSB)
        procs_to_print.append("Zmetric")

    mlt.print_begin()
    mlt.print_latex_yield_table(yld_dict_1,key_order=procs_to_print,subkey_order=cats_to_print,tag=tag1)
    mlt.print_latex_yield_table(yld_dict_2,key_order=procs_to_print,subkey_order=cats_to_print,tag=tag2)
    mlt.print_latex_yield_table(pdiff_dict,key_order=procs_to_print,subkey_order=cats_to_print,tag=f"Percent diff between {tag1} and {tag2}")
    mlt.print_latex_yield_table(diff_dict, key_order=procs_to_print,subkey_order=cats_to_print,tag=f"Diff between {tag1} and {tag2}")
    mlt.print_end()


# Dump the counts dict to a latex table
def print_counts(counts_dict):

    cats_to_print = ["all_events", "4l_presel", "sr_4l_sf_A", "sr_4l_sf_B", "sr_4l_sf_C", "sr_4l_of_1", "sr_4l_of_2", "sr_4l_of_3", "sr_4l_of_4"]

    # Print the yields directly
    mlt.print_latex_yield_table(
        counts_dict,
        tag="Raw MC counts (ewkcoffea)",
        key_order=counts_dict.keys(),
        subkey_order=cats_to_print,
        print_begin_info=True,
        print_end_info=True,
        column_variable="subkeys",
    )


# This should maybe be in a different script
################### Hist manipulation and plotting ###################


# Get the list of categories on the sparese axis
def get_axis_cats(histo,axis_name):
    process_list = [x for x in histo.axes[axis_name]]
    return process_list


# Merges the last bin (overflow) into the second to last bin, zeros the content of the last bin, returns a new hist
def merge_overflow(hin):
    hout = copy.deepcopy(hin)
    for cat_idx,arr in enumerate(hout.values(flow=True)):
        hout.values(flow=True)[cat_idx][-2] += hout.values(flow=True)[cat_idx][-1]
        hout.values(flow=True)[cat_idx][-1] = 0
        hout.variances(flow=True)[cat_idx][-2] += hout.variances(flow=True)[cat_idx][-1]
        hout.variances(flow=True)[cat_idx][-1] = 0
    return hout


# Rebin according to https://github.com/CoffeaTeam/coffea/discussions/705
def rebin(histo,factor):
    return histo[..., ::hist.rebin(factor)]


# Regroup categories (e.g. processes)
def group(h, oldname, newname, grouping):

    # Build up a grouping dict that drops any proc that is not in our h
    grouping_slim = {}
    #proc_lst = get_axis_cats(h,"process")
    proc_lst = get_axis_cats(h,oldname)
    for grouping_name in grouping.keys():
        for proc in grouping[grouping_name]:
            if proc in proc_lst:
                if grouping_name not in grouping_slim:
                    grouping_slim[grouping_name] = []
                grouping_slim[grouping_name].append(proc)
            #else:
            #    print(f"WARNING: process {proc} not in this hist")

    # From Nick: https://github.com/CoffeaTeam/coffea/discussions/705#discussioncomment-4604211
    hnew = hist.Hist(
        hist.axis.StrCategory(grouping_slim, name=newname),
        *(ax for ax in h.axes if ax.name != oldname),
        storage=h.storage_type(),
    )
    for i, indices in enumerate(grouping_slim.values()):
        hnew.view(flow=True)[i] = h[{oldname: indices}][{oldname: sum}].view(flow=True)

    return hnew


# Takes a mc hist and data hist and plots both
def make_cr_fig(histo_mc,histo_data=None,title="test",unit_norm_bool=False):

    # Create the figure
    fig, (ax, rax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(7,7),
        gridspec_kw={"height_ratios": (3, 1)},
        sharex=True
    )
    fig.subplots_adjust(hspace=.07)

    # Plot the mc
    histo_mc.plot1d(
        stack=True,
        histtype="fill",
        color=CLR_LST,
        ax=ax,
    )
    # Plot the data
    if histo_data is not None:
        histo_data.plot1d(
            stack=False,
            histtype="errorbar",
            color="k",
            ax=ax,
            w2=histo_data.variances(),
            w2method="sqrt",
        )
    # Plot a dummy hist on rax to get the label to show up
    histo_mc.plot1d(alpha=0, ax=rax)

    ### Get the errs on MC and plot them by hand ###
    histo_mc_sum = histo_mc[{"process_grp":sum}]
    mc_arr = histo_mc_sum.values()
    mc_err_arr = np.sqrt(histo_mc_sum.variances())
    err_p = np.append(mc_arr + mc_err_arr, 0)
    err_m = np.append(mc_arr - mc_err_arr, 0)
    bin_edges_arr = histo_mc_sum.axes[0].edges
    bin_centers_arr = histo_mc_sum.axes[0].centers
    ax.fill_between(bin_edges_arr,err_m,err_p, step='post', facecolor='none', edgecolor='gray', alpha=0.5, linewidth=0.0, label='MC stat', hatch='/////')

    ### Get the errs on data and ratios and plot them by hand ###
    if histo_data is not None:
        histo_data_sum = histo_data[{"process_grp":sum}]

        data_arr = histo_data_sum.values()
        data_err_arr = np.sqrt(histo_data_sum.variances())

        err_ratio_p = np.append(1+mc_err_arr/mc_arr,1)
        err_ratio_m = np.append(1-mc_err_arr/mc_arr,1)

        data_ratio_err_p = (data_arr + data_err_arr)/mc_arr
        data_ratio_err_m = (data_arr - data_err_arr)/mc_arr

        rax.fill_between(bin_edges_arr,err_ratio_m,err_ratio_p,step='post', facecolor='none',edgecolor='gray', label='MC stat', linewidth=0.0, hatch='/////',alpha=0.5)
        rax.scatter(bin_centers_arr,data_arr/mc_arr,facecolor='black',edgecolor='black',marker="o")
        rax.vlines(bin_centers_arr,data_ratio_err_p,data_ratio_err_m,color='k')

    # Scale the y axis and labels
    ax.legend(fontsize="12")
    ax.set_title(title)
    ax.autoscale(axis='y')
    ax.set_xlabel(None)
    rax.set_ylabel('Ratio')
    rax.set_ylim(0.0,2.0)
    rax.axhline(1.0,linestyle="-",color="k",linewidth=1)
    ax.tick_params(axis='y', labelsize=16)
    rax.tick_params(axis='x', labelsize=16)
    #ax.set_yscale('log')

    return fig


# Plots a hist
def make_single_fig(histo_mc,title,unit_norm_bool=False):
    #print("\nPlotting values:",histo.values())
    fig, ax = plt.subplots(1, 1, figsize=(7,7))

    # Guard against accidently looking at data before unblind
    proc_lst = get_axis_cats(histo_mc,"process_grp")
    for proc_name in proc_lst:
        if "data" in proc_name:
            raise Exception(f"CAUTION: Are you sure you want to look at this data in {proc_name}?")

    # Plot the mc
    histo_mc.plot1d(
        stack=True,
        histtype="fill",
        #color=["red","yellow","blue","green","orange","grey"],
        color=CLR_LST,
        yerr=True,
    )

    plt.title(title)
    ax.autoscale(axis='y')
    plt.legend()
    return fig


# Main function for making CR plots
def make_plots(histo_dict,save_dir_path):

    for var_name in histo_dict.keys():
        print(f"\n{var_name}")
        if "counts" in var_name: continue
        #if var_name != "njets": continue # TMP
        histo = histo_dict[var_name]

        # Rebin if continous variable
        if var_name not in ["njets","nbtagsl","nleps"]:
            histo = rebin(histo,6)

        # Group SR procs together
        #grouping_sr_procs = {"sr_4l_sf":["sr_4l_sf_A","sr_4l_sf_B","sr_4l_sf_C"],"sr_4l_of":["sr_4l_of_1","sr_4l_of_2","sr_4l_of_3","sr_4l_of_4"]}
        #histo = group(histo,"category","category",grouping_sr_procs)

        # Loop over categories and make plots for each
        for cat_name in histo.axes["category"]:
            if cat_name not in ["cr_4l_sf","cr_4l_of"]: continue # TMP
            print(cat_name)

            histo_cat = histo[{"category":cat_name}]

            # Group the mc samples
            grouping_mc = sample_dict
            histo_grouped_mc = group(histo_cat,"process","process_grp",grouping_mc)

            # Group the data samples
            grouping_data = {'data': ["UL16APV_data","UL16_data","UL17_data","UL18_data"]}
            histo_grouped_data = group(histo_cat,"process","process_grp",grouping_data)

            #####
            #if cat_name == "cr_4l_of" and var_name == "nleps":
            #    print("mc\n",histo_grouped_mc)
            #    print("data\n",histo_grouped_data)
            #    print("val mc\n",histo_grouped_mc.values(flow=True))
            #    print("val data\n",histo_grouped_data.values(flow=True))
            #    print("var mc\n",(histo_grouped_mc.variances(flow=True)))
            #    print("var data\n",(histo_grouped_data.variances(flow=True)))
            #continue
            #print("\Yields")
            #print(type(histo_cat.values(flow=True)))
            #print(sum(histo_cat.values(flow=True)))
            #print(sum(sum(histo_cat.values(flow=True))))
            #exit()
            ####

            # Merge overflow into last bin (so it shows up in the plot)
            histo_grouped_data = merge_overflow(histo_grouped_data)
            histo_grouped_mc = merge_overflow(histo_grouped_mc)


            # Make figure
            title = f"{cat_name}_{var_name}"
            if "cr" in title:
                fig = make_cr_fig(histo_grouped_mc,histo_grouped_data,title=title)
            else:
                fig = make_cr_fig(histo_grouped_mc,title=title)

            # Save
            save_dir_path_cat = os.path.join(save_dir_path,cat_name)
            if not os.path.exists(save_dir_path_cat): os.mkdir(save_dir_path_cat)
            fig.savefig(os.path.join(save_dir_path_cat,title+".pdf"))
            fig.savefig(os.path.join(save_dir_path_cat,title+".png"))

            make_html(os.path.join(os.getcwd(),save_dir_path_cat))



################### Main ###################

def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_file_path", help = "The path to the pkl file")
    parser.add_argument("-o", "--output-path", default=".", help = "The path the output files should be saved to")
    parser.add_argument('-y', "--get-yields", action='store_true', help = "Get yields from the pkl file")
    parser.add_argument('-p', "--make-plots", action='store_true', help = "Make plots from the pkl file")
    args = parser.parse_args()

    # Get the counts from the input hiso
    histo_dict = pickle.load(gzip.open(args.pkl_file_path))

    # Wrapper around the code for getting the raw counts and dump to latex table
    #counts_dict = get_yields(histo_dict,raw_counts=True)
    #print_counts(counts_dict)
    #exit()

    out_path = "plots" # Could make this an argument

    # Wrapper around the code for getting the yields for sr and bkg samples
    if args.get_yields:
        yld_dict = get_yields(histo_dict)
        put_s_over_root_b(yld_dict)
        print_yields(yld_dict)

        # Dump yield dict to json
        json_name = "process_yields.json" # Could be an argument
        json_name = os.path.join(out_path,json_name)
        with open(json_name,"w") as out_file: json.dump(yld_dict, out_file, indent=4)
        print(f"\nSaved json file: {json_name}\n")

    # Make plots
    if args.make_plots:
        make_plots(histo_dict,save_dir_path=out_path)




if __name__ == "__main__":
    main()

