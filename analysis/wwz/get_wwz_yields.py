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

BDT_INPUT_LST = [
    "mll_wl0_wl1",
    "dphi_4l_met",
    "dphi_zleps_met",
    "dphi_wleps_met",
    "dr_wl0_wl1",
    "dr_zl0_zl1",
    "dr_wleps_zleps",
    "met",
    "mt2",
    "ptl4",
    "scalarptsum_lepmet",
    "scalarptsum_lepmetjet",
    "z_lep0_pt",
    "z_lep1_pt",
    "w_lep0_pt",
    "w_lep1_pt",
]

TMP_LST = [
    "j0pt",
    "njets",
    "nbtagsl",
    "nleps",
    "met",
    "l0pt",
]

# No SFs (e.g. oct27_new_sr_presel_cats.pkl.gz)
#EWK_REF = {'WWZ': {'all_events': (143.03255106410325, 0.0043850520761986814), '4l_presel': (25.390131740285142, 0.0007803646139015601), 'cr_4l_sf': (1.0234111444715381, 3.140601082099961e-05), 'cr_4l_of': (0.45467763889973867, 1.3605910434064165e-05), 'sr_4l_sf_C': (0.5507712043909123, 1.7292573108649655e-05), 'sr_4l_sf_presel': (4.887572965786603, 0.0001514161312930492), 'sr_4l_of_2': (0.7141527369931282, 2.1956767077417794e-05), 'sr_4l_of_presel': (8.977618631754012, 0.0002774286898025886), 'sr_4l_sf_A': (2.248077894200833, 6.983757537364668e-05), 'sr_4l_of_4': (4.9765614088337315, 0.00015460696855228265), 'sr_4l_sf_B': (1.9092815210624394, 5.8682199683671906e-05), 'sr_4l_of_1': (0.6239818002213724, 1.907498091594866e-05), 'sr_4l_of_3': (1.433722815978399, 4.4002888617222127e-05), 'sr_sf_all': [4.708130619654185, 0.00014581234816596824], 'sr_of_all': [7.748418762026631, 0.00023964160516287122], 'sr_all': [12.456549381680816, 0.0003854539533288395]}, 'ZH': {'all_events': (137.65612493509434, 0.008261347318550601), '4l_presel': (19.875361077707566, 0.0010595854754236352), 'cr_4l_sf': (0.06583264991149917, 2.963363355469417e-06), 'cr_4l_of': (0.3745634370106927, 2.4116204796741056e-05), 'sr_4l_sf_C': (0.6304575557069256, 2.7550956936472485e-05), 'sr_4l_sf_presel': (3.124085131260472, 0.00019582211706122878), 'sr_4l_of_2': (1.2643676003663131, 6.131047557490629e-05), 'sr_4l_of_presel': (7.463838904360273, 0.0003795543171914658), 'sr_4l_sf_A': (0.8808398754063091, 6.696664634914036e-05), 'sr_4l_of_4': (0.14063188295949658, 1.0594258701507094e-05), 'sr_4l_sf_B': (1.4252621539799293, 9.199921721242461e-05), 'sr_4l_of_1': (2.8618885583600786, 0.0001379722619070393), 'sr_4l_of_3': (0.32845933757653256, 1.7735717802349422e-05), 'sr_sf_all': [2.936559585093164, 0.00018651682049803745], 'sr_of_all': [4.595347379262421, 0.0002276127139858021], 'sr_all': [7.531906964355585, 0.0004141295344838396]}, 'ZZ': {'all_events': (20839.935963596385, 12.449965101821629), '4l_presel': (2996.3448773944438, 1.6969037274479533), 'cr_4l_sf': (1683.1491276815104, 0.9502221381958456), 'cr_4l_of': (0.8202006802384858, 0.0004966314810163349), 'sr_4l_sf_C': (2.4055145616730442, 0.0013489537268409199), 'sr_4l_sf_presel': (11.824718278345244, 0.00681902139812232), 'sr_4l_of_2': (0.562517372865841, 0.00033284359074704046), 'sr_4l_of_presel': (19.06702512973061, 0.011167948980671838), 'sr_4l_sf_A': (1.19626292139219, 0.0006944572452723625), 'sr_4l_of_4': (0.4632858518671128, 0.0002795386629397082), 'sr_4l_sf_B': (4.06302963554117, 0.002343257632045054), 'sr_4l_of_1': (0.5727219310065266, 0.00033610000247904766), 'sr_4l_of_3': (0.3522159432868648, 0.00020917560040396212), 'sr_sf_all': [7.664807118606404, 0.004386668604158336], 'sr_of_all': [1.9507410990263452, 0.0011576578565697584], 'sr_all': [9.61554821763275, 0.005544326460728095]}, 'ttZ': {'all_events': (6400.391729545197, 30.51472802648067), '4l_presel': (155.00430985842831, 0.5408661444449773), 'cr_4l_sf': (0.4853159709600732, 0.0016168707289103397), 'cr_4l_of': (58.69226784154307, 0.2055967200446931), 'sr_4l_sf_C': (0.23643477854784578, 0.0008312534314345273), 'sr_4l_sf_presel': (2.415647323941812, 0.008086969492348196), 'sr_4l_of_2': (0.2673212867230177, 0.0008240181824348471), 'sr_4l_of_presel': (3.787581167416647, 0.012204073027489382), 'sr_4l_sf_A': (1.0880661539267749, 0.003780754297758552), 'sr_4l_of_4': (2.0460971421562135, 0.006911022257176301), 'sr_4l_sf_B': (0.9408668631222099, 0.0030715775891913986), 'sr_4l_of_1': (0.26424446457531303, 0.0007768689629915706), 'sr_4l_of_3': (0.7172879233257845, 0.0020985737304703837), 'sr_sf_all': [2.2653677955968305, 0.007683585318384478], 'sr_of_all': [3.2949508167803288, 0.010610483133073102], 'sr_all': [5.560318612377159, 0.018294068451457583]}, 'tWZ': {'all_events': (457.53748542949324, 0.30642806343796253), '4l_presel': (21.40333431190811, 0.014331435286133691), 'cr_4l_sf': (0.1605230116401799, 0.00010696442888256418), 'cr_4l_of': (7.203838749264833, 0.00482214366093501), 'sr_4l_sf_C': (0.08661178010515869, 5.8415004415396626e-05), 'sr_4l_sf_presel': (0.7771313883713447, 0.0005206595096742938), 'sr_4l_of_2': (0.12694021221250296, 8.386069092457832e-05), 'sr_4l_of_presel': (1.3504042699933052, 0.0009034388377762367), 'sr_4l_sf_A': (0.3471746818977408, 0.00023228077732035136), 'sr_4l_of_4': (0.7259955864865333, 0.00048648538614210764), 'sr_4l_sf_B': (0.311000743182376, 0.00020852835627886666), 'sr_4l_of_1': (0.1106245769187808, 7.408093602573077e-05), 'sr_4l_of_3': (0.21680891863070428, 0.00014520287282259782), 'sr_sf_all': [0.7447872051852755, 0.0004992241380146146], 'sr_of_all': [1.1803692942485213, 0.0007896298859150146], 'sr_all': [1.9251564994337969, 0.001288854023929629]}, 'other': {'all_events': (2676923.753788651, 4742941.134583146), '4l_presel': (52.9465736518614, 6.320437286998802), 'cr_4l_sf': (11.619387141661718, 0.5389615805972202), 'cr_4l_of': (2.3255568619351834, 0.05312613938916803), 'sr_4l_sf_C': (0.3666284556966275, 0.0348884565099036), 'sr_4l_sf_presel': (2.418255006778054, 0.22476817034763386), 'sr_4l_of_2': (0.1325544456485659, 0.014538473315149195), 'sr_4l_of_presel': (3.819963646121323, 0.32266726616108504), 'sr_4l_sf_A': (0.6415139838354662, 0.04076051110880008), 'sr_4l_of_4': (1.3523445471655577, 0.07366978238918248), 'sr_4l_sf_B': (1.228122082655318, 0.14388323140421685), 'sr_4l_of_1': (0.6906138394260779, 0.06387016633097309), 'sr_4l_of_3': (0.08836240391246974, 0.0339883363117061), 'sr_sf_all': [2.236264522187412, 0.21953219902292054], 'sr_of_all': [2.2638752361526713, 0.18606675834701086], 'sr_all': [4.500139758340083, 0.4055989573699314]}, '$S/\\sqrt{B}$': {'sr_4l_sf_A': [1.7294976252152816, None], 'sr_4l_sf_B': [1.3036088753027868, None], 'sr_4l_sf_C': [0.6714139887147356, None], 'sr_4l_of_1': [2.7234984522148564, None], 'sr_4l_of_2': [1.8956585995614936, None], 'sr_4l_of_3': [1.5029715209596433, None], 'sr_4l_of_4': [2.3890939134622236, None], 'sr_sf_all': [2.267455595956694, None], 'sr_of_all': [4.356335439961683, None], 'sr_all': [4.911111212862266, None]}, '$S/\\sqrt{S+B}$': {'sr_4l_sf_A': [1.236626367090044, None], 'sr_4l_sf_B': [1.060990525966462, None], 'sr_4l_sf_C': [0.5712075260687947, None], 'sr_4l_of_1': [1.5399388843353916, None], 'sr_4l_of_2': [1.1295961393615663, None], 'sr_4l_of_3': [0.9949549440590747, None], 'sr_4l_of_4': [1.6426155209771396, None], 'sr_sf_all': [1.7266220506557137, None], 'sr_of_all': [2.708416566189788, None], 'sr_all': [3.2119688668824025, None]}, 'Sig': {'sr_4l_sf_A': [3.128917769607142, 0.00013680422172278702], 'sr_4l_sf_B': [3.3345436750423687, 0.00015068141689609653], 'sr_4l_sf_C': [1.181228760097838, 4.484353004512214e-05], 'sr_4l_of_1': [3.485870358581451, 0.00015704724282298796], 'sr_4l_of_2': [1.9785203373594413, 8.326724265232408e-05], 'sr_4l_of_3': [1.7621821535549316, 6.173860641957155e-05], 'sr_4l_of_4': [5.117193291793228, 0.00016520122725378975], 'sr_sf_all': [7.644690204747349, 0.0003323291686640057], 'sr_of_all': [12.343766141289052, 0.0004672543191486733], 'sr_all': [19.9884563460364, 0.0007995834878126791]}, 'Bkg': {'sr_4l_sf_A': [3.273017741052172, 0.045468003429151346], 'sr_4l_sf_B': [6.543019324501074, 0.1495065949817322], 'sr_4l_sf_C': [3.095189576022676, 0.03712707867259444], 'sr_4l_of_1': [1.6382048119266983, 0.06505721623246943], 'sr_4l_of_2': [1.0893333174499276, 0.01577919577925566], 'sr_4l_of_3': [1.3746751891558233, 0.03644128851540304], 'sr_4l_of_4': [4.587723127675417, 0.0813468286954406], 'sr_sf_all': [12.911226641575922, 0.23210167708347798], 'sr_of_all': [8.689936446207867, 0.19862452922256874], 'sr_all': [21.60116308778379, 0.43072620630604674]}, 'Zmetric': {'sr_4l_sf_A': [1.5271304041317981, None], 'sr_4l_sf_B': [1.211362434502514, None], 'sr_4l_sf_C': [0.6343417608521935, None], 'sr_4l_of_1': [2.171342366726005, None], 'sr_4l_of_2': [1.5478844248069015, None], 'sr_4l_of_3': [1.285097369234633, None], 'sr_4l_of_4': [2.0756701298826195, None], 'sr_sf_all': [2.0498574800479834, None], 'sr_of_all': [3.61532233730269, None], 'sr_all': [4.1560162765692406, None]}}

# Lep SFs top MVA (e.g. nov04_lepSFTight_CRttZ1l.pkl.gz)
EWK_REF = {'WWZ': {'sr_4l_sf_B': (1.74868254069881, 4.917508504242763e-05), 'sr_4l_sf_presel': (4.485213322035126, 0.0001274834371272613), 'all_events': (140.49884896972947, 0.0042351610608142515), '4l_presel': (23.282545120181545, 0.0006555966586360433), 'cr_4l_sf': (0.9351660713566458, 2.6186025483711817e-05), 'sr_4l_of_4': (4.622167446493546, 0.00013320708323522686), 'sr_4l_of_presel': (8.22255374678742, 0.00023250123569042185), 'cr_4l_of': (0.4185447171302684, 1.1512598140844779e-05), 'sr_4l_of_3': (1.3010239895778046, 3.616291039660801e-05), 'sr_4l_sf_A': (2.0713672516110853, 5.9315776828571735e-05), 'sr_4l_sf_C': (0.5017651491147224, 1.4345213220363006e-05), 'sr_4l_of_2': (0.6392456712632927, 1.7562125977397695e-05), 'sr_4l_of_1': (0.5547224892334729, 1.504718134311537e-05), 'sr_sf_all': [4.321814941424617, 0.00012283607509136236], 'sr_of_all': [7.117159596568117, 0.00020197930095234795], 'sr_all': [11.438974537992735, 0.0003248153760437103]}, 'ZH': {'sr_4l_sf_B': (1.2788690998095447, 7.608744194847245e-05), 'sr_4l_sf_presel': (2.8082857463062836, 0.00016194320312425535), 'all_events': (135.1955275044606, 0.00804450770608423), '4l_presel': (17.867805208533284, 0.0008812736600852035), 'cr_4l_sf': (0.05997741084194117, 2.476135332980253e-06), 'sr_4l_of_4': (0.12922648540757373, 8.987281304726395e-06), 'sr_4l_of_presel': (6.700939767087803, 0.0003149247951256664), 'cr_4l_of': (0.33829148431526557, 2.0302748603509765e-05), 'sr_4l_of_3': (0.29687829180933717, 1.4754310441059466e-05), 'sr_4l_sf_A': (0.7998454691442916, 5.547575575722742e-05), 'sr_4l_sf_C': (0.5617188863544369, 2.263170941530228e-05), 'sr_4l_of_2': (1.1344083935660307, 5.067890254114451e-05), 'sr_4l_of_1': (2.55994128557563, 0.00011380052967514463), 'sr_sf_all': [2.640433455308273, 0.00015419490712100214], 'sr_of_all': [4.120454456358572, 0.000188221023962075], 'sr_all': [6.7608879116668446, 0.00034241593108307717]}, 'ZZ': {'sr_4l_sf_B': (3.67848669195266, 0.0019192733711653476), 'sr_4l_sf_presel': (10.644875305946407, 0.005518689817435923), 'all_events': (20480.140682549172, 12.056295896759893), '4l_presel': (2720.7554111197423, 1.3958887964699902), 'cr_4l_sf': (1534.2599060524597, 0.787760489758707), 'sr_4l_of_4': (0.42259493833710116, 0.0002322511185316642), 'sr_4l_of_presel': (16.656698245005934, 0.008506865772611031), 'cr_4l_of': (0.720524285849299, 0.0003826803013148129), 'sr_4l_of_3': (0.3154593548227072, 0.00016727607892548892), 'sr_4l_sf_A': (1.0859978774794223, 0.0005727572547785321), 'sr_4l_sf_C': (2.1505048033201524, 0.0010757501569306244), 'sr_4l_of_2': (0.49655303148241503, 0.00025851613979802616), 'sr_4l_of_1': (0.5022183542971218, 0.0002573535069119564), 'sr_sf_all': [6.914989372752235, 0.0035677807828745044], 'sr_of_all': [1.7368256789393453, 0.0009153968441671357], 'sr_all': [8.65181505169158, 0.00448317762704164]}, 'ttZ': {'sr_4l_sf_B': (0.8623534921040487, 0.0026072132614722385), 'sr_4l_sf_presel': (2.2199135626257247, 0.006866412477031506), 'all_events': (6383.087986846862, 30.391735907350938), '4l_presel': (142.05677820627022, 0.4559282075895868), 'cr_4l_sf': (0.4439092328322823, 0.0013620522379399643), 'sr_4l_of_4': (1.8944468392895115, 0.0059636123355944455), 'sr_4l_of_presel': (3.4678790244953532, 0.010304685806516485), 'cr_4l_of': (53.768419315835885, 0.1731114682795787), 'sr_4l_of_3': (0.653330679680918, 0.001748962496855188), 'sr_4l_sf_A': (1.0055752564436247, 0.0032314494211373073), 'sr_4l_sf_C': (0.2161743800251384, 0.0006931474238017385), 'sr_4l_of_2': (0.24119583487607194, 0.0006679596224336652), 'sr_4l_of_1': (0.23368665546995349, 0.0006164622415794652), 'sr_sf_all': [2.0841031285728118, 0.006531810106411284], 'sr_of_all': [3.022660009316455, 0.008996996696462763], 'sr_all': [5.106763137889266, 0.015528806802874048]}, 'tWZ': {'sr_4l_sf_B': (0.2859664381300148, 0.00017667278493666024), 'sr_4l_sf_presel': (0.7152586082876218, 0.00044209172801275276), 'all_events': (455.4749377694133, 0.30381106357458043), '4l_presel': (19.656311394705504, 0.01211181234852329), 'cr_4l_sf': (0.14715859224582475, 9.008463241592307e-05), 'sr_4l_of_4': (0.6759082105012022, 0.0004220498538698098), 'sr_4l_of_presel': (1.2397287775225816, 0.0007629024552210294), 'cr_4l_of': (6.611493527295336, 0.00406995273997062), 'sr_4l_of_3': (0.1972787183086039, 0.00012040479187796226), 'sr_4l_sf_A': (0.3201974562306746, 0.00019807125613647463), 'sr_4l_sf_C': (0.07923669398284293, 4.903932175832338e-05), 'sr_4l_of_2': (0.11362721564154872, 6.739094537387002e-05), 'sr_4l_of_1': (0.09882170912487395, 5.9356989736279064e-05), 'sr_sf_all': [0.6854005883435323, 0.0004237833628314583], 'sr_of_all': [1.0856358535762287, 0.000669202580857921], 'sr_all': [1.7710364419197613, 0.0010929859436893793]}, 'other': {'sr_4l_sf_B': (1.1039758356047003, 0.1153498845316028), 'sr_4l_sf_presel': (2.175467418271284, 0.17904304942867552), 'all_events': (2676885.562879981, 4742821.485251619), '4l_presel': (47.720825052669255, 5.0050931433288754), 'cr_4l_sf': (10.618131090570161, 0.45610101045872087), 'sr_4l_of_4': (1.2560581627118375, 0.06480727751973193), 'sr_4l_of_presel': (3.424158982788239, 0.25620604445088613), 'cr_4l_of': (2.1181532489905788, 0.04497175017953798), 'sr_4l_of_3': (0.07670962255594513, 0.02764394673567354), 'sr_4l_sf_A': (0.587462637715449, 0.034215018598536415), 'sr_4l_sf_C': (0.3208864176982962, 0.02532342224375908), 'sr_4l_of_2': (0.11683187903784763, 0.011819071436898526), 'sr_4l_of_1': (0.6131084627263095, 0.04833481876350218), 'sr_sf_all': [2.0123248910184452, 0.17488832537389828], 'sr_of_all': [2.0627081270319394, 0.15260511445580618], 'sr_all': [4.0750330180503855, 0.3274934398297046]}, '$S/\\sqrt{B}$': {'sr_4l_sf_A': [1.6579073236239172, None], 'sr_4l_sf_B': [1.243184430438196, None], 'sr_4l_sf_C': [0.6393547723449631, None], 'sr_4l_of_1': [2.5885205017133126, None], 'sr_4l_of_2': [1.8025386843695004, None], 'sr_4l_of_3': [1.4333537035004176, None], 'sr_4l_of_4': [2.305033579902631, None], 'sr_sf_all': [2.1686260043470016, None], 'sr_of_all': [4.161402016178655, None], 'sr_all': [4.692569177645207, None]}, '$S/\\sqrt{S+B}$': {'sr_4l_sf_A': [1.185031303185296, None], 'sr_4l_sf_B': [1.0115280535666447, None], 'sr_4l_sf_C': [0.5433946453973049, None], 'sr_4l_of_1': [1.4581754520985202, None], 'sr_4l_of_2': [1.0711397102996822, None], 'sr_4l_of_3': [0.9480669557798651, None], 'sr_4l_of_4': [1.5837625997173428, None], 'sr_sf_all': [1.6500805838891106, None], 'sr_of_all': [2.584714849454142, None], 'sr_all': [3.0665154143288724, None]}, 'Sig': {'sr_4l_sf_A': [2.871212720755377, 0.00011479153258579916], 'sr_4l_sf_B': [3.0275516405083547, 0.0001252625269909001], 'sr_4l_sf_C': [1.0634840354691593, 3.697692263566529e-05], 'sr_4l_of_1': [3.114663774809103, 0.00012884771101826], 'sr_4l_of_2': [1.7736540648293233, 6.824102851854221e-05], 'sr_4l_of_3': [1.5979022813871417, 5.0917220837667476e-05], 'sr_4l_of_4': [4.75139393190112, 0.00014219436453995326], 'sr_sf_all': [6.962248396732891, 0.00027703098221236455], 'sr_of_all': [11.237614052926688, 0.0003902003249144229], 'sr_all': [18.199862449659577, 0.0006672313071267874]}, 'Bkg': {'sr_4l_sf_A': [2.999233227869171, 0.03821729653058873], 'sr_4l_sf_B': [5.930782457791424, 0.12005304394917705], 'sr_4l_sf_C': [2.76680229502643, 0.027141359146249766], 'sr_4l_of_1': [1.4478351816182586, 0.04926799150172988], 'sr_4l_of_2': [0.9682079610378834, 0.012812938144504089], 'sr_4l_of_3': [1.2427783753681743, 0.029680590103332176], 'sr_4l_of_4': [4.249008150839652, 0.07142519082772786], 'sr_sf_all': [11.696817980687026, 0.18541169962601556], 'sr_of_all': [7.907829668863968, 0.163186710577294], 'sr_all': [19.604647649550998, 0.34859841020330956]}, 'Zmetric': {'sr_4l_sf_A': [1.4637115657112032, None], 'sr_4l_sf_B': [1.1550907652954854, None], 'sr_4l_sf_C': [0.6038316021520365, None], 'sr_4l_of_1': [2.0601831062065683, None], 'sr_4l_of_2': [1.4700171813207716, None], 'sr_4l_of_3': [1.2251301329720645, None], 'sr_4l_of_4': [2.0020793789373705, None], 'sr_sf_all': [1.9599231687556398, None], 'sr_of_all': [3.45174892290265, None], 'sr_all': [3.9693663794344753, None]}}

# Lep SFs topMVA and btag corrections (e.g. nov19_lepSFTight_CRttZ1l_btagSFsingleWP_withUL16APVF_corlibBtagSF_newBtagEFFbd4162f.pkl.gz)
#EWK_REF = {'WWZ': {'sr_4l_sf_A': (1.9328528313343771, 5.213755600250673e-05), 'sr_4l_sf_B': (1.6453567291371223, 4.388333364131873e-05), 'sr_4l_sf_C': (0.4761728793852009, 1.2995780148850957e-05), 'sr_4l_of_1': (0.5232237891276427, 1.34827875241377e-05), 'sr_4l_of_2': (0.6034339692979168, 1.576638237951547e-05), 'sr_4l_of_3': (1.2253236255671685, 3.23119614163263e-05), 'sr_4l_of_4': (4.34639790601154, 0.00011872493338036402), 'sr_4l_sf_presel': (4.2059161040137925, 0.00011304061115484224), 'sr_4l_of_presel': (7.744366932601156, 0.00020780212769824956), 'all_events': (154.73522662092805, 0.005647558222033173), '4l_presel': (23.506105509185, 0.0006969176653045498), 'cr_4l_of': (2.719997349356926, 0.00010234799809708152), 'cr_4l_sf': (0.8826387009123812, 2.3492256993059274e-05), 'sr_sf_all': [4.054382439856701, 0.0001090166697926764], 'sr_of_all': [6.698379290004268, 0.0001802860647003435], 'sr_all': [10.752761729860968, 0.00028930273449301986]}, 'ZH': {'sr_4l_sf_A': (0.7488679706238699, 4.7678868181546e-05), 'sr_4l_sf_B': (1.211665390873752, 6.72379729411185e-05), 'sr_4l_sf_C': (0.5355616177276739, 2.0290587523138885e-05), 'sr_4l_of_1': (2.433313092376997, 0.00010041310383502371), 'sr_4l_of_2': (1.0788351300897419, 4.512231360582697e-05), 'sr_4l_of_3': (0.28174983232668194, 1.3033933078325302e-05), 'sr_4l_of_4': (0.12237884045258428, 7.932386116297714e-06), 'sr_4l_sf_presel': (2.6514237376142438, 0.00014179180443829814), 'sr_4l_of_presel': (6.367442568621571, 0.0002788470806316566), 'all_events': (148.35524375314753, 0.010281545149311684), '4l_presel': (18.01854265138008, 0.0009278156064695831), 'cr_4l_of': (1.918142139692452, 0.00014479361131858044), 'cr_4l_sf': (0.057239810557850615, 2.259675458538832e-06), 'sr_sf_all': [2.4960949792252958, 0.00013520742864580338], 'sr_of_all': [3.916276895246005, 0.0001665017366354737], 'sr_all': [6.4123718744713, 0.00030170916528127705]}, 'ZZ': {'sr_4l_sf_A': (1.016690630867433, 0.0005039304203131004), 'sr_4l_sf_B': (3.3825657961164852, 0.0016343260971195938), 'sr_4l_sf_C': (1.9903822665283606, 0.0009283558449268911), 'sr_4l_of_1': (0.4747548751517215, 0.00023153334037399474), 'sr_4l_of_2': (0.4699339390755087, 0.00023338714214114095), 'sr_4l_of_3': (0.2979325959045335, 0.0001499326592357303), 'sr_4l_of_4': (0.40183271849912555, 0.0002109213257393783), 'sr_4l_sf_presel': (9.90901400195508, 0.004818257468338145), 'sr_4l_of_presel': (15.942973111892883, 0.007826042075432008), 'all_events': (22510.162130849723, 15.858978389565964), '4l_presel': (2750.3427774826228, 1.4780420922210347), 'cr_4l_of': (4.194991716194073, 0.0029413111574327793), 'cr_4l_sf': (1469.077011626406, 0.7264159103039762), 'sr_sf_all': [6.389638693512278, 0.0030666123623595853], 'sr_of_all': [1.6444541286308891, 0.0008257744674902443], 'sr_all': [8.034092822143167, 0.00389238682984983]}, 'ttZ': {'sr_4l_sf_A': (1.0040904162898063, 0.003174164118648877), 'sr_4l_sf_B': (0.8644839242267756, 0.002627522256113313), 'sr_4l_sf_C': (0.21876828039843882, 0.0007184064251811959), 'sr_4l_of_1': (0.23622393602033473, 0.0006408656684392281), 'sr_4l_of_2': (0.24545506437456024, 0.0007084272257368942), 'sr_4l_of_3': (0.6670764706438969, 0.0018062633088032733), 'sr_4l_of_4': (1.898767666101388, 0.005953020451609651), 'sr_4l_sf_presel': (2.2294353861152296, 0.006871342423870653), 'sr_4l_of_presel': (3.5029047410250254, 0.010455439688981667), 'all_events': (6732.614511395827, 37.5506666045899), '4l_presel': (143.20231527299708, 0.48353752370665865), 'cr_4l_of': (58.53269960779905, 0.19862214416068574), 'cr_4l_sf': (0.43789838384843677, 0.0013672675495066703), 'sr_sf_all': [2.0873426209150208, 0.006520092799943387], 'sr_of_all': [3.04752313714018, 0.009108576654589046], 'sr_all': [5.1348657580552, 0.015628669454532434]}, 'tWZ': {'sr_4l_sf_A': (0.3151479764540492, 0.0001956475153507406), 'sr_4l_sf_B': (0.28262821736325916, 0.00017550439905985452), 'sr_4l_sf_C': (0.07862374972712925, 4.8976386938752056e-05), 'sr_4l_of_1': (0.09656565498504384, 5.7523984922841484e-05), 'sr_4l_of_2': (0.11151339673621262, 6.619464194934773e-05), 'sr_4l_of_3': (0.19378064624962915, 0.00011778900840292277), 'sr_4l_of_4': (0.6672360288263123, 0.0004176260567433958), 'sr_4l_sf_presel': (0.7054981850916262, 0.00043782367788507376), 'sr_4l_of_presel': (1.2242457705058638, 0.0007557616277358607), 'all_events': (481.4677892777096, 0.37290934855350727), '4l_presel': (19.80510777197338, 0.012832275215475131), 'cr_4l_of': (7.642521036504865, 0.004986257914654622), 'cr_4l_sf': (0.14568836262943768, 8.990692673261788e-05), 'sr_sf_all': [0.6763999435444376, 0.00042012830134934714], 'sr_of_all': [1.0690957267971979, 0.0006591336920185077], 'sr_all': [1.7454956703416356, 0.001079261993367855]}, 'other': {'sr_4l_sf_A': (0.5690752960390437, 0.03105833484186301), 'sr_4l_sf_B': (1.0741286377497536, 0.0927840461373301), 'sr_4l_sf_C': (0.29824033849416226, 0.023373605387380765), 'sr_4l_of_1': (0.5510847696019162, 0.03959363958073122), 'sr_4l_of_2': (0.09538458600770297, 0.011042215526065423), 'sr_4l_of_3': (0.08089019268582368, 0.023890108872888356), 'sr_4l_of_4': (1.1768197791978092, 0.05488650128467712), 'sr_4l_sf_presel': (2.0826611598475586, 0.15048406599050587), 'sr_4l_of_presel': (3.1709478761658714, 0.2187377941695375), 'all_events': (2860728.18394136, 6017867.169534823), '4l_presel': (47.88942795718206, 4.38712317244549), 'cr_4l_of': (3.417519193830976, 0.13915612452607556), 'cr_4l_sf': (9.825664693920267, 0.388140956221059), 'sr_sf_all': [1.9414442722829597, 0.14721598636657388], 'sr_of_all': [1.9041793274932521, 0.12941246526436212], 'sr_all': [3.8456235997762116, 0.276628451630936]}, '$S/\\sqrt{B}$': {'sr_4l_sf_A': [1.5734037088293156, None], 'sr_4l_sf_B': [1.2069020767748808, None], 'sr_4l_sf_C': [0.6291454631538783, None], 'sr_4l_of_1': [2.536488064796381, None], 'sr_4l_of_2': [1.751710879889635, None], 'sr_4l_of_3': [1.3535666345294313, None], 'sr_4l_of_4': [2.1950499282717715, None], 'sr_sf_all': [2.080393200260103, None], 'sr_of_all': [4.019035895746074, None], 'sr_all': [4.52555912556492, None]}, '$S/\\sqrt{S+B}$': {'sr_4l_sf_A': [1.1345794389997155, None], 'sr_4l_sf_B': [0.9822163334368799, None], 'sr_4l_sf_C': [0.5333976768387407, None], 'sr_4l_of_1': [1.4232609752341936, None], 'sr_4l_of_2': [1.0423860982624902, None], 'sr_4l_of_3': [0.9093364269383806, None], 'sr_4l_of_4': [1.5226515189827519, None], 'sr_sf_all': [1.592649525390964, None], 'sr_of_all': [2.5015197315018383, None], 'sr_all': [2.9654904278754786, None]}, 'Sig': {'sr_4l_sf_A': [2.681720801958247, 9.981642418405273e-05], 'sr_4l_sf_B': [2.8570221200108743, 0.00011112130658243723], 'sr_4l_sf_C': [1.0117344971128748, 3.328636767198984e-05], 'sr_4l_of_1': [2.95653688150464, 0.00011389589135916141], 'sr_4l_of_2': [1.6822690993876588, 6.088869598534244e-05], 'sr_4l_of_3': [1.5070734578938505, 4.53458944946516e-05], 'sr_4l_of_4': [4.468776746464124, 0.00012665731949666173], 'sr_sf_all': [6.550477419081996, 0.00024422409843847983], 'sr_of_all': [10.614656185250274, 0.0003467878013358172], 'sr_all': [17.16513360433227, 0.000591011899774297]}, 'Bkg': {'sr_4l_sf_A': [2.9050043196503323, 0.03493207689617573], 'sr_4l_sf_B': [5.603806575456273, 0.09722139888962286], 'sr_4l_sf_C': [2.586014635148091, 0.025069344044427604], 'sr_4l_of_1': [1.3586292357590164, 0.040523562574467285], 'sr_4l_of_2': [0.9222869861939845, 0.012050224535892806], 'sr_4l_of_3': [1.2396799054838832, 0.02596409384933028], 'sr_4l_of_4': [4.144656192624635, 0.06146806911876954], 'sr_sf_all': [11.094825530254697, 0.15722281983022618], 'sr_of_all': [7.665252320061519, 0.14000595007845992], 'sr_all': [18.760077850316215, 0.29722876990868613]}, 'Zmetric': {'sr_4l_sf_A': [1.3941054315112151, None], 'sr_4l_sf_B': [1.121470334733581, None], 'sr_4l_sf_C': [0.593646320326605, None], 'sr_4l_of_1': [2.015105967190228, None], 'sr_4l_of_2': [1.4294612189280722, None], 'sr_4l_of_3': [1.1646040222160488, None], 'sr_4l_of_4': [1.9141453826105042, None], 'sr_sf_all': [1.88511050591042, None], 'sr_of_all': [3.3353060593828427, None], 'sr_all': [3.831175815497084, None]}}

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
    #years = ["UL16"]
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
            var = sum(sum(histos_dict[dense_axis][{"category":cat_name,"process":sample_dict[proc_name]}].variances(flow=True)))
            yld_dict[proc_name][cat_name] = (val,var)

    # Print to screen
    if not quiet:
        for proc in yld_dict.keys():
            print(f"\n{proc}:")
            for cat in yld_dict[proc].keys():
                val = yld_dict[proc][cat]
                print(f"\t{cat}: {val}")

    return yld_dict


# Gets the process sums for S and B and gets metrics e.g. S/sqrt(B) and puts it into the dict
# Hard coded for the summed values (e.g. looking for "ZH" not "GluGluZH","qqToZHToZTo2L")
def put_proc_row_sums(yld_dict):
    sig_lst = ["WWZ","ZH"]
    bkg_lst = ["ZZ","ttZ","tWZ","other"]
    sig_sum = {"sr_4l_sf_A":[0,0], "sr_4l_sf_B":[0,0], "sr_4l_sf_C":[0,0], "sr_4l_of_1":[0,0], "sr_4l_of_2":[0,0], "sr_4l_of_3":[0,0], "sr_4l_of_4":[0,0]}
    bkg_sum = {"sr_4l_sf_A":[0,0], "sr_4l_sf_B":[0,0], "sr_4l_sf_C":[0,0], "sr_4l_of_1":[0,0], "sr_4l_of_2":[0,0], "sr_4l_of_3":[0,0], "sr_4l_of_4":[0,0]}
    for proc in yld_dict.keys():
        print(proc)
        for cat in yld_dict[proc].keys():
            if cat not in sig_sum: continue
            val,var = yld_dict[proc][cat]
            print("   ",cat,val)
            if proc in sig_lst:
                sig_sum[cat][0] += val
                sig_sum[cat][1] += var
            if proc in bkg_lst:
                bkg_sum[cat][0] += val
                bkg_sum[cat][1] += var

    yld_dict[SOVERROOTB] = {}
    yld_dict[SOVERROOTSPLUSB] = {}
    yld_dict["Sig"] = {}
    yld_dict["Bkg"] = {}
    yld_dict["Zmetric"] = {}
    for cat in sig_sum.keys():
        s = sig_sum[cat][0]
        b = bkg_sum[cat][0]
        s_var = sig_sum[cat][1]
        b_var = bkg_sum[cat][1]
        yld_dict[SOVERROOTB][cat]      = [s/math.sqrt(b) , None]
        yld_dict[SOVERROOTSPLUSB][cat] = [s/math.sqrt(s+b) , None]
        yld_dict["Zmetric"][cat] = [math.sqrt(2 * ((s + b) * math.log(1 + s / b) - s)), None] # Eq 18 https://cds.cern.ch/record/2203244/files/1087459_109-114.pdf
        yld_dict["Sig"][cat] = [s, s_var]
        yld_dict["Bkg"][cat] = [b, b_var]

# Gets the sums of categoreis (assumed to be columns in the input dict) and puts them into the dict
# Special handling for rows that are metrics (e.g. s/sqrt(b)), sums these in quadrature
def put_cat_col_sums(yld_dict,metrics_names_lst=["Zmetric",SOVERROOTB,SOVERROOTSPLUSB]):

    # Hard coded names we expect the columns to be
    sr_sf_lst = ["sr_4l_sf_A","sr_4l_sf_B","sr_4l_sf_C"]
    sr_of_lst = ["sr_4l_of_1","sr_4l_of_2","sr_4l_of_3","sr_4l_of_4"]
    sr_lst = sr_sf_lst + sr_of_lst

    # Loop over rows (processes) and sum columns together, fill the result into new_dict
    new_dict = {}
    for proc in yld_dict:
        sr_sf_val = 0
        sr_sf_var = 0
        sr_of_val = 0
        sr_of_var = 0
        sr_val = 0
        sr_var = 0
        for cat in yld_dict[proc]:
            val = yld_dict[proc][cat][0]
            var = yld_dict[proc][cat][1]
            if cat in sr_sf_lst:
                if proc in metrics_names_lst:
                    sr_sf_val += val*val
                else:
                    sr_sf_val += val
                    sr_sf_var += var
            if cat in sr_of_lst:
                if proc in metrics_names_lst:
                    sr_of_val += val*val
                else:
                    sr_of_val += val
                    sr_of_var += var
            if cat in sr_lst:
                if proc in metrics_names_lst:
                    sr_val += val*val
                else:
                    sr_val += val
                    sr_var += var

        # Fill our new_dict with what we've computed
        new_dict[proc] = {}
        if proc in metrics_names_lst:
            new_dict[proc]["sr_sf_all"] = [np.sqrt(sr_sf_val),None]
            new_dict[proc]["sr_of_all"] = [np.sqrt(sr_of_val),None]
            new_dict[proc]["sr_all"]    = [np.sqrt(sr_val),None]
        else:
            new_dict[proc]["sr_sf_all"] = [sr_sf_val, sr_sf_var]
            new_dict[proc]["sr_of_all"] = [sr_of_val, sr_of_var]
            new_dict[proc]["sr_all"]    = [sr_val, sr_var]

    # Put the columns into the yld_dict
    for proc in new_dict:
        yld_dict[proc]["sr_sf_all"] = new_dict[proc]["sr_sf_all"]
        yld_dict[proc]["sr_of_all"] = new_dict[proc]["sr_of_all"]
        yld_dict[proc]["sr_all"] = new_dict[proc]["sr_all"]


# Print yields
def print_yields(yld_dict_in,print_fom=True):

    # Get err from var
    def get_err_from_var(in_dict):
        out_dict = {}
        for proc in in_dict:
            out_dict[proc] = {}
            for cat in in_dict[proc]:
                if in_dict[proc][cat][1] is None: var = None
                else: var = np.sqrt(in_dict[proc][cat][1])
                out_dict[proc][cat] = [in_dict[proc][cat][0],var]
        return out_dict

    yld_dict = get_err_from_var(yld_dict_in)

    # Dump the yields to dict for latex table
    yld_dict_for_printing = {}
    cats_to_print = [ "sr_4l_sf_A", "sr_4l_sf_B", "sr_4l_sf_C", "sr_sf_all", "sr_4l_of_1", "sr_4l_of_2", "sr_4l_of_3", "sr_4l_of_4", "sr_of_all", "sr_all"]
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
        print_errs=False, # If you want to turn this on, figure out where/when you want to sqrt the var
        column_variable="subkeys",
    )
    #exit()


    # Compare with other yields, print comparison

    #tag1 = "ewkcoffea"
    #tag2 = "VVVNanoLooper"
    tag1 = "lepTightSF btagSFtopeft"
    tag2 = "lepTightSF " # Hard coded ref

    #yld_dict_comp = utils.put_none_errs(KEEGAN_YIELDS)
    yld_dict_comp = get_err_from_var(EWK_REF)

    yld_dict_1 = yld_dict_for_printing
    yld_dict_2 = yld_dict_comp

    pdiff_dict = utils.get_diff_between_nested_dicts(yld_dict_1,yld_dict_2,difftype="percent_diff",inpercent=True)
    diff_dict  = utils.get_diff_between_nested_dicts(yld_dict_1,yld_dict_2,difftype="absolute_diff")

    procs_to_print = ["WWZ","ZH","Sig","ZZ","ttZ","tWZ","other","Bkg",SOVERROOTB,SOVERROOTSPLUSB,"Zmetric"]
    hlines = [1,2,6,7]
    #vlines = "cccc|c|cccc|c|c" # Not impliemnted in MLT

    mlt.print_begin()
    mlt.print_latex_yield_table(yld_dict_1,key_order=procs_to_print,subkey_order=cats_to_print,tag=tag1,hz_line_lst=hlines,print_errs=True,small=True)
    mlt.print_latex_yield_table(yld_dict_2,key_order=procs_to_print,subkey_order=cats_to_print,tag=tag2,hz_line_lst=hlines,print_errs=True,small=True)
    mlt.print_latex_yield_table(pdiff_dict,key_order=procs_to_print,subkey_order=cats_to_print,tag=f"Percent diff between {tag1} and {tag2}",hz_line_lst=hlines,small=True)
    mlt.print_latex_yield_table(diff_dict, key_order=procs_to_print,subkey_order=cats_to_print,tag=f"Diff between {tag1} and {tag2}",hz_line_lst=hlines,small=True)
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
        if var_name == "nbtagsm": continue # TMP
        #if var_name not in BDT_INPUT_LST and "bdt" not in var_name: continue # TMP
        if var_name not in TMP_LST: continue # TMP
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
            #if cat_name not in ["cr_4l_sf","cr_4l_of","sr_4l_of_presel","sr_4l_sf_presel"]: continue # TMP
            print(cat_name)

            histo_cat = histo[{"category":cat_name}]

            # Group the mc samples
            grouping_mc = sample_dict
            histo_grouped_mc = group(histo_cat,"process","process_grp",grouping_mc)

            # Group the data samples
            grouping_data = {'data': ["UL16APV_data","UL16_data","UL17_data","UL18_data"]}
            #grouping_data = {'data': ["UL16APV_data"]}
            histo_grouped_data = group(histo_cat,"process","process_grp",grouping_data)

            ######
            #if cat_name == "cr_4l_sf" and var_name == "nleps":
            #    print("mc\n",histo_grouped_mc)
            #    print("data\n",histo_grouped_data)
            #    print("val mc\n",sum(histo_grouped_mc.values(flow=True)))
            #    print("val data\n",sum(histo_grouped_data.values(flow=True)))
            #    #print("var mc\n",(histo_grouped_mc.variances(flow=True)))
            #    #print("var data\n",(histo_grouped_data.variances(flow=True)))
            #continue
            #print("\Yields")
            #print(type(histo_cat.values(flow=True)))
            #print(sum(histo_cat.values(flow=True)))
            #print(sum(sum(histo_cat.values(flow=True))))
            #exit()
            #####

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
        put_proc_row_sums(yld_dict)
        put_cat_col_sums(yld_dict)
        #print(yld_dict)
        #exit()
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

