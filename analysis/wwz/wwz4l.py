#!/usr/bin/env python
#import sys
import math
import coffea
import numpy as np
import awkward as ak
import copy
np.seterr(divide='ignore', invalid='ignore', over='ignore')
from coffea import processor
import hist
from hist import axis
from coffea.analysis_tools import PackedSelection
from coffea.lumi_tools import LumiMask

from topcoffea.modules.paths import topcoffea_path
import topcoffea.modules.event_selection as es_tc
import topcoffea.modules.object_selection as os_tc
import topcoffea.modules.corrections as cor_tc

from ewkcoffea.modules.paths import ewkcoffea_path as ewkcoffea_path
import ewkcoffea.modules.selection_wwz as es_ec
import ewkcoffea.modules.objects_wwz as os_ec
import ewkcoffea.modules.corrections as cor_ec

from topcoffea.modules.get_param_from_jsons import GetParam
get_tc_param = GetParam(topcoffea_path("params/params.json"))
get_ec_param = GetParam(ewkcoffea_path("params/params.json"))

import hist.dask as hda


# Small helper function for creating the list of systematics
# Append "Up" and "Down" to all base strings in a given syst list
def append_up_down_to_sys_base(sys_lst_in):
    sys_lst_out = []
    for s in sys_lst_in:
        sys_lst_out.append(f"{s}Up")
        sys_lst_out.append(f"{s}Down")
    return sys_lst_out


from memory_profiler import profile

import objgraph

import sys

class AnalysisProcessor(processor.ProcessorABC):

    def __init__(self, samples, wc_names_lst=[], hist_lst=None, ecut_threshold=None, do_errors=False, do_systematics=False, split_by_lepton_flavor=False, skip_signal_regions=False, skip_control_regions=False, muonSyst='nominal', dtype=np.float32):

        self._samples = samples
        self._wc_names_lst = wc_names_lst
        self._dtype = dtype

        # Create the dense axes for the histograms
        self._dense_axes_dict = {
            "mt2"   : axis.Regular(180, 0, 100, name="mt2",  label="mt2"),
            "met"   : axis.Regular(180, 0, 300, name="met",  label="met"),
            "metphi": axis.Regular(180, -4, 4, name="metphi", label="met phi"),
            "ptl4"  : axis.Regular(180, 0, 500, name="ptl4", label="ptl4"),
            "scalarptsum_lep" : axis.Regular(180, 0, 500, name="scalarptsum_lep", label="S_T"),
            "scalarptsum_lepmet" : axis.Regular(180, 0, 600, name="scalarptsum_lepmet", label="S_T + metpt"),
            "scalarptsum_lepmetjet" : axis.Regular(180, 0, 1100, name="scalarptsum_lepmetjet", label="S_T + metpt + H_T"),
            "mll_01": axis.Regular(180, 0, 200, name="mll_01",  label="mll_l0_l1"),
            "mllll": axis.Regular(180, 0, 600, name="mllll",  label="mllll"),
            "l0pt"  : axis.Regular(180, 0, 500, name="l0pt", label="l0pt"),
            "j0pt"  : axis.Regular(180, 0, 500, name="j0pt", label="j0pt"),

            "w_lep0_pt"  : axis.Regular(180, 0, 300, name="w_lep0_pt", label="Leading W lep pt"),
            "w_lep1_pt"  : axis.Regular(180, 0, 300, name="w_lep1_pt", label="Subleading W lep pt"),
            "z_lep0_pt"  : axis.Regular(180, 0, 300, name="z_lep0_pt", label="Leading Z lep pt"),
            "z_lep1_pt"  : axis.Regular(180, 0, 300, name="z_lep1_pt", label="Subleading Z lep pt"),
            "w_lep0_eta" : axis.Regular(180, -3, 3, name="w_lep0_eta", label="Leading W lep eta"),
            "w_lep1_eta" : axis.Regular(180, -3, 3, name="w_lep1_eta", label="Subleading W lep eta"),
            "z_lep0_eta" : axis.Regular(180, -3, 3, name="z_lep0_eta", label="Leading Z lep eta"),
            "z_lep1_eta" : axis.Regular(180, -3, 3, name="z_lep1_eta", label="Subleading Z lep eta"),
            "w_lep0_phi" : axis.Regular(180, -4, 4, name="w_lep0_phi", label="Leading W lep phi"),
            "w_lep1_phi" : axis.Regular(180, -4, 4, name="w_lep1_phi", label="Subleading W lep phi"),
            "z_lep0_phi" : axis.Regular(180, -4, 4, name="z_lep0_phi", label="Leading Z lep phi"),
            "z_lep1_phi" : axis.Regular(180, -4, 4, name="z_lep1_phi", label="Subleading Z lep phi"),
            "mll_wl0_wl1" : axis.Regular(180, 0, 200, name="mll_wl0_wl1", label="mll(W lep0, W lep1)"),
            "mll_zl0_zl1" : axis.Regular(180, 0, 200, name="mll_zl0_zl1", label="mll(Z lep0, Z lep1)"),

            "pt_zl0_zl1" : axis.Regular(180, 0, 300, name="pt_zl0_zl1", label="pt(Zl0 + Zl1)"),
            "pt_wl0_wl1" : axis.Regular(180, 0, 300, name="pt_wl0_wl1", label="pt(Wl0 + Wl1)"),
            "dr_zl0_zl1" : axis.Regular(180, 0, 5, name="dr_zl0_zl1", label="dr(Zl0,Zl1)"),
            "dr_wl0_wl1" : axis.Regular(180, 0, 5, name="dr_wl0_wl1", label="dr(Wl0,Wl1)"),
            "dr_wleps_zleps" : axis.Regular(180, 0, 5, name="dr_wleps_zleps", label="dr((Wl0+Wl1),(Zl0,Zl1))"),

            "absdphi_zl0_zl1" : axis.Regular(180, 0, 4, name="absdphi_zl0_zl1", label="abs dphi(Zl0,Zl1)"),
            "absdphi_wl0_wl1" : axis.Regular(180, 0, 4, name="absdphi_wl0_wl1", label="abs dphi(Wl0,Wl1)"),
            "absdphi_z_ww"    : axis.Regular(180, 0, 4, name="absdphi_z_ww", label="abs dphi((Zl0+Zl1),(Wl0+Wl1+met))"),
            "dphi_4l_met"    : axis.Regular(180, -4, 4, name="dphi_4l_met", label="dphi((Zl0+Zl1+Wl0+Wl1),met)"),
            "dphi_zleps_met" : axis.Regular(180, -4, 4, name="dphi_zleps_met", label="dphi((Zl0+Zl1),met)"),
            "dphi_wleps_met" : axis.Regular(180, -4, 4, name="dphi_wleps_met", label="dphi((Wl0+Wl1),met)"),

            "absdphi_min_afas" : axis.Regular(180, 0, 4, name="absdphi_min_afas",  label="min(abs(delta phi of all pairs))"),
            "absdphi_min_afos" : axis.Regular(180, 0, 4, name="absdphi_min_afos",  label="min(abs(delta phi of OS pairs))"),
            "absdphi_min_sfos" : axis.Regular(180, 0, 4, name="absdphi_min_sfos",  label="min(abs(delta phi of SFOS pairs))"),
            "mll_min_afas" : axis.Regular(180, 0, 150, name="mll_min_afas",  label="min mll of all pairs"),
            "mll_min_afos" : axis.Regular(180, 0, 150, name="mll_min_afos",  label="min mll of OF pairs"),
            "mll_min_sfos" : axis.Regular(180, 0, 150, name="mll_min_sfos",  label="min mll of SFOF pairs"),

            "mlb_min" : axis.Regular(180, 0, 300, name="mlb_min",  label="min mass(b+l)"),
            "mlb_max" : axis.Regular(180, 0, 500, name="mlb_max",  label="max mass(b+l)"),

            "njets"   : axis.Regular(8, 0, 8, name="njets",   label="Jet multiplicity"),
            "nleps"   : axis.Regular(5, 0, 5, name="nleps",   label="Lep multiplicity"),
            "nbtagsl" : axis.Regular(6, 0, 6, name="nbtagsl", label="Loose btag multiplicity"),
            "nbtagsm" : axis.Regular(4, 0, 4, name="nbtagsm", label="Medium btag multiplicity"),

            "njets_counts"   : axis.Regular(30, 0, 30, name="njets_counts",   label="Jet multiplicity counts"),
            "nleps_counts"   : axis.Regular(30, 0, 30, name="nleps_counts",   label="Lep multiplicity counts"),
            "nbtagsl_counts" : axis.Regular(30, 0, 30, name="nbtagsl_counts", label="Loose btag multiplicity counts"),

            "bdt_of_wwz_raw": axis.Regular(180, -3.5, 3.5, name="bdt_of_wwz_raw", label="Raw score bdt_of_wwz"),
            "bdt_sf_wwz_raw": axis.Regular(180, -3.5, 3.5, name="bdt_sf_wwz_raw", label="Raw score bdt_sf_wwz"),
            "bdt_of_zh_raw" : axis.Regular(180, -3.5, 3.5, name="bdt_of_zh_raw", label="Raw score bdt_of_zh"),
            "bdt_sf_zh_raw" : axis.Regular(180, -3.5, 3.5, name="bdt_sf_zh_raw", label="Raw score bdt_sf_zh"),
            "bdt_of_wwz": axis.Regular(180, -1, 1, name="bdt_of_wwz", label="Score bdt_of_wwz"),
            "bdt_sf_wwz": axis.Regular(180, -1, 1, name="bdt_sf_wwz", label="Score bdt_sf_wwz"),
            "bdt_of_zh" : axis.Regular(180, -1, 1, name="bdt_of_zh", label="Score bdt_of_zh"),
            "bdt_sf_zh" : axis.Regular(180, -1, 1, name="bdt_sf_zh", label="Score bdt_sf_zh"),

        }

        # Set the list of hists to fill
        if hist_lst is None:
            # If the hist list is none, assume we want to fill all hists
            self._hist_lst = list(self._dense_axes_dict.keys())
        else:
            # Otherwise, just fill the specified subset of hists
            for hist_to_include in hist_lst:
                if hist_to_include not in self._dense_axes_dict.keys():
                    raise Exception(f"Error: Cannot specify hist \"{hist_to_include}\", it is not defined in the processor.")
            self._hist_lst = hist_lst # Which hists to fill

        # Set the energy threshold to cut on
        self._ecut_threshold = ecut_threshold

        # Set the booleans
        self._do_errors = do_errors # Whether to calculate and store the w**2 coefficients
        self._do_systematics = do_systematics # Whether to process systematic samples
        self._split_by_lepton_flavor = split_by_lepton_flavor # Whether to keep track of lepton flavors individually
        self._skip_signal_regions = skip_signal_regions # Whether to skip the SR categories
        self._skip_control_regions = skip_control_regions # Whether to skip the CR categories


    # Main function: run on a given dataset
    #@profile
    def process(self, events):

        # Loop over samples and fill histos
        #hout = {}
        #for events in events_dict.values():

        # Dataset parameters
        dataset = events.metadata["dataset"]

        # If we pass the root files instead of events object
        #from coffea.nanoevents import NanoEventsFactory
        #from coffea.nanoevents import NanoAODSchema
        #events = NanoEventsFactory.from_root(
        #    {fpath: "/Events" for fpath in fpaths},
        #    schemaclass=NanoAODSchema,
        #    metadata={"dataset": dataset},
        #).events()

        isData             = self._samples[dataset]["isData"]
        histAxisName       = self._samples[dataset]["histAxisName"]
        year               = self._samples[dataset]["year"]
        xsec               = self._samples[dataset]["xsec"]
        sow                = self._samples[dataset]["nSumOfWeights"]

        # Get up down weights from input dict
        if (self._do_systematics and not isData):
            lhe_sow = self._samples[dataset]["nSumOfLheWeights"]
            # This assumes we have an NLO xsec, so for these systs we will have e.g. xsec_NLO*(N_pass_up/N_gen_up)
            # Thus these systs should only affect acceptance and effeciency and shape
            # The uncty on xsec comes from NLO and is applied as a rate uncty in the text datacard
            if lhe_sow == []:
                sow_renormDown     = sow
                sow_factDown       = sow
                sow_factUp         = sow
                sow_renormUp       = sow
            elif len(lhe_sow) == 9:
                sow_renormDown     = lhe_sow[1]
                sow_factDown       = lhe_sow[3]
                sow_factUp         = lhe_sow[5]
                sow_renormUp       = lhe_sow[7]
            elif len(lhe_sow) == 8:
                sow_renormDown     = lhe_sow[1]
                sow_factDown       = lhe_sow[3]
                sow_factUp         = lhe_sow[4]
                sow_renormUp       = lhe_sow[6]
            else: raise Exception("ERROR: Unknown LHE weights length {len(lhe_sow)}")
        else:
            sow_renormUp       = -1
            sow_renormDown     = -1
            sow_factUp         = -1
            sow_factDown       = -1


        datasets = ["SingleMuon", "SingleElectron", "EGamma", "MuonEG", "DoubleMuon", "DoubleElectron", "DoubleEG"]
        for d in datasets:
            if d in dataset: dataset = dataset.split('_')[0]

        # Initialize objects
        #met  = events.MET
        met  = events.PuppiMET
        ele  = events.Electron
        mu   = events.Muon
        tau  = events.Tau
        jets = events.Jet

        # An array of lenght events that is just 1 for each event
        # Probably there's a better way to do this, but we use this method elsewhere so I guess why not..
        events.nom = ak.ones_like(met.pt)

        # Get the lumi mask for data
        if year == "2016" or year == "2016APV":
            golden_json_path = topcoffea_path("data/goldenJsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt")
        elif year == "2017":
            golden_json_path = topcoffea_path("data/goldenJsons/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt")
        elif year == "2018":
            golden_json_path = topcoffea_path("data/goldenJsons/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt")
        else:
            raise ValueError(f"Error: Unknown year \"{year}\".")
        lumi_mask = LumiMask(golden_json_path)(events.run,events.luminosityBlock)


        ################### Lepton selection ####################

        # Do the object selection for the WWZ eleectrons
        ele_presl_mask = os_ec.is_presel_wwz_ele(ele,tight=True)
        ele["topmva"] = os_ec.get_topmva_score_ele(events, year)
        ele["is_tight_lep_for_wwz"] = ((ele.topmva > get_tc_param("topmva_wp_t_e")) & ele_presl_mask)

        # Do the object selection for the WWZ muons
        mu_presl_mask = os_ec.is_presel_wwz_mu(mu)
        mu["topmva"] = os_ec.get_topmva_score_mu(events, year)
        mu["is_tight_lep_for_wwz"] = ((mu.topmva > get_tc_param("topmva_wp_t_m")) & mu_presl_mask)

        # Get tight leptons for WWZ selection
        ele_wwz_t = ele[ele.is_tight_lep_for_wwz]
        mu_wwz_t = mu[mu.is_tight_lep_for_wwz]

        # Attach the lepton SFs to the electron and muons collections
        cor_ec.AttachElectronSF(ele_wwz_t,year=year)
        cor_ec.AttachMuonSF(mu_wwz_t,year=year)

        l_wwz_t = ak.with_name(ak.concatenate([ele_wwz_t,mu_wwz_t],axis=1),'PtEtaPhiMCandidate')
        l_wwz_t = l_wwz_t[ak.argsort(l_wwz_t.pt, axis=-1,ascending=False)] # Sort by pt


        # For WWZ: Compute pair invariant masses
        llpairs_wwz = ak.combinations(l_wwz_t, 2, fields=["l0","l1"])
        os_pairs_mask = (llpairs_wwz.l0.pdgId*llpairs_wwz.l1.pdgId < 0)   # Maks for opposite-sign pairs
        sfos_pairs_mask = (llpairs_wwz.l0.pdgId == -llpairs_wwz.l1.pdgId) # Mask for same-flavor-opposite-sign pairs
        ll_absdphi_pairs = abs(llpairs_wwz.l0.delta_phi(llpairs_wwz.l1))
        ll_mass_pairs = (llpairs_wwz.l0+llpairs_wwz.l1).mass            # The mll for each ll pair
        absdphi_min_afas = ak.min(ll_absdphi_pairs,axis=-1)
        absdphi_min_afos = ak.min(ll_absdphi_pairs[os_pairs_mask],axis=-1)
        absdphi_min_sfos = ak.min(ll_absdphi_pairs[sfos_pairs_mask],axis=-1)
        mll_min_afas = ak.min(ll_mass_pairs,axis=-1)
        mll_min_afos = ak.min(ll_mass_pairs[os_pairs_mask],axis=-1)
        mll_min_sfos = ak.min(ll_mass_pairs[sfos_pairs_mask],axis=-1)
        events["min_mll_afos"] = mll_min_afos # Attach this one to the event info since we need it for selection

        # For WWZ
        l_wwz_t_padded = ak.pad_none(l_wwz_t, 4)
        l0 = l_wwz_t_padded[:,0]
        l1 = l_wwz_t_padded[:,1]
        l2 = l_wwz_t_padded[:,2]
        l3 = l_wwz_t_padded[:,3]

        nleps = ak.num(l_wwz_t)

        # Put njets and l_fo_conept_sorted into events and get 4l event selection mask
        events["l_wwz_t"] = l_wwz_t
        es_ec.add4lmask_wwz(events, year, isData, histAxisName)


        ######### Normalization and weights ###########


        # These weights can go outside of the outside sys loop since they do not depend on pt of mu or jets
        # We only calculate these values if not isData
        # Note: add() will generally modify up/down weights, so if these are needed for any reason after this point, we should instead pass copies to add()
        # Note: Here we will to the weights object the SFs that do not depend on any of the forthcoming loops
        weights_obj_base = coffea.analysis_tools.Weights(None,storeIndividual=True)
        if not isData:
            genw = events["genWeight"]

            # If it's an EFT sample, just take SM piece
            sm_wgt = 1.0
            eft_coeffs = ak.to_numpy(events["EFTfitCoefficients"]) if hasattr(events, "EFTfitCoefficients") else None
            if eft_coeffs is not None:
                sm_wgt = eft_coeffs[:,0]

            # Normalize by (xsec/sow)*genw where genw is 1 for EFT samples
            # Note that for theory systs, will need to multiply by sow/sow_wgtUP to get (xsec/sow_wgtUp)*genw and same for Down
            lumi = 1000.0*get_tc_param(f"lumi_{year}")
            weights_obj_base.add("norm",(xsec/sow)*genw*lumi*sm_wgt)


            # Scale weights
            cor_tc.AttachPSWeights(events)
            cor_tc.AttachScaleWeights(events)
            # FSR/ISR weights
            # For now only consider variations in the numerator
            weights_obj_base.add('ISR', events.nom, events.ISRUp, events.ISRDown)
            weights_obj_base.add('FSR', events.nom, events.FSRUp, events.FSRDown)
            # Renorm/fact scale
            weights_obj_base.add('renorm', events.nom, events.renormUp*(sow/sow_renormUp), events.renormDown*(sow/sow_renormDown))
            weights_obj_base.add('fact', events.nom, events.factUp*(sow/sow_factUp), events.factDown*(sow/sow_factDown))

            # Misc other experimental SFs and systs
            weights_obj_base.add('PreFiring', events.L1PreFiringWeight.Nom,  events.L1PreFiringWeight.Up,  events.L1PreFiringWeight.Dn)
            weights_obj_base.add('PU', cor_tc.GetPUSF((events.Pileup.nTrueInt), year), cor_tc.GetPUSF(events.Pileup.nTrueInt, year, 'up'), cor_tc.GetPUSF(events.Pileup.nTrueInt, year, 'down'))

            # Lepton SFs and systs
            weights_obj_base.add("lepSF_muon", events.sf_4l_muon, copy.copy(events.sf_4l_hi_muon), copy.copy(events.sf_4l_lo_muon))
            weights_obj_base.add("lepSF_elec", events.sf_4l_elec, copy.copy(events.sf_4l_hi_elec), copy.copy(events.sf_4l_lo_elec))


        # Set up the list of systematics that are handled via event weight variations
        wgt_correction_syst_lst = [
            "btagSFlight_correlated", "btagSFbc_correlated", f"btagSFlight_uncorrelated_{year}", f"btagSFbc_uncorrelated_{year}",
            "lepSF_elec", "lepSF_muon", "PreFiring", "PU",
            "renorm", "fact", "ISR", "FSR",
        ]
        wgt_correction_syst_lst = append_up_down_to_sys_base(wgt_correction_syst_lst)


        ######### The rest of the processor is inside this loop over systs that affect object kinematics  ###########

        obj_correction_systs = [] # Will have e.g. jes etc

        # If we're doing systematics and this isn't data, we will loop over the obj correction syst lst list
        if self._do_systematics and not isData: obj_corr_syst_var_list = ["nominal"] + obj_correction_systs
        # Otherwise loop juse once, for nominal
        else: obj_corr_syst_var_list = ['nominal']

        # Loop over the list of systematic variations (that impact object kinematics) that we've constructed
        for obj_corr_syst_var in obj_corr_syst_var_list:
            # Make a copy of the base weights object, so that each time through the loop we do not double count systs
            # In this loop over systs that impact kinematics, we will add to the weights objects the SFs that depend on the object kinematics
            weights_obj_base_for_kinematic_syst = copy.copy(weights_obj_base) # TODO do we need copy here?


            #################### Jets ####################

            # Clean with dr (though another option is to use jetIdx)
            cleanedJets = os_ec.get_cleaned_collection(l_wwz_t,jets)

            # Selecting jets and cleaning them
            jetptname = "pt_nom" if hasattr(cleanedJets, "pt_nom") else "pt"
            cleanedJets["is_good"] = os_tc.is_tight_jet(getattr(cleanedJets, jetptname), cleanedJets.eta, cleanedJets.jetId, pt_cut=20., eta_cut=get_ec_param("wwz_eta_j_cut"), id_cut=get_ec_param("wwz_jet_id_cut"))
            goodJets = cleanedJets[cleanedJets.is_good]

            # Count jets
            njets = ak.num(goodJets)
            ht = ak.sum(goodJets.pt,axis=-1)
            j0 = goodJets[ak.argmax(goodJets.pt,axis=-1,keepdims=True)]


            # Loose DeepJet WP
            btagger = "btag" # For deep flavor WPs
            #btagger = "btagcsv" # For deep CSV WPs
            if year == "2017":
                btagwpl = get_tc_param(f"{btagger}_wp_loose_UL17")
                btagwpm = get_tc_param(f"{btagger}_wp_medium_UL17")
            elif year == "2018":
                btagwpl = get_tc_param(f"{btagger}_wp_loose_UL18")
                btagwpm = get_tc_param(f"{btagger}_wp_medium_UL18")
            elif year=="2016":
                btagwpl = get_tc_param(f"{btagger}_wp_loose_UL16")
                btagwpm = get_tc_param(f"{btagger}_wp_medium_UL16")
            elif year=="2016APV":
                btagwpl = get_tc_param(f"{btagger}_wp_loose_UL16APV")
                btagwpm = get_tc_param(f"{btagger}_wp_medium_UL16APV")
            else:
                raise ValueError(f"Error: Unknown year \"{year}\".")

            if btagger == "btag":
                isBtagJetsLoose = (goodJets.btagDeepFlavB > btagwpl)
                isBtagJetsMedium = (goodJets.btagDeepFlavB > btagwpm)
            if btagger == "btagcsv":
                isBtagJetsLoose = (goodJets.btagDeepB > btagwpl)
                isBtagJetsMedium = (goodJets.btagDeepB > btagwpm)

            isNotBtagJetsLoose = np.invert(isBtagJetsLoose)
            nbtagsl = ak.num(goodJets[isBtagJetsLoose])

            isNotBtagJetsMedium = np.invert(isBtagJetsMedium)
            nbtagsm = ak.num(goodJets[isBtagJetsMedium])


            ######### Apply SFs #########

            if not isData:

                ### Evaluate btag weights ###
                jets_light = goodJets[goodJets.hadronFlavour==0]
                jets_bc    = goodJets[goodJets.hadronFlavour>0]

                # Workaround to use UL16APV SFs for UL16 for light jets
                year_light = year
                if year == "2016": year_light = "2016APV"

                btag_sf_light = cor_tc.btag_sf_eval(jets_light, "L",year_light,"deepJet_incl","central")
                btag_sf_bc    = cor_tc.btag_sf_eval(jets_bc,    "L",year,      "deepJet_comb","central")

                btag_eff_light = cor_ec.btag_eff_eval(jets_light,"L",year)
                btag_eff_bc = cor_ec.btag_eff_eval(jets_bc,"L",year)

                wgt_light = cor_tc.get_method1a_wgt_singlewp(btag_eff_light,btag_sf_light, jets_light.btagDeepFlavB>btagwpl)
                wgt_bc    = cor_tc.get_method1a_wgt_singlewp(btag_eff_bc,   btag_sf_bc,    jets_bc.btagDeepFlavB>btagwpl)

                wgt_btag_nom = wgt_light*wgt_bc
                weights_obj_base_for_kinematic_syst.add("btagSF", wgt_btag_nom)

                # Put the btagging up and down weight variations into the weights object
                if self._do_systematics:
                    for btag_sys in ["correlated", "uncorrelated"]:
                        year_tag = f"_{year}"
                        if btag_sys == "correlated": year_tag = ""

                        btag_sf_light_up   = cor_tc.btag_sf_eval(jets_light, "L",year_light,"deepJet_incl",f"up_{btag_sys}")
                        btag_sf_light_down = cor_tc.btag_sf_eval(jets_light, "L",year_light,"deepJet_incl",f"down_{btag_sys}")
                        btag_sf_bc_up      = cor_tc.btag_sf_eval(jets_bc,    "L",year,      "deepJet_comb",f"up_{btag_sys}")
                        btag_sf_bc_down    = cor_tc.btag_sf_eval(jets_bc,    "L",year,      "deepJet_comb",f"down_{btag_sys}")

                        wgt_light_up   = cor_tc.get_method1a_wgt_singlewp(btag_eff_light,btag_sf_light_up, jets_light.btagDeepFlavB>btagwpl)
                        wgt_bc_up      = cor_tc.get_method1a_wgt_singlewp(btag_eff_bc,   btag_sf_bc_up,    jets_bc.btagDeepFlavB>btagwpl)
                        wgt_light_down = cor_tc.get_method1a_wgt_singlewp(btag_eff_light,btag_sf_light_down, jets_light.btagDeepFlavB>btagwpl)
                        wgt_bc_down    = cor_tc.get_method1a_wgt_singlewp(btag_eff_bc,   btag_sf_bc_down,    jets_bc.btagDeepFlavB>btagwpl)

                        # Note, up and down weights scaled by 1/wgt_btag_nom so that don't double count the central btag correction (i.e. don't apply it also in the case of up and down variations)
                        weights_obj_base_for_kinematic_syst.add(f"btagSFlight_{btag_sys}{year_tag}", events.nom, wgt_light_up*wgt_bc/wgt_btag_nom, wgt_light_down*wgt_bc/wgt_btag_nom)
                        weights_obj_base_for_kinematic_syst.add(f"btagSFbc_{btag_sys}{year_tag}",    events.nom, wgt_light*wgt_bc_up/wgt_btag_nom, wgt_light*wgt_bc_down/wgt_btag_nom)


            ######### Masks we need for the selection ##########

            # Pass trigger mask
            pass_trg = es_tc.trg_pass_no_overlap(events,isData,dataset,str(year),dataset_dict=es_ec.dataset_dict,exclude_dict=es_ec.exclude_dict)
            pass_trg = (pass_trg & es_ec.trg_matching(events,year))

            # b jet masks
            bmask_atleast1med_atleast2loose = ((nbtagsm>=1)&(nbtagsl>=2)) # Used for 2lss and 4l
            bmask_exactly0loose = (nbtagsl==0) # Used for 4l WWZ SR
            bmask_exactly0med = (nbtagsm==0) # Used for 3l CR and 2los Z CR
            bmask_exactly1med = (nbtagsm==1) # Used for 3l SR and 2lss CR
            bmask_exactly2med = (nbtagsm==2) # Used for CRtt
            bmask_atleast2med = (nbtagsm>=2) # Used for 3l SR
            bmask_atmost2med  = (nbtagsm< 3) # Used to make 2lss mutually exclusive from tttt enriched
            bmask_atleast3med = (nbtagsm>=3) # Used for tttt enriched
            bmask_atleast1med = (nbtagsm>=1)
            bmask_atleast1loose = (nbtagsl>=1)
            bmask_atleast2loose = (nbtagsl>=2)


            ######### WWZ event selection stuff #########

            # Get some preliminary things we'll need
            es_ec.attach_wwz_preselection_mask(events,l_wwz_t_padded[:,0:4]) # Attach preselection sf and of flags to the events
            leps_from_z_candidate_ptordered, leps_not_z_candidate_ptordered = es_ec.get_wwz_candidates(l_wwz_t_padded[:,0:4]) # Get ahold of the leptons from the Z and from the W

            w_lep0 = leps_not_z_candidate_ptordered[:,0]
            w_lep1 = leps_not_z_candidate_ptordered[:,1]
            mll_wl0_wl1 = (w_lep0 + w_lep1).mass

            # Make masks for the SF regions
            w_candidates_mll_far_from_z = ak.fill_none(abs(mll_wl0_wl1 - get_ec_param("zmass")) > 10.0,False) # Will enforce this for SF in the PackedSelection
            ptl4 = (l0+l1+l2+l3).pt
            sf_A = ak.fill_none(met.pt >= 120.0,False) # This should never be None, but just keep syntax same as other categories
            sf_B = ak.fill_none((met.pt >= 65.0) & (met.pt < 120.0) & (ptl4 >= 70.0),False)
            sf_C = ak.fill_none((met.pt >= 65.0) & (met.pt < 120.0) & (ptl4 >= 40.0) & (ptl4 < 70.0),False)

            # Make masks for the OF regions
            of_1 = ak.fill_none((mll_wl0_wl1 >= 0.0)  & (mll_wl0_wl1 < 40.0),False)
            of_2 = ak.fill_none((mll_wl0_wl1 >= 40.0) & (mll_wl0_wl1 < 60.0),False)
            of_3 = ak.fill_none((mll_wl0_wl1 >= 60.0) & (mll_wl0_wl1 < 100.0),False)
            of_4 = ak.fill_none((mll_wl0_wl1 >= 100.0),False)

            # Mask for mt2 cut
            mt2_val = es_ec.get_mt2(w_lep0,w_lep1,met)
            mt2_mask = ak.fill_none(mt2_val>25.0,False)

            ######### Get variables #########

            l0pt = l0.pt
            j0pt = ak.flatten(j0.pt) # Flatten to go from [[j0pt],[j0pt],...] -> [j0pt,j0pt,...]
            mll_01 = (l0+l1).mass
            mllll = (l0+l1+l2+l3).mass
            scalarptsum_lep = l0.pt + l1.pt + l2.pt + l3.pt
            scalarptsum_lepmet = l0.pt + l1.pt + l2.pt + l3.pt + met.pt
            scalarptsum_lepmetjet = l0.pt + l1.pt + l2.pt + l3.pt + met.pt + ak.sum(goodJets.pt,axis=-1)

            # Get lep from Z
            z_lep0 = leps_from_z_candidate_ptordered[:,0]
            z_lep1 = leps_from_z_candidate_ptordered[:,1]

            mll_zl0_zl1 = (z_lep0 + z_lep1).mass

            pt_zl0_zl1 = (z_lep0 + z_lep1).pt
            pt_wl0_wl1 = (w_lep0 + w_lep1).pt

            dr_zl0_zl1 = z_lep0.delta_r(z_lep1)
            dr_wl0_wl1 = w_lep0.delta_r(w_lep1)
            dr_wleps_zleps = (w_lep0 + w_lep1).delta_r(z_lep0 + z_lep1)

            absdphi_zl0_zl1 = abs(z_lep0.delta_phi(z_lep1))
            absdphi_wl0_wl1 = abs(w_lep0.delta_phi(w_lep1))
            absdphi_z_ww = abs((z_lep0 + z_lep1).delta_phi(w_lep0 + w_lep1 + met))
            dphi_4l_met = (z_lep0 + z_lep1 + w_lep0 + w_lep1).delta_phi(met)
            dphi_wleps_met = (w_lep0 + w_lep1).delta_phi(met)
            dphi_zleps_met = (z_lep0 + z_lep1).delta_phi(met)

            # lb pairs (i.e. always one lep, one bjet)
            bjets = goodJets[isBtagJetsLoose]
            #lb_pairs = ak.cartesian({"l":l_wwz_t,"j":bjets})
            #mlb_min = ak.min((lb_pairs["l"] + lb_pairs["j"]).mass,axis=-1)
            #mlb_max = ak.max((lb_pairs["l"] + lb_pairs["j"]).mass,axis=-1)

            # Get BDT values
            bdt_feat_lst = [ "m_ll", "dPhi_4Lep_MET", "dPhi_Zcand_MET", "dPhi_WW_MET", "dR_Wcands", "dR_Zcands", "dR_WW_Z", "MET", "MT2", "Pt4l", "STLepHad", "STLep", "leading_Zcand_pt", "subleading_Zcand_pt", "leading_Wcand_pt", "subleading_Wcand_pt"]
            bdt_var_dict = {
                "m_ll"               : ak.fill_none(mll_wl0_wl1,-9999),
                "dPhi_4Lep_MET"      : ak.fill_none(dphi_4l_met,-9999),
                "dPhi_Zcand_MET"     : ak.fill_none(dphi_zleps_met,-9999),
                "dPhi_WW_MET"        : ak.fill_none(dphi_wleps_met,-9999),
                "dR_Wcands"          : ak.fill_none(dr_wl0_wl1,-9999),
                "dR_Zcands"          : ak.fill_none(dr_zl0_zl1,-9999),
                "dR_WW_Z"            : ak.fill_none(dr_wleps_zleps,-9999),
                "MET"                : ak.fill_none(met.pt,-9999),
                "MT2"                : ak.fill_none(mt2_val,-9999),
                "Pt4l"               : ak.fill_none(ptl4,-9999),
                "STLepHad"           : ak.fill_none(scalarptsum_lepmet,-9999),
                "STLep"              : ak.fill_none(scalarptsum_lepmetjet,-9999),
                "leading_Zcand_pt"   : ak.fill_none(z_lep0.pt,-9999),
                "subleading_Zcand_pt": ak.fill_none(z_lep1.pt,-9999),
                "leading_Wcand_pt"   : ak.fill_none(w_lep0.pt,-9999),
                "subleading_Wcand_pt": ak.fill_none(w_lep1.pt,-9999),
            }

            bdt_of_wwz_raw = os_ec.xgb_eval_wrapper(bdt_feat_lst,bdt_var_dict,ewkcoffea_path("data/wwz_zh_bdt/of_WWZ.json"))
            bdt_sf_wwz_raw = os_ec.xgb_eval_wrapper(bdt_feat_lst,bdt_var_dict,ewkcoffea_path("data/wwz_zh_bdt/sf_WWZ.json"))
            bdt_of_zh_raw  = os_ec.xgb_eval_wrapper(bdt_feat_lst,bdt_var_dict,ewkcoffea_path("data/wwz_zh_bdt/of_ZH.json"))
            bdt_sf_zh_raw  = os_ec.xgb_eval_wrapper(bdt_feat_lst,bdt_var_dict,ewkcoffea_path("data/wwz_zh_bdt/sf_ZH.json"))

            # Match TMVA's scaling https://root.cern.ch/doc/v606/MethodBDT_8cxx_source.html
            bdt_of_wwz = (2.0*((1.0+math.e**(-2*bdt_of_wwz_raw))**(-1))) - 1.0
            bdt_sf_wwz = (2.0*((1.0+math.e**(-2*bdt_sf_wwz_raw))**(-1))) - 1.0
            bdt_of_zh  = (2.0*((1.0+math.e**(-2*bdt_of_zh_raw))**(-1))) - 1.0
            bdt_sf_zh  = (2.0*((1.0+math.e**(-2*bdt_sf_zh_raw))**(-1))) - 1.0

            ### BDT SRs ###
            # SF BDT SRs
            sf_wwz_sr1 = ( (bdt_sf_wwz > 0.9) & (bdt_sf_zh  >  0.8))
            sf_wwz_sr2 = ( (bdt_sf_wwz > 0.9) & (bdt_sf_zh  > -0.6) & (bdt_sf_zh  < 0.8))
            sf_zh_sr1  = ( (bdt_sf_wwz < 0.9) & (bdt_sf_wwz >  0.7) & (bdt_sf_zh  > 0.85))
            sf_zh_sr2  = ( (bdt_sf_wwz < 0.7) & (bdt_sf_wwz >  0.6) & (bdt_sf_zh  > 0.85))
            sf_any     = ( sf_wwz_sr1 | sf_wwz_sr2 | sf_zh_sr1 | sf_zh_sr2)
            sf_wwz_sr3 = ( ~sf_any & ((bdt_sf_zh > 0.5) & (bdt_sf_wwz > 0.35)))
            sf_wwz_sr4 = ( ~(sf_any | sf_wwz_sr3) & ( (bdt_sf_zh > 0.85) & (bdt_sf_wwz > -0.5)))
            sf_zh_sr3  = ( ~(sf_any | sf_wwz_sr3 | sf_wwz_sr4) & ( bdt_sf_wwz > 0.8 ) )

            # OF BDT SRs
            of_wwz_sr1 = ( (bdt_of_wwz > 0.7) & (bdt_of_zh < -0.3) )
            of_wwz_sr2 = ( (bdt_of_wwz < 0.7) & (bdt_of_wwz > 0.4) & (bdt_of_zh < -0.6) )
            of_zh_sr1  = ( (bdt_of_wwz > 0.5) & (bdt_of_zh > 0.7) )
            of_zh_sr2  = ( (bdt_of_wwz < 0.5) & (bdt_of_wwz > -0.2) & (bdt_of_zh > 0.7) )
            of_any     = ( of_wwz_sr1 | of_wwz_sr2 | of_zh_sr1 | of_zh_sr2 )
            of_wwz_sr3 = ( ~of_any & (bdt_of_wwz > 0.0) & (bdt_of_zh < (0.8*(bdt_of_wwz-1.))) )
            of_wwz_sr4 = ( (~of_any & ~of_wwz_sr3) & (bdt_of_wwz > 0.0) )
            of_zh_sr3  = ( (~of_any & ~of_wwz_sr3 & ~of_wwz_sr4) & (bdt_of_zh > 0.5) )
            of_zh_sr4  = ( (~of_any & ~of_wwz_sr3 & ~of_wwz_sr4 & ~of_zh_sr3) & (bdt_of_zh > 0.0) & (bdt_of_wwz > -0.5) )



            ######### Store boolean masks with PackedSelection ##########

            selections = PackedSelection(dtype='uint64')

            # Lumi mask (for data)
            selections.add("is_good_lumi",lumi_mask)

            zeroj = (njets==0)

            # For WWZ selection
            trigger_and_basic_cuts = (lumi_mask & pass_trg & events.is4lWWZ) if isData else (pass_trg & events.is4lWWZ)
            tabc_with_bmask0loose = trigger_and_basic_cuts & bmask_exactly0loose
            tabc_with_bmask1loose = trigger_and_basic_cuts & bmask_atleast1loose
            tabc_b0l_wwz_presel_sf = tabc_with_bmask0loose & events.wwz_presel_sf
            tabc_b0l_sf_far_z = tabc_b0l_wwz_presel_sf & w_candidates_mll_far_from_z

            tabc_b0l_wwz_presel_of = tabc_with_bmask0loose & events.wwz_presel_of
            tabc_b0l_wwz_of_mt2 = tabc_b0l_wwz_presel_of & mt2_mask

            tabc_b1l_wwz_presel_of = tabc_with_bmask1loose & events.wwz_presel_of
            tabc_b1l_wwz_presel_sf = tabc_with_bmask1loose & events.wwz_presel_sf

            tabc_b1l_sf_far_z = tabc_b1l_wwz_presel_sf & w_candidates_mll_far_from_z
            
            selections.add("sr_4l_sf_A", (tabc_b0l_sf_far_z & sf_A))
            selections.add("sr_4l_sf_B", (tabc_b0l_sf_far_z & sf_B))
            selections.add("sr_4l_sf_C", (tabc_b0l_sf_far_z & sf_C))
            selections.add("sr_4l_of_1", (tabc_b0l_wwz_of_mt2 & of_1))
            selections.add("sr_4l_of_2", (tabc_b0l_wwz_of_mt2 & of_2))
            selections.add("sr_4l_of_3", (tabc_b0l_wwz_of_mt2 & of_3))
            selections.add("sr_4l_of_4", (tabc_b0l_wwz_presel_of & of_4))

            selections.add("all_events", (ak.ones_like(events.is4lWWZ)))
            selections.add("4l_presel", (events.is4lWWZ)) # This matches the VVV looper selection (object selection and event selection)
            
            selections.add("sr_4l_sf_presel", (tabc_b0l_sf_far_z & (met.pt > 65.0)))
            selections.add("sr_4l_of_presel", (tabc_b0l_wwz_presel_of))

            # For BDT SRs
            selections.add("sr_4l_bdt_sf_wwz_sr1", (tabc_b0l_wwz_presel_of & sf_wwz_sr1))
            selections.add("sr_4l_bdt_sf_wwz_sr2", (tabc_b0l_wwz_presel_of & sf_wwz_sr2))
            selections.add("sr_4l_bdt_sf_wwz_sr3", (tabc_b0l_wwz_presel_of & sf_wwz_sr3))
            selections.add("sr_4l_bdt_sf_wwz_sr4", (tabc_b0l_wwz_presel_of & sf_wwz_sr4))
            selections.add("sr_4l_bdt_sf_zh_sr1", (tabc_b0l_wwz_presel_of & sf_zh_sr1))
            selections.add("sr_4l_bdt_sf_zh_sr2", (tabc_b0l_wwz_presel_of & sf_zh_sr2))
            selections.add("sr_4l_bdt_sf_zh_sr3", (tabc_b0l_wwz_presel_of & sf_zh_sr3))

            selections.add("sr_4l_bdt_of_wwz_sr1", (tabc_b0l_wwz_presel_of & of_wwz_sr1))
            selections.add("sr_4l_bdt_of_wwz_sr2", (tabc_b0l_wwz_presel_of & of_wwz_sr2))
            selections.add("sr_4l_bdt_of_wwz_sr3", (tabc_b0l_wwz_presel_of & of_wwz_sr3))
            selections.add("sr_4l_bdt_of_wwz_sr4", (tabc_b0l_wwz_presel_of & of_wwz_sr4))
            selections.add("sr_4l_bdt_of_zh_sr1", (tabc_b0l_wwz_presel_of & of_zh_sr1))
            selections.add("sr_4l_bdt_of_zh_sr2", (tabc_b0l_wwz_presel_of & of_zh_sr2))
            selections.add("sr_4l_bdt_of_zh_sr3", (tabc_b0l_wwz_presel_of & of_zh_sr3))
            selections.add("sr_4l_bdt_of_zh_sr4", (tabc_b0l_wwz_presel_of & of_zh_sr4))

            # CRs
            ww_ee = ((abs(w_lep0.pdgId) == 11) & (abs(w_lep1.pdgId) == 11))
            ww_mm = ((abs(w_lep0.pdgId) == 13) & (abs(w_lep1.pdgId) == 13))
            ww_em = ((abs(w_lep0.pdgId) == 11) & (abs(w_lep1.pdgId) == 13))
            ww_me = ((abs(w_lep0.pdgId) == 13) & (abs(w_lep1.pdgId) == 11))
            selections.add("cr_4l_btag_of",            (tabc_b1l_wwz_presel_of))
            selections.add("cr_4l_btag_sf",            (tabc_b1l_wwz_presel_sf))
            selections.add("cr_4l_btag_sf_offZ",       (tabc_b1l_sf_far_z))
            selections.add("cr_4l_btag_sf_offZ_met80", (tabc_b1l_sf_far_z & (met.pt > 80.0)))
            selections.add("cr_4l_sf", (tabc_b0l_wwz_presel_sf & (~w_candidates_mll_far_from_z)))

            bdt_sr_names = [
                "sr_4l_bdt_sf_wwz_sr1",
                "sr_4l_bdt_sf_wwz_sr2",
                "sr_4l_bdt_sf_wwz_sr3",
                "sr_4l_bdt_sf_wwz_sr4",
                "sr_4l_bdt_sf_zh_sr1",
                "sr_4l_bdt_sf_zh_sr2",
                "sr_4l_bdt_sf_zh_sr3",

                "sr_4l_bdt_of_wwz_sr1",
                "sr_4l_bdt_of_wwz_sr2",
                "sr_4l_bdt_of_wwz_sr3",
                "sr_4l_bdt_of_wwz_sr4",
                "sr_4l_bdt_of_zh_sr1",
                "sr_4l_bdt_of_zh_sr2",
                "sr_4l_bdt_of_zh_sr3",
                "sr_4l_bdt_of_zh_sr4",
            ]

            cat_dict = {
                "lep_chan_lst" : [
                    "sr_4l_sf_A","sr_4l_sf_B","sr_4l_sf_C","sr_4l_of_1","sr_4l_of_2","sr_4l_of_3","sr_4l_of_4",
                    "sr_4l_sf_presel", "sr_4l_of_presel",
                    "all_events","4l_presel",
                    "cr_4l_btag_of","cr_4l_sf", "cr_4l_btag_sf", "cr_4l_btag_sf_offZ", "cr_4l_btag_sf_offZ_met80",
                ] + bdt_sr_names
            }

            ######### Fill histos #########

            dense_variables_dict = {
                "mt2" : mt2_val,
                "met" : met.pt,
                "metphi" : met.phi,
                "ptl4" : ptl4,
                "scalarptsum_lep" : scalarptsum_lep,
                "scalarptsum_lepmet" : scalarptsum_lepmet,
                "scalarptsum_lepmetjet" : scalarptsum_lepmetjet,
                "mll_01" : mll_01,
                "mllll" : mllll,
                "l0pt" : l0pt,
                "j0pt" : j0pt,

                "z_lep0_pt" : z_lep0.pt,
                "z_lep1_pt" : z_lep1.pt,
                "w_lep0_pt" : w_lep0.pt,
                "w_lep1_pt" : w_lep1.pt,
                "z_lep0_eta" : z_lep0.eta,
                "z_lep1_eta" : z_lep1.eta,
                "w_lep0_eta" : w_lep0.eta,
                "w_lep1_eta" : w_lep1.eta,
                "z_lep0_phi" : z_lep0.phi,
                "z_lep1_phi" : z_lep1.phi,
                "w_lep0_phi" : w_lep0.phi,
                "w_lep1_phi" : w_lep1.phi,

                "mll_wl0_wl1" : mll_wl0_wl1,
                "mll_zl0_zl1" : mll_zl0_zl1,

                "pt_zl0_zl1" : pt_zl0_zl1,
                "pt_wl0_wl1" : pt_wl0_wl1,
                "dr_zl0_zl1" : dr_zl0_zl1,
                "dr_wl0_wl1" : dr_wl0_wl1,
                "dr_wleps_zleps" : dr_wleps_zleps,
                "absdphi_zl0_zl1" : absdphi_zl0_zl1,
                "absdphi_wl0_wl1" : absdphi_wl0_wl1,
                "absdphi_z_ww" : absdphi_z_ww,
                "dphi_4l_met" : dphi_4l_met,
                "dphi_zleps_met" : dphi_zleps_met,
                "dphi_wleps_met" : dphi_wleps_met,

                "nleps" : nleps,
                "njets" : njets,
                "nbtagsl" : nbtagsl,

                "nleps_counts" : nleps,
                "njets_counts" : njets,
                "nbtagsl_counts" : nbtagsl,

                "absdphi_min_afas" : absdphi_min_afas,
                "absdphi_min_afos" : absdphi_min_afos,
                "absdphi_min_sfos" : absdphi_min_sfos,
                "mll_min_afas" : mll_min_afas,
                "mll_min_afos" : mll_min_afos,
                "mll_min_sfos" : mll_min_sfos,

                #"mlb_min" : mlb_min,
                #"mlb_max" : mlb_max,

                "bdt_of_wwz_raw": bdt_of_wwz_raw,
                "bdt_sf_wwz_raw": bdt_sf_wwz_raw,
                "bdt_of_zh_raw" : bdt_of_zh_raw,
                "bdt_sf_zh_raw" : bdt_sf_zh_raw,
                "bdt_of_wwz": bdt_of_wwz,
                "bdt_sf_wwz": bdt_sf_wwz,
                "bdt_of_zh" : bdt_of_zh,
                "bdt_sf_zh" : bdt_sf_zh,
            }

            dense_variables_array = ak.zip(
                dense_variables_dict,
                depth_limit=1,
            )

            dense_variables_array_with_cuts = {}
            dense_variables_array_with_cuts["lep_chan_lst"] = {}
            cached_masks = {"lep_chan_lst": {}}
            for sr_cat in cat_dict["lep_chan_lst"]:
                mask = selections.all(sr_cat)
                cached_masks["lep_chan_lst"][sr_cat] = mask
                dense_variables_array_with_cuts["lep_chan_lst"][sr_cat] = dense_variables_array[mask]
                
            
            # List the hists that are only defined for some categories
            analysis_cats = ["sr_4l_sf_A","sr_4l_sf_B","sr_4l_sf_C","sr_4l_of_1","sr_4l_of_2","sr_4l_of_3","sr_4l_of_4"] + bdt_sr_names
            exclude_var_dict = {
                "mt2" : ["all_events"],
                "ptl4" : ["all_events"],
                "j0pt" : ["all_events", "4l_presel", "sr_4l_sf_presel", "sr_4l_of_presel", "cr_4l_sf"] + analysis_cats,
                "l0pt" : ["all_events"],
                "mll_01" : ["all_events"],
                "mllll" : ["all_events"],
                "scalarptsum_lep" : ["all_events"],
                "scalarptsum_lepmet" : ["all_events"],
                "scalarptsum_lepmetjet" : ["all_events"],
                "w_lep0_pt"  : ["all_events"],
                "w_lep1_pt"  : ["all_events"],
                "z_lep0_pt"  : ["all_events"],
                "z_lep1_pt"  : ["all_events"],
                "w_lep0_eta" : ["all_events"],
                "w_lep1_eta" : ["all_events"],
                "z_lep0_eta" : ["all_events"],
                "z_lep1_eta" : ["all_events"],
                "w_lep0_phi" : ["all_events"],
                "w_lep1_phi" : ["all_events"],
                "z_lep0_phi" : ["all_events"],
                "z_lep1_phi" : ["all_events"],
                "mll_wl0_wl1" : ["all_events"],
                "mll_zl0_zl1" : ["all_events"],

                "pt_zl0_zl1" : ["all_events"],
                "pt_wl0_wl1" : ["all_events"],
                "dr_zl0_zl1" : ["all_events"],
                "dr_wl0_wl1" : ["all_events"],
                "dr_wleps_zleps" : ["all_events"],
                "absdphi_zl0_zl1" : ["all_events"],
                "absdphi_wl0_wl1" : ["all_events"],
                "absdphi_z_ww" : ["all_events"],
                "dphi_4l_met" : ["all_events"],
                "dphi_zleps_met" : ["all_events"],
                "dphi_wleps_met" : ["all_events"],

                "absdphi_min_afas" : ["all_events"],
                "absdphi_min_afos" : ["all_events"],
                "absdphi_min_sfos" : ["all_events"],
                "mll_min_afas" : ["all_events"],
                "mll_min_afos" : ["all_events"],
                "mll_min_sfos" : ["all_events"],

                "mlb_min" : ["all_events","4l_presel", "sr_4l_sf_presel", "sr_4l_of_presel", "cr_4l_sf"] + analysis_cats,
                "mlb_max" : ["all_events","4l_presel", "sr_4l_sf_presel", "sr_4l_of_presel", "cr_4l_sf"] + analysis_cats,

                "bdt_of_wwz_raw": ["all_events"],
                "bdt_sf_wwz_raw": ["all_events"],
                "bdt_of_zh_raw" : ["all_events"],
                "bdt_sf_zh_raw" : ["all_events"],
                "bdt_of_wwz": ["all_events"],
                "bdt_sf_wwz": ["all_events"],
                "bdt_of_zh" : ["all_events"],
                "bdt_sf_zh" : ["all_events"],
            }

            # Set up the list of weight fluctuations to loop over
            # For now the syst do not depend on the category, so we can figure this out outside of the filling loop
            wgt_var_lst = ["nominal"]
            if self._do_systematics:
                if not isData:
                    if (obj_corr_syst_var != "nominal"):
                        # In this case, we are dealing with systs that change the kinematics of the objs (e.g. JES)
                        # So we don't want to loop over up/down weight variations here
                        wgt_var_lst = [obj_corr_syst_var]
                    else:
                        # Otherwise we want to loop over the up/down weight variations
                        wgt_var_lst = wgt_var_lst + wgt_correction_syst_lst
                


            # Loop over the hists we want to fill
            hout = {} # This is what we'll eventually return
            masks_cache = {}
            weights_cache = {} # work around for Weights object making a new task graph each time!
            masked_weights_cache = {} # ditto            
            #print("growth before histogram loop")
            #objgraph.show_growth(limit=50, shortnames=False)
            for dense_axis_name, dense_axis_vals in dense_variables_dict.items():
                if dense_axis_name not in self._hist_lst:
                    print(f"Skipping \"{dense_axis_name}\", it is not in the list of hists to include.")
                    continue
                #print("\ndense_axis_name,vals",dense_axis_name)
                #print("dense_axis_name,vals",dense_axis_vals)

                # Create the hist for this dense axis variable
                hout[dense_axis_name] = hda.Hist(
                    hist.axis.StrCategory([], growth=True, name="process", label="process"),
                    hist.axis.StrCategory([], growth=True, name="category", label="category"),
                    hist.axis.StrCategory([], growth=True, name="systematic", label="systematic"),
                    self._dense_axes_dict[dense_axis_name],
                    storage="weight", # Keeps track of sumw2
                    name="Counts",
                )

                # Loop over weight fluctuations
                for wgt_fluct in wgt_var_lst:

                    # Get the appropriate weight fluctuation
                    if (wgt_fluct == "nominal") or (wgt_fluct in obj_corr_syst_var_list):
                        # In the case of "nominal", no weight systematic variation is used
                        weight = weights_obj_base_for_kinematic_syst.weight(None)
                    else:
                        # Otherwise get the weight from the Weights object
                        if wgt_fluct not in weights_cache:
                            weights_cache[wgt_fluct] = weights_obj_base_for_kinematic_syst.weight(wgt_fluct)
                        weight = weights_cache[wgt_fluct]

                    # Loop over categories
                    for sr_cat in cat_dict["lep_chan_lst"]:

                        # Skip filling if this variable is not relevant for this selection
                        if (dense_axis_name in exclude_var_dict) and (sr_cat in exclude_var_dict[dense_axis_name]): continue

                        # If this is a counts hist, forget the weights and just fill with unit weights
                        if isData or dense_axis_name.endswith("_counts"): weight = None #events.nom
                        #else: weights = weights_obj_base_for_kinematic_syst.partial_weight(include=["norm"]) # For testing
                        #else: weights = weights_obj_base_for_kinematic_syst.weight(None) # For testing

                        # Make the cuts mask
                        cuts_lst = [sr_cat] 
                        #if isData: cuts_lst.append("is_good_lumi") # Apply golden json requirements if this is data
                        all_cuts_mask = None #selections.all(*cuts_lst)

                        masked_weights = None
                        cache_key = tuple(cuts_lst + [wgt_fluct])
                        # Do not recalculate the mask if we've already computed it
                        if cache_key not in masks_cache:
                            weights_cache[cache_key] = None if weight is None else weight[cached_masks["lep_chan_lst"][sr_cat]]
                        masked_weights = weights_cache[cache_key]
                            
                        #run = events.run[all_cuts_mask]
                        #luminosityBlock = events.luminosityBlock[all_cuts_mask]
                        #event = events.event[all_cuts_mask]
                        #w = weights[all_cuts_mask]
                        #if dense_axis_name == "njets":
                        #    print("\nSTARTPRINT")
                        #    for i,j in enumerate(w):
                        #        out_str = f"PRINTTAG {i} {dense_axis_name} {year} {sr_cat} {event[i]} {run[i]} {luminosityBlock[i]} {w[i]}"
                        #        print(out_str,file=sys.stderr,flush=True)
                        #    print("ENDPRINT\n")
                        #print("\ndense_axis_name",dense_axis_name)
                        #print("sr_cat",sr_cat)
                        #print("dense_axis_vals[all_cuts_mask]",dense_axis_vals[all_cuts_mask])
                        #print("end")

                        # Fill the histos
                        axes_fill_info_dict = {
                            dense_axis_name : dense_variables_array_with_cuts["lep_chan_lst"][sr_cat][dense_axis_name],
                            "weight"        : masked_weights,
                            "process"       : histAxisName,
                            "category"      : sr_cat,
                            "systematic"    : wgt_fluct,
                        }
                        hout[dense_axis_name].fill(**axes_fill_info_dict)
                    #print("growth after regions loop")
                    #objgraph.show_growth(limit=50, shortnames=False)
                    #objgraph.show_refs(
                    #    objgraph.by_type("hist.dask.hist.Hist"),
                    #    filename="dask-hist-refs.png",
                    #)
                    #objgraph.show_backrefs(
		    #   objgraph.by_type("hist.dask.hist.Hist"),
                    #    filename="dask-hist-backrefs.png",
                    #)
                        

                    #sys.exit()
        print(sum(len(hout[key].dask.layers) for key in hout))
        import dask
        print(len(dask.optimize(hout)[0]["met"].dask.layers))
        return hout

    def postprocess(self, accumulator):
        return accumulator

