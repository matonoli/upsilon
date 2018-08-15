/*
 * LechLabels.h
 *
 *  Created on: 13-12-2013
 *      Author: Khaless
 */

#ifndef LECHLABELS_H_
#define LECHLABELS_H_

	// event / run

	const char * label_Run_Id =  "; Run Id [1]";
	const char * label_Run_Number =  "; Run Number [1]";
	const char * label_Event_Id =  "; Event Id [1]";
	const char * label_Event_Number =  "; Event Number [1]";

	const char * label_vz =  "; V_{z} [cm]";
	const char * label_vpd_vz =  "; V^{VPD}_{z} [cm]";
	const char * label_vzdiff =  "; V^{Rec-VPD}_{z} [cm]";
	const char * label_Event_nVtx =  "; N_{VTX} [1]";
	const char * label_Event_nVtxId =  "; VTX_{Id} [1]";

	const char * label_Event_ranking =  "; ranking [1]";

	const char * label_Event_refMult =  "; refMult [1]";
	const char * label_Event_grefMult =  "; grefMult [1]";
	const char * label_Event_nTrack =  "; N_{trk} [1]";
	const char * label_Event_nTrackMatch =  "; matched N_{trk} [1]";

	const char * label_dN_vz =  "; #frac{dN}{dV_{z}} [#frac{1}{cm}]";
	const char * label_dN_vpd_vz =  "; #frac{dN}{dV^{VPD}_{z}} [#frac{1}{cm}]";
	const char * label_dN_vzdiff =  "; #frac{dN}{dV^{Rec-VPD}_{z}} [#frac{1}{cm}]";
	const char * label_dN_Event_nVtx =  "; #frac{dN}{dN_{VTX}} [1]";
	const char * label_dN_Event_nVtxId =  "; #frac{dN}{dVTX_{Id}} [1]";

	const char * label_dN_Run_Id =  "; #frac{dN}{dId} [1]";
	const char * label_dN_Run_Number =  "; #frac{dN}{d#} [1]";
	const char * label_dN_Event_Id =  "; #frac{dN}{dId} [1]";
	const char * label_dN_Event_Number =  "; #frac{dN}{d#} [1]";

	const char * label_dN_Event_ranking =  "; #frac{dN}{dranking} [1]";

	const char * label_dN_Event_refMult =  "; #frac{dN}{drefMult} [1]";
	const char * label_dN_Event_grefMult =  "; #frac{dN}{dgrefMult} [1]";
	const char * label_dN_Event_nTrack =  "; #frac{dN}{dN_{trk}} [1]";


	const char * label_Event_adcMax =  "; ADC_{max} [1]";
	const char * label_Event_dsmMaxOnl =  "; online dsmADC_{max} [1]";
	const char * label_Event_dsmMaxOff =  "; offline dsmADC_{max} [1]";
	const char * label_Event_L0matched_p =  "; p^{L0}_{match} [#frac{GeV}{c}]";

	const char * label_dN_Event_adcMax =  "; #frac{dN}{dADC_{max}} [1]";
	const char * label_dN_Event_dsmMax =  "; #frac{dN}{ddsmADC_{max}} [1]";
	const char * label_dN_Event_L0matched_p =  "; #frac{dN^{L0}_{match}}{dADC_{max}} [#frac{c}{GeV}]";


	const char * label_Event_zdcRate =  "; ZDC coincidence rate [Hz]";
	const char * label_Event_bbcRate =  "; BBC coincidence rate [Hz]";

	const char * label_dN_Event_zdcRate =  "; #frac{dN}{dZDC} [#frac{1}{Hz}]";
	const char * label_dN_Event_bbcRate =  "; #frac{dN}{dBBC} [#frac{1}{Hz}]";

	// Upsilon analysis

	const char * label_Event_nCand =  "; candidates N_{cand} [1]";
	const char * label_Event_nVtxCand =  "; vertices with candidates N_{vtx} [1]";
	const char * label_Event_nL0 =  "; L0 towers N_{L0} [1]";

	const char * label_dN_Event_nCand =  "; #frac{dN}{dN_{cand}} [1]";
	const char * label_dN_Event_nL0 =  "; #frac{dN}{dN_{L0}} [1]";

	// particle / track

	const char * label_id =  "; id [1]";
	const char * label_q =  "; q [e]";
	const char * label_p =  "; p [#frac{GeV}{c}]";
	const char * label_pt =  "; p_{T} [#frac{GeV}{c}]";
	const char * label_pz =  "; p_{z} [#frac{GeV}{c}]";
	const char * label_eta =  "; #eta [1]";
	const char * label_phi =  "; #phi [rad]";
	const char * label_y =  "; y [1]";
	const char * label_mass =  "; m_{ee} [#frac{GeV}{c^{2}}]";
	const char * label_vr =  "; r [#frac{mm}{c}]";
	const char * label_dEdx =  "; #frac{dE}{dx} [#frac{keV}{cm}]";
	const char * label_nSigmaE =  "; n#sigmae [1]";
	const char * label_nSigmaPi =  "; n#sigma#pi [1]";

	const char * label_dN_id =  "; #frac{dN}{did} [1]";
	const char * label_dN_q =  "; #frac{dN}{dq} [#frac{1}{e}]";
	const char * label_dN_p =  "; #frac{dN}{dp} [#frac{c}{GeV}]";
	const char * label_dN_pt =  "; #frac{dN}{dp_{T}} [#frac{c}{GeV}]";
	const char * label_dN_pz =  "; #frac{dN}{dp_{z}} [#frac{c}{GeV}]";
	const char * label_dN_eta =  "; #frac{dN}{d#eta} [1]";
	const char * label_dN_phi =  "; #frac{dN}{d#phi} [1]";
	const char * label_dN_y =  "; #frac{dN}{dy} [1]";
	const char * label_dN_mass =  "; #frac{dN}{dm_{ee}} [#frac{c^{2}}{GeV}]";
	const char * label_dN_vr =  "; #frac{dN}{dr} [#frac{c}{mm}]";
	const char * label_dN_nSigmaE =  "; #frac{dN}{dn#sigmae} [1]";
	const char * label_dN_nSigmaPi =  "; #frac{dN}{dn#sigma#pi} [1]";

	const char * label_dN_eta_phi =  "; #frac{d^{2}N}{d#etad#phi} [1]";
	const char * label_dN_y_pt =  "; #frac{d^{2}N}{dydp_{T}} [#frac{c}{GeV}]";
	const char * label_dN_mass_pt =  "; #frac{d^{2}N}{dm_{ee}dp_{T}} [#frac{c^{3}}{GeV^{2}}]";
	const char * label_dN_y_pt_mass =  "; #frac{d^{3}N}{dydp_{T}dm_{ee}} [#frac{c^{3}}{GeV^{2}}]";


	const char * label_dcaG =  "; DCA_{Glob} [cm]";
	const char * label_nFitPts =  "; N_{fitPts} [1]";
	const char * label_nHitsMax =  "; N_{hitMax} [1]";
	const char * label_nFitPtRatio =  "; #frac{N_{fitPts}}{N_{hitMax}} [1]";

	const char * label_dN_dcaG =  "; #frac{dN}{dDCA_{Glob}} [#frac{1}{cm}]";
	const char * label_dN_nFitPts =  "; #frac{dN}{dN_{fitPts}} [1]";
	const char * label_dN_nHitsMax =  "; #frac{dN}{dN_{hitMax}} [1]";
	const char * label_dN_nFitPtRatio =  "; #frac{dN}{d#frac{N_{fitPts}}{N_{hitMax}}} [1]";

	// EEMC

	const char * label_stat =  "; stat [1]";
	const char * label_energy =  "; E [GeV]";
	const char * label_pedestal =  "; Pedestal [1]";
	const char * label_adc =  "; ADC [1]";
	const char * label_gain =  "; Gain [1]";

	const char * label_E =  "; E [GeV]";
	const char * label_ET =  "; E_{T} [GeV]";
	const char * label_Etow =  "; E_{tow} [GeV]";
	const char * label_Eclu =  "; E_{clu} [GeV]";
	const char * label_Efrac =  "; #frac{E_{tow}}{E_{clu}} [1]";
	const char * label_EoP =  "; #frac{E}{p} [c]";

	const char * label_Rbemc =  "; R_{TOW} [1]";
	const char * label_Rsmd =  "; R_{SMD} [1]";

	const char * label_dN_adc =  "; #frac{dN}{dADC} [1]";

	const char * label_dN_id_stat =  "; #frac{d^{2}N}{diddstat} [1]";
	const char * label_dN_id_energy =  "; #frac{d^{2}N}{diddE} [#frac{1}{GeV}]";
	const char * label_dN_id_energy_t =  "; #frac{d^{2}N}{diddE_{T}} [#frac{1}{GeV}]";
	const char * label_dN_id_pedestal =  "; #frac{d^{2}N}{diddPED} [1]";
	const char * label_dN_id_adc =  "; #frac{d^{2}N}{diddADC} [1]";
	const char * label_dN_id_gain =  "; #frac{d^{2}N}{diddGain} [1]";

	const char * label_dN_adc_ET =  "; #frac{d^{2}N}{dADCdE_{T}} [#frac{1}{GeV}]";

	const char * label_dN_E =  "; #frac{dN}{dE} [#frac{1}{GeV}]";
	const char * label_dN_ET =  "; #frac{dN}{dE_{T}} [#frac{1}{GeV}]";
	const char * label_dN_Etow =  "; #frac{dN}{dE_{tow}} [#frac{1}{GeV}]";
	const char * label_dN_Eclu =  "; #frac{dN}{dE_{clu}} [#frac{1}{GeV}]";
	const char * label_dN_Efrac =  "; #frac{dN}{dE_{tow}/E_{clu}} [1]";
	const char * label_dN_EoP =  "; #frac{dN}{dE/p} [#frac{1}{c}]";

	const char * label_dN_Rbemc =  "; #frac{dN}{dR_{TOW}} [1]";
	const char * label_dN_Rsmd =  "; #frac{dN}{dR_{SMD}} [1]";

	// Upsilon analysis

	const char * label_ups_pt =  "; p_{T}(#Upsilon) [#frac{GeV}{c}]";

	const char * label_softId =  "; softId [1]";
	const char * label_DSMadc =  "; dsmADC [1]";
	const char * label_DSMmatch =  "; #dsm match flag N_{DSM} [1]";
	const char * label_cosTheta =  "; cos(#theta) [1]";

	const char * label_dN_ups_pt =  "; #frac{dN}{dp_{T}(#Upsilon)} [#frac{c}{GeV}]";
	const char * label_eff_ups_pt =  "; #epsilon(#Upsilon)";

	const char * label_dN_softId =  "; #frac{dN}{dsoftId} [1]";
	const char * label_dN_DSMadc =  "; #frac{dN}{ddsmADC} [1]";
	const char * label_dN_DSMmatch =  "; #frac{dN}{dN_{DSM}} [1]";
	const char * label_dN_cosTheta =  "; #frac{dN}{dcos(#theta)} [1]";

	// Parameters

	const char * label_yield_ele =  "; N_{e} [1]";
	const char * label_yield_Pi =  "; N_{#pi} [1]";
	const char * label_yield_P =  "; N_{P} [1]";

	const char * label_mean_ele =  "; #mu_{e} [1]";
	const char * label_mean_Pi =  "; #mu_{#pi} [1]";
	const char * label_mean_P =  "; #mu_{P} [1]";

	const char * label_sigma_ele =  "; #sigma_{e} [1]";
	const char * label_sigma_Pi =  "; #sigma_{#pi} [1]";
	const char * label_sigma_P =  "; #sigma_{P} [1]";

	// Efficiency

	const char * label_eff =  "; #epsilon [1]";


#endif /* LECHLABELS_H_ */
