#include "ClusterFinder.hh"

#include <iostream>
using namespace std;


uint16_t getPeakBinOf5(uint16_t et[5], uint16_t etSum) {
#pragma HLS PIPELINE II=6
#pragma HLS ARRAY_PARTITION variable=et complete dim=0
  uint16_t iEtSum = 
    (et[0] >> 1)                +  // 0.5xet[0]
    (et[1] >> 1) + et[1]        +  // 1.5xet[1]
    (et[2] >> 1) + (et[2] << 1) +  // 2.5xet[2]
    (et[3] << 2) - (et[3] >> 1) +  // 3.5xet[3]
    (et[4] << 2) + (et[4] >> 1) ;  // 4.5xet[4]
  uint16_t iAve = 0xBEEF;
  if(     iEtSum <= etSum) iAve = 0;
  else if(iEtSum <= (etSum << 1)) iAve = 1;
  else if(iEtSum <= (etSum + (etSum << 1))) iAve = 2;
  else if(iEtSum <= (etSum << 2)) iAve = 3;
  else iAve = 4;
  return iAve;
}



bool getClustersInTower(uint16_t crystals[NCrystalsPerEtaPhi][NCrystalsPerEtaPhi],
                        uint16_t *peakEta,
                        uint16_t *peakPhi,
                        uint16_t *towerET,
                        uint16_t *clusterET,
						uint16_t *clusteret2x5) {
#pragma HLS PIPELINE II=6
#pragma HLS ARRAY_PARTITION variable=crystals complete dim=0
  uint16_t phiStripSum[NCrystalsPerEtaPhi];
#pragma HLS ARRAY_PARTITION variable=phiStripSum complete dim=0
  for(int phi = 0; phi < NCrystalsPerEtaPhi; phi++) {
#pragma HLS UNROLL
    phiStripSum[phi] = 0;
    for(int eta = 0; eta < NCrystalsPerEtaPhi; eta++) {
#pragma HLS UNROLL
      phiStripSum[phi] += crystals[eta][phi];
    }
  }
  uint16_t etaStripSum[NCrystalsPerEtaPhi];
#pragma HLS ARRAY_PARTITION variable=etaStripSum complete dim=0
  for(int eta = 0; eta < NCrystalsPerEtaPhi; eta++) {
#pragma HLS UNROLL
    etaStripSum[eta] = 0;
    for(int phi = 0; phi < NCrystalsPerEtaPhi; phi++) {
#pragma HLS UNROLL
      etaStripSum[eta] += crystals[eta][phi];
    }
  }
  // Large cluster ET is the ET of the full tower
  *towerET = 0;
  for(int phi = 0; phi < NCrystalsPerEtaPhi; phi++) {
#pragma HLS UNROLL
    *towerET += phiStripSum[phi];
  }
  *peakEta = getPeakBinOf5(etaStripSum, *towerET);
  *peakPhi = getPeakBinOf5(phiStripSum, *towerET);
  // Small cluster ET is just the 3x5 around the peak
  *clusterET = 0;
  for(int dEta = -1; dEta <= 1; dEta++) {
#pragma HLS UNROLL
      int eta = *peakEta + dEta;
        if(eta >= 0 && eta < NCrystalsPerEtaPhi) {
        	 *clusterET += etaStripSum[eta];
      }
  }
//subCluster 2X5L
  uint16_t clusterETL ;
  clusterETL =0;
	for(int dEtaL = -1; dEtaL <= 0; dEtaL++) {
#pragma HLS UNROLL
	  int eta = *peakEta + dEtaL;
	  if(eta >= 0 && eta < NCrystalsPerEtaPhi){
	    clusterETL += etaStripSum[eta];
	  }
	}
	//subCluster 2X5R
	uint16_t clusterETR ;
	clusterETR =0;
	for(int dEtaR = 0; dEtaR <= 1; dEtaR++) {
#pragma HLS UNROLL
	  int eta = *peakEta + dEtaR;
	  if(eta >= 0 && eta < NCrystalsPerEtaPhi){
	    clusterETR += etaStripSum[eta];
	  }
	}
	//cluster2X5 is equal to max of 2X5L or 2X5R
	*clusteret2x5 = biggerLR(clusterETL, clusterETR);

	return true;
}

uint16_t biggerLR(uint16_t clusterETL, uint16_t clusterETR){


//
#pragma HLS PIPELINE II=6

  uint16_t clusterf = 0;

  if(clusterETL>clusterETR)
    clusterf = clusterETL;

  else
    clusterf = clusterETR;

  return clusterf;
}


/*
bool isolation9x15(uint16_t crystals[NCrystalsPerEtaPhi][NCrystalsPerEtaPhi],
		   uint16_t *iso9x15ET,
                   uint16_t *peakEta,
		   uint16_t *peakPhi
		  )
{

#pragma HLS PIPELINE II=6
#pragma HLS ARRAY_PARTITION variable=crystals complete dim=0
  uint16_t etaStripSum1[9];

#pragma HLS ARRAY_PARTITION variable=etaStripSum1 complete dim=0

  for(int dEta = -4; dEta <= 4; dEta++) {
#pragma HLS UNROLL
    etaStripSum1[dEta+4]=0;
    for(int dPhi = -7; dPhi <= 7; dPhi++) {
#pragma HLS UNROLL
      int phi = *peakPhi + dPhi;
      int eta = *peakEta + dEta;
      if(eta >= 0 && eta < NCaloLayer1Eta && phi >= 0 && phi < NCaloLayer1Phi) {
	//	cout<<"eta=="<<eta<<"  phi=="<<phi<<"  crystals[eta][phi]=="<<crystals[eta][phi]<<"    etaStripSum1[dEta+4==]"<<etaStripSum1[dEta+4]<<endl;
	etaStripSum1[dEta+4] += crystals[eta][phi];
      }
    }
  }
  *iso9x15ET=0;
  //  cout<<"iso9x15ET=="<<*iso9x15ET<<endl;
  for(int dEta = -4; dEta <= 4; dEta++)
  {
#pragma HLS UNROLL
    //    cout<<"iso9x15ET=="<<*iso9x15ET<<"dEta: "<< dEta<<"etaStripSum1[dEta]=="<<etaStripSum1[dEta+4]<<endl;
    *iso9x15ET += etaStripSum1[dEta+4];

  }
    cout<<"iso9x15ET=="<<*iso9x15ET<<endl;
  return true;
}
*/
bool isolation10x15(uint16_t crystals[NCrystalsPerEtaPhi][NCrystalsPerEtaPhi],
		   uint16_t *iso10x15ET,
                   uint16_t *peakEta,
		   uint16_t *peakPhi,
		   uint16_t *towerET
		  )
{

#pragma HLS PIPELINE II=6
#pragma HLS ARRAY_PARTITION variable=crystals complete dim=0
  uint16_t etaStripSum1[3];

#pragma HLS ARRAY_PARTITION variable=etaStripSum1 complete dim=0
  *iso10x15ET = 0;
  for(int dEta = -1; dEta <= 1; dEta++) {
#pragma HLS UNROLL
      int eta = *towerET + dEta;
      uint16_t etaStripSum1[eta] = 0;
      if(eta >= 0 && eta < NCaloLayer1Eta) {
	//	cout<<"eta=="<<eta<<"  phi=="<<phi<<"  crystals[eta][phi]=="<<crystals[eta][phi]<<"    etaStripSum1[dEta+4==]"<<etaStripSum1[dEta+4]<<endl;
    	  *iso10x15ET += etaStripSum1[eta];
      //}
    //}
  //subCluster 10X15L
    uint16_t iso10x15ETL ;
    iso10x15ETL =0;
  	for(int dEtaL = -1; dEtaL <= 0; dEtaL++) {
  #pragma HLS UNROLL
  	  int eta = *towerET + dEtaL;
  	  if(eta >= 0 && eta < NCrystalsPerEtaPhi){
  		iso10x15ETL += etaStripSum1[eta];
  	  }
  	}
  	//subCluster 2X5R
  	uint16_t iso10x15ETR ;
  	iso10x15ETR =0;
  	for(int dEtaR = 0; dEtaR <= 1; dEtaR++) {
  #pragma HLS UNROLL
  	  int eta = *towerET + dEtaR;
  	  if(eta >= 0 && eta < NCrystalsPerEtaPhi){
  		iso10x15ETR += etaStripSum1[eta];
  	  }
  	}


  	//cluster2X5 is equal to max of 2X5L or 2X5R
  	*iso10x15ET = iso10x15biggerLR(iso10x15ETL, iso10x15ETR);
  }
  }


    cout<<"iso10x15ET=="<<*iso10x15ET<<endl;
  return true;
}
uint16_t iso10x15biggerLR(uint16_t iso10x15ETL, uint16_t iso10x15ETR){


//
#pragma HLS PIPELINE II=6

  uint16_t iso10x15f = 0;

  if(iso10x15ETL>iso10x15ETR)
	  iso10x15f = iso10x15ETL;

  else
	  iso10x15f = iso10x15ETR;

  return iso10x15f;
}
bool mergeClusters(uint16_t ieta1, uint16_t iphi1, uint16_t itet1, uint16_t icet1,uint16_t ipet1,
                   uint16_t ieta2, uint16_t iphi2, uint16_t itet2, uint16_t icet2,uint16_t ipet2,
                   uint16_t *eta1, uint16_t *phi1, uint16_t *tet1, uint16_t *cet1,uint16_t *pet1,
                   uint16_t *eta2, uint16_t *phi2, uint16_t *tet2, uint16_t *cet2,uint16_t *pet2) {
  // Check that the clusters are neighbors in eta or phi
  if((ieta1 == ieta2) || (iphi1 == iphi2)) {
    if(icet1 > icet2) {
      // Merge 2 in to 1, and set 2 to remnant energy centered in tower
      *eta1 = ieta1;
      *phi1 = iphi1;
      *cet1 = icet1 + icet2;
      *tet1 = itet1 + icet2;
      *pet1 = ipet1 + ipet2;
      *eta2 = 2;
      *phi2 = 2;
      *cet2 = 0;
      *tet2 = itet2 - icet2;
      *pet2 = 0;
    }
    else {
      // Merge 1 in to 2, and set 1 to remnant energy centered in tower
      *eta2 = ieta2;
      *phi2 = iphi2;
      *cet2 = icet2 + icet1;
      *tet2 = itet2 + icet1;
      *pet2 = ipet2 + ipet1;
      *eta1 = 2;
      *phi1 = 2;
      *cet1 = 0;
      *tet1 = itet1 - icet1;
      *pet1  = 0;
    }
  }
  else {
    *eta1 = ieta1;
    *phi1 = iphi1;
    *cet1 = icet1;
    *tet1 = itet1;
    *pet1 = ipet1;
    *eta2 = ieta2;
    *phi2 = iphi2;
    *cet2 = icet2;
    *tet2 = itet2;
    *pet2 = ipet2;
  }
  return true;
}

bool getClustersInCard(uint16_t crystals[NCaloLayer1Eta][NCaloLayer1Phi][NCrystalsPerEtaPhi][NCrystalsPerEtaPhi],
                       uint16_t peakEta[NCaloLayer1Eta][NCaloLayer1Phi],
                       uint16_t peakPhi[NCaloLayer1Eta][NCaloLayer1Phi],
                       uint16_t towerET[NCaloLayer1Eta][NCaloLayer1Phi],
                       uint16_t clusterET[NCaloLayer1Eta][NCaloLayer1Phi],
		       uint16_t clusterET2x5[NCaloLayer1Eta][NCaloLayer1Phi],
		       uint16_t iso10x15ET[NCaloLayer1Eta][NCaloLayer1Phi]) {
#pragma HLS PIPELINE II=6
#pragma HLS ARRAY_PARTITION variable=crystals complete dim=0
#pragma HLS ARRAY_PARTITION variable=peakEta complete dim=0
#pragma HLS ARRAY_PARTITION variable=peakPhi complete dim=0
#pragma HLS ARRAY_PARTITION variable=towerET complete dim=0
#pragma HLS ARRAY_PARTITION variable=clusterET complete dim=0
#pragma HLS ARRAY_PARTITION variable=clusterET2x5 complete dim=0
#pragma HLS ARRAY_PARTITION variable=iso10x15ET complete dim=0
  uint16_t preMergePeakEta[NCaloLayer1Eta][NCaloLayer1Phi];
  uint16_t preMergePeakPhi[NCaloLayer1Eta][NCaloLayer1Phi];
  uint16_t preMergeTowerET[NCaloLayer1Eta][NCaloLayer1Phi];
  uint16_t preMergeClusterET[NCaloLayer1Eta][NCaloLayer1Phi];
  uint16_t preMergeclusterET2x5[NCaloLayer1Eta][NCaloLayer1Phi];
  uint16_t preMergeiso10x15ET[NCaloLayer1Eta][NCaloLayer1Phi];
#pragma HLS ARRAY_PARTITION variable=preMergePeakEta complete dim=0
#pragma HLS ARRAY_PARTITION variable=preMergePeakPhi complete dim=0
#pragma HLS ARRAY_PARTITION variable=preMergeTowerET complete dim=0
#pragma HLS ARRAY_PARTITION variable=preMergeClusterET complete dim=0
#pragma HLS ARRAY_PARTITION variable=preMergeclusterET2x5 complete dim=0
#pragma HLS ARRAY_PARTITION variable=preMergeiso10x15ET complete dim=0
  for(int tEta = 0; tEta < NCaloLayer1Eta; tEta++) {
#pragma HLS UNROLL
    for(int tPhi = 0; tPhi < NCaloLayer1Phi; tPhi++) {
#pragma HLS UNROLL
      preMergePeakEta[tEta][tPhi] = 999;
      preMergePeakPhi[tEta][tPhi] = 999;
      preMergeTowerET[tEta][tPhi] = 0;
      preMergeClusterET[tEta][tPhi] = 0;
      preMergeclusterET2x5[tEta][tPhi] = 0;
      preMergeiso10x15ET[tEta][tPhi] = 0;
      if(!getClustersInTower(crystals[tEta][tPhi], 
                             &preMergePeakEta[tEta][tPhi],
                             &preMergePeakPhi[tEta][tPhi],
                             &preMergeTowerET[tEta][tPhi],
                             &preMergeClusterET[tEta][tPhi],
			     &preMergeclusterET2x5[tEta][tPhi])
         ) {
        return false;
      }
    }
  }
  // Merge neighboring split-clusters here
  for(int tEta = 0; tEta < NCaloLayer1Eta; tEta++) {
#pragma HLS UNROLL
    for(int tPhi = 0; tPhi < NCaloLayer1Phi; tPhi++) {
#pragma HLS UNROLL
      peakEta[tEta][tPhi] = preMergePeakEta[tEta][tPhi];
      peakPhi[tEta][tPhi] = preMergePeakPhi[tEta][tPhi];
      towerET[tEta][tPhi] = preMergeTowerET[tEta][tPhi];
      clusterET[tEta][tPhi] = preMergeClusterET[tEta][tPhi];
      clusterET2x5[tEta][tPhi] = preMergeclusterET2x5[tEta][tPhi];
      preMergeiso10x15ET[tEta][tPhi] = 0;
      int nEta = -1;
      int nPhi = -1;
      if(preMergePeakEta[tEta][tPhi] == 0 && tEta != 0) nEta = tEta - 1;
      if(preMergePeakEta[tEta][tPhi] == NCaloLayer1Phi && tEta != 16) nEta = tEta + 1;
      if(preMergePeakPhi[tEta][tPhi] == 0 && tPhi != 0) nPhi = tPhi - 1;
      if(preMergePeakPhi[tEta][tPhi] == NCaloLayer1Phi && tPhi != 3) nPhi = tPhi + 1;
      //std::cout<<"Before merging tEta/tPhi/peakEta/peakPhi/"<<tEta<<"/"<<tPhi<<"/"<<peakEta[tEta][tPhi]<<"/"<<peakPhi[tEta][tPhi]<<endl;
      if(nEta >= 0 && nEta < NCaloLayer1Eta && nPhi >= 0 && nPhi < NCaloLayer1Phi) {
	if(!mergeClusters(preMergePeakEta[tEta][tPhi],
			  preMergePeakPhi[tEta][tPhi],
			  preMergeTowerET[tEta][tPhi],
			  preMergeClusterET[tEta][tPhi],
			  preMergeclusterET2x5[tEta][tPhi],
			  preMergePeakEta[nEta][nPhi],
			  preMergePeakPhi[nEta][nPhi],
			  preMergeTowerET[nEta][nPhi],
			  preMergeClusterET[nEta][nPhi],
			  preMergeclusterET2x5[nEta][nPhi],
			  &peakEta[tEta][tPhi],
			  &peakPhi[tEta][tPhi],
			  &towerET[tEta][tPhi],
			  &clusterET[tEta][tPhi],
			  &clusterET2x5[tEta][tPhi],
			  &peakEta[nEta][nPhi],
			  &peakPhi[nEta][nPhi],
			  &towerET[nEta][nPhi],
			  &clusterET[nEta][nPhi],
			  &clusterET2x5[nEta][nPhi]))
	  return false;

      }

      //else {
      //continue;
      //}

      //std::cout<<"After merging tEta/tPhi/peakEta/peakPhi/"<<tEta<<"/"<<tPhi<<"/"<<peakEta[tEta][tPhi]<<"/"<<peakPhi[tEta][tPhi]<<endl;
	if(!isolation10x15(crystals[tEta][tPhi], &preMergeiso10x15ET[tEta][tPhi], &peakEta[tEta][tPhi], &peakPhi[tEta][tPhi])){
	  return false;
	}
	iso10x15ET[tEta][tPhi]= preMergeiso10x15ET[tEta][tPhi];

    }
  }


  //Calculate isolation
//  for(int tEta = 0; tEta < NCaloLayer1Eta; tEta++) {
//#pragma HLS UNROLL
//    for(int tPhi = 0; tPhi < NCaloLayer1Phi; tPhi++) {
//#pragma HLS UNROLL
//
//      preMergeiso9x15ET[tEta][tPhi] = 0;
//      std::cout<<"peakEta[tEta][tPhi]"<<peakEta[tEta][tPhi]<<"peakEta[tEta][tPhi]"<<peakEta[tEta][tPhi]<<std::endl;
//      if(!isolation9x15(crystals[tEta][tPhi], &preMergeiso9x15ET[tEta][tPhi], &peakEta[tEta][tPhi], &peakPhi[tEta][tPhi])){
//	return false;
//      }
//      iso9x15ET[tEta][tPhi]= preMergeiso9x15ET[tEta][tPhi];
//
//    }
//  }
  return true;
}














