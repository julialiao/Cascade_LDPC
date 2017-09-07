#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>
#include <time.h>
#include "SIM_PAR.h"
#include "CASCADE_GC_LDPC.h"
#include "COMM_BOX.h"

using namespace std;
int main() {

	// initialize random seed
	srand(time(NULL));
	CASCADE_GC_LDPC casGcLdpc;
	COMM_BOX commBox;

	vector<float> snrRange;
	vector<float> uBER;
	vector<float> cBER;
	ofstream  ofs; 
	char fname[200];

	casGcLdpc.generatorMatrixGen();

vector<int> N(CASCADE_NUM_LOCAL+1);
vector<int> K(CASCADE_NUM_LOCAL+1);
#ifdef N_MASK_COL_LOCAL
	for(int p = 0; p <CASCADE_NUM_LOCAL; ++p){
		N[p] = casGcLdpc.cascadeQC[p].N - casGcLdpc.nMaskColLocal[p];
		K[p] = casGcLdpc.cascadeQC[p].K - casGcLdpc.nMaskColLocal[p];  
	}
#else
	for(int p = 0; p <CASCADE_NUM_LOCAL; ++p){
		N[p] = casGcLdpc.cascadeQC[p].N;
		K[p] = casGcLdpc.cascadeQC[p].K;  
	}
#endif
//#ifdef N_MASK_COL_GLOBAL
//	N[CASCADE_NUM_LOCAL] = casGcLdpc.cascadeQC[CASCADE_NUM_LOCAL].N - N_MASK_COL_GLOBAL;
//	K[CASCADE_NUM_LOCAL] = casGcLdpc.cascadeQC[CASCADE_NUM_LOCAL].K - N_MASK_COL_GLOBAL;  
	
//#else
	N[CASCADE_NUM_LOCAL] = casGcLdpc.cascadeQC[CASCADE_NUM_LOCAL].N;
	K[CASCADE_NUM_LOCAL] = casGcLdpc.cascadeQC[CASCADE_NUM_LOCAL].K;  
	
//#endif

	int totalK = 0;
	for(int p=0; p < CASCADE_NUM_LOCAL; ++p){
		totalK +=  K[p];
	} 
	
	float rate;	
	int effectiveK = 0;
	for(int p=0; p <=CASCADE_NUM_LOCAL; ++p ){
		rate = (float)K[p]/(float)N[p];
		cout << "LDPC # " << p << ", (N, K) = (" <<  N[p] << "," << K[p] << "), "	<< "Rate = " << rate << endl;
		effectiveK += K[p];
	}
	float eRate;
	if (MAX_DECODE_ITERATION == 0) {
		eRate = (K[0]) / (float)N[0];
		sprintf(fname,"simResults_%d_local_%d_%d.txt",CASCADE_NUM_LOCAL, N[0] ,K[0]);
		ofs.open(fname);
	}
	else {
		eRate = (effectiveK - K[CASCADE_NUM_LOCAL]) / (float)N[CASCADE_NUM_LOCAL];
		sprintf(fname,"simResults_%d_local_%d_%d.txt",CASCADE_NUM_LOCAL,N[CASCADE_NUM_LOCAL],(effectiveK - K[CASCADE_NUM_LOCAL]) );
		ofs.open(fname);
		
	}
	cout << "effective Rate = " << eRate << endl;
	
//	ofs << "LDPC (N, K) = (" << gcLdpc.N << "," << gcLdpc.dataLength << "), "	<< "Rate = " << rate << endl;

	vector< float > iteratoinCnt;
	
	double vecSNR[10] = {4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9};
	//double vecSNR[8] = {3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3};
	//double vecSNR[10] = {4.8,5,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8};// { 3.8, 3.9, 4, 4.1,4.2,4.3};
#ifdef PRINT_FAIL_CHECK_NUM
	ofstream  ofs_cfn; 
	ofs_cfn.open("fail_check_cnt.txt");
#endif
	
	for (int snrIdx = 0; snrIdx < 10; snrIdx++ ) {
		double snrOffset = -0.3;
		double SNR = snrOffset + vecSNR[snrIdx] +10 * log10(eRate);
		cout << "Eb/N0 = " << snrOffset+vecSNR[snrIdx]  << " ,SNR = " << SNR <<endl;
		
		unsigned uNoErrorBits = 0;
		unsigned cNoErrorBits = 0;
		unsigned pktSim = 1E7;
		unsigned pkNum = 0;		unsigned nErrorPkt = 0;

		int sumItrCntLocal = 0; 
		int sumItrCntGlobal = 0;
		int breakFlg = 0;
		for (int pkt = 0; pkt < pktSim; ++pkt) {

			++pkNum;
			// generate test data
			vector<vector<bool> > data2D(CASCADE_NUM_LOCAL);
			for(int p=0; p < CASCADE_NUM_LOCAL; ++p){
				commBox.GenData(2, casGcLdpc.cascadeQC[p].dataLength, data2D[p]);
			}
			
			// encoding
			vector<bool> codeword;
			casGcLdpc.cascadeEncoding(data2D, codeword);
	 
			//BPSK 
			vector<float> modSig;
			commBox.BPSK(codeword, modSig);
			//cout << "BPSK" << modSig.size() << endl;

			//AWGN channel
			vector<float> rx;
			commBox.AWGN2(SNR, modSig, rx);
			//cout << "AWGN" << rx.size() << endl;

			// 
			vector<bool> demod;
			commBox.BPSK_DeMod(rx, demod);

			uNoErrorBits += commBox.ComputeBER(demod, codeword);

			
			vector<vector<bool> > decoded2D;


			casGcLdpc.cascadeDecoding( rx, decoded2D);

			sumItrCntLocal  += casGcLdpc.iterationCntLocal;
			sumItrCntGlobal += casGcLdpc.iterationCntGlobal;

			unsigned nErr;
			for(int p=0; p < CASCADE_NUM_LOCAL; ++p){
				nErr = commBox.ComputeBER(decoded2D[p], data2D[p], casGcLdpc.cascadeQC[p].maskBits);				
				cNoErrorBits += nErr;

			}
			nErrorPkt += (nErr >0);
#ifdef PRINT_FAIL_CHECK_NUM
			for(int p=0; p <= CASCADE_NUM_LOCAL; ++p){
				for(int j=0; j <casGcLdpc.cascadeQC[p].failCheckNum.size(); ++j){
					ofs_cfn << casGcLdpc.cascadeQC[p].failCheckNum[j] << "\t";
				}
			}
			ofs_cfn << endl;
			ofs_cfn.flush();
#endif				 

					 
			
			if (cNoErrorBits > 1000 && pkNum > 1000){ break; }
		}
		
		float avgItrCntLocal = (float)sumItrCntLocal/(float)pkNum/(float)CASCADE_NUM_LOCAL;
		float avgItrCntGlobal = (float)sumItrCntGlobal/(float)pkNum;
		iteratoinCnt.push_back(avgItrCntLocal);		
		cout << "average local decode iteration " << avgItrCntLocal << endl;
		iteratoinCnt.push_back(avgItrCntGlobal);		
		cout << "average global decode iteration " << avgItrCntGlobal << endl;


		uBER.push_back((float)uNoErrorBits / (float)pkNum / (float)N[CASCADE_NUM_LOCAL]);
		cBER.push_back((float)cNoErrorBits / (float)pkNum / (float)totalK);

		cout << SNR << "\t" << (float)uNoErrorBits / (float)pkNum / (float)N[CASCADE_NUM_LOCAL] << "\t" << (float)cNoErrorBits / (float)pkNum / (float)totalK
			<< "\t" << (float)nErrorPkt / (float)pkNum << "\t"<< pkNum << endl;

		ofs << snrOffset+vecSNR[snrIdx]  << "\t"<< SNR << "\t" << (float)uNoErrorBits / (float)pkNum / (float)N[CASCADE_NUM_LOCAL] << "\t" << (float)cNoErrorBits / (float)pkNum / totalK
			<< "\t" << (float)nErrorPkt/(float)pkNum 
			<< "\t" <<avgItrCntLocal <<   "\t"<<avgItrCntGlobal << "\t"  << pkNum << endl;

		ofs.flush();
		
		if(breakFlg) break;
	}

	ofs.close();
#ifdef PRINT_FAIL_CHECK_NUM
	ofs_cfn.close();
#endif
	return 0;
}
