#include "CASCADE_GC_LDPC.h"

CASCADE_GC_LDPC::CASCADE_GC_LDPC(){}

CASCADE_GC_LDPC::~CASCADE_GC_LDPC(){}

void CASCADE_GC_LDPC::baseMatrixGen(){

	cascadeQC.resize(CASCADE_NUM_LOCAL+1);


	// define base matrix for each local qc ldpc
	for(int p=0; p < CASCADE_NUM_LOCAL; ++p){
		// memory allocation		
		int D1 = baseB_localD[p];
		int N1 = baseB_localN[p];

		QC_LDPC* currQC = &cascadeQC[p];
		

		currQC->baseH.resize(D1);
		
		for(int i=0; i < D1; ++i){
			currQC->baseH[i].resize(N1,-1);
			for(int j=0; j <N1; ++j){
				int idxRot = ((i + 1)*j*BASE_MATRIX_P_LOCAL)%CPM_SIZE;
				currQC->baseH[i][j]= idxRot;	
				currQC->nRowBlockH = baseB_localD[p];
				currQC->nColBlockH = baseB_localN[p]; 
			}
		}

		// masking for the local qc ldpc base matrix
#ifdef SQUARE_MASKING
		//for (int i = 0; i < D1; i += 2) {
			for (int j = 0; j < N1; j+=2) {
				currQC->baseH[0][j] = -1;
			}
		//}
#endif

		currQC->cnDegree.resize(D1);
		for (int i = 0; i < D1; ++i) {
			int cnDeg = 0;
			for (int j = 0; j < N1; ++j) {
				cnDeg += (currQC->baseH[i][j] >= 0);
			}
			currQC->cnDegree[i] = cnDeg;
		}

		for(int i=0; i < currQC->baseH.size(); ++i){
			for(int j=0; j <currQC->baseH[i].size(); ++j ){
				cout << currQC->baseH[i][j] << ",";
			}
			cout << endl;
		}

	}

	

	

	// define base matrix for the globally connect qc ldpc
	QC_LDPC* currQC = &cascadeQC[CASCADE_NUM_LOCAL];
	currQC->baseH.resize(CASCADE_GLOBAL_D);
	for(int i=0; i < CASCADE_GLOBAL_D; ++i){
		currQC->baseH[i].resize(CASCADE_GLOBAL_N,-1);
		for(int j=0; j <CASCADE_GLOBAL_N; ++j){
			int idxRot = ((i + 1)*j*BASE_MATRIX_P_GLOBAL)%CPM_SIZE;
			currQC->baseH[i][j]= idxRot;	 
		}
	}
	currQC->nRowBlockH = CASCADE_GLOBAL_D;
	currQC->nColBlockH = CASCADE_GLOBAL_N; 

	// masking for the  qc ldpc base matrix
#ifdef BASE_MATRIX_GLOBAL_DOWNSAMPLE
	for (int i = 0; i < currQC->baseH.size(); i += 2) {
		for (int j = 0; j < currQC->baseH[i].size(); j+=2) {
			currQC->baseH[i][j] = -1;
		}
	}
#endif
	// update cnDegree
	currQC->cnDegree.resize(CASCADE_GLOBAL_D);
	for (int i = 0; i < CASCADE_GLOBAL_D; ++i) {
		int cnDeg = 0;
		for (int j = 0; j < CASCADE_GLOBAL_N; ++j) {
			cnDeg += (currQC->baseH[i][j] >= 0);
		}
		currQC->cnDegree[i] = cnDeg;
	}
	
	for(int i=0; i < currQC->baseH.size(); ++i){
		for(int j=0; j <currQC->baseH[i].size(); ++j ){
			cout << currQC->baseH[i][j] << ",";
		}
		cout << endl;
	}


}

void CASCADE_GC_LDPC::generatorMatrixGen(){



	// generate RS-based QC LDPC base matrix 
	baseMatrixGen();

	// CPM, masking H, H to G
	cascadeQC.resize(CASCADE_NUM_LOCAL+1);
	localM.resize(CASCADE_NUM_LOCAL);
	localN.resize(CASCADE_NUM_LOCAL);
	localK.resize(CASCADE_NUM_LOCAL);
	for(int p=0; p < CASCADE_NUM_LOCAL; ++p){
		cascadeQC[p].M = baseB_localD[p]*CPM_SIZE;
		cascadeQC[p].N = baseB_localN[p]*CPM_SIZE;
		cascadeQC[p].rowMask.resize(cascadeQC[p].M,false);
		cascadeQC[p].colMask.resize(cascadeQC[p].N,false);
#ifdef N_MASK_ROW_LOCAL
		cascadeQC[p].maskRowH();
#endif
#ifdef N_MASK_COL_LOCAL
		cascadeQC[p].maskColH();
#endif
		cascadeQC[p].B2H();
		cascadeQC[p].generatorMatrixGen();

		localN[p] = cascadeQC[p].N;
		localM[p] = cascadeQC[p].M;
		localK[p] = cascadeQC[p].K;
	}

	cascadeQC[CASCADE_NUM_LOCAL].M = CASCADE_GLOBAL_D*CPM_SIZE;
	cascadeQC[CASCADE_NUM_LOCAL].N = CASCADE_GLOBAL_N*CPM_SIZE;
	cascadeQC[CASCADE_NUM_LOCAL].rowMask.resize(cascadeQC[CASCADE_NUM_LOCAL].M,false);
	cascadeQC[CASCADE_NUM_LOCAL].colMask.resize(cascadeQC[CASCADE_NUM_LOCAL].N,false);

#ifdef N_MASK_ROW_GLOBAL
	cascadeQC[CASCADE_NUM_LOCAL].maskRowH();
#endif

#ifdef N_MASK_COL_GLOBAL
	cascadeQC[CASCADE_NUM_LOCAL].maskColH();
#endif

	cascadeQC[CASCADE_NUM_LOCAL].B2H();
	cascadeQC[CASCADE_NUM_LOCAL].generatorMatrixGen();

	globalN = cascadeQC[CASCADE_NUM_LOCAL].N;
	globalM = cascadeQC[CASCADE_NUM_LOCAL].M;
	globalK = cascadeQC[CASCADE_NUM_LOCAL].K;

}

void CASCADE_GC_LDPC::cascadeEncoding(vector< vector <bool>> &data2D, vector<bool> &casCodeword){

	vector<bool> codeword;
	for(int p=0; p < CASCADE_NUM_LOCAL; ++p){
		vector<bool> cword;
		cascadeQC[p].encode(data2D[p],cword);
		codeword.insert(codeword.end(),cword.begin(),cword.end());
		
	}

	//pad zero
	int nPadd = cascadeQC[CASCADE_NUM_LOCAL].K - codeword.size();
	for(int n=0; n < nPadd; ++n){
		codeword.push_back(false);
	}
	cascadeQC[CASCADE_NUM_LOCAL].nPadding = nPadd;
	

	cascadeQC[CASCADE_NUM_LOCAL].encode(codeword,casCodeword);
	/* note that the local codewords are permutated after global encoding*/

}

void CASCADE_GC_LDPC::cascadeDecoding( vector<float> &llr, vector<vector<bool> > &decoded2D){


	decoded2D.resize(CASCADE_NUM_LOCAL);
	
	for(int p=0; p < CASCADE_NUM_LOCAL; ++p){
		decoded2D[p].resize(cascadeQC[p].N);
	}

	int globalK = cascadeQC[CASCADE_NUM_LOCAL].dataLength;
	int globalN = cascadeQC[CASCADE_NUM_LOCAL].N;
	
	//assign local and global LLR pointer 
	vector<float*> globalLLR(globalN);
	vector<float*> localLLR(globalK);
	
	for(int i=0; i < llr.size(); ++i){
		globalLLR[i] = &llr[i];
	}
	//vector<int>::iterator globalPermIdx = cascadeQC[CASCADE_NUM_LOCAL].columnPermIdx.begin();
	int startIdx = globalN-globalK;
	for(int j=0; j < globalK; ++j){
		int pIdx = cascadeQC[CASCADE_NUM_LOCAL].columnPermIdx[j + startIdx];
		localLLR[j] = &llr[pIdx];
	}	

	// perform local decoding at least once
	//vector<float>::iterator itLocalLLR = localLLR.begin();
	int ptr = 0;
	bool allSuccess = true;
	bool allSwitch2Phase = true;
	iterationCntLocal = 0;
	for(int p=0; p < CASCADE_NUM_LOCAL; ++p){	
		cascadeQC[p].failCheckNum.clear();					  
		cascadeQC[p].minSumDecode(MAX_DECODE_ITERATION_LOCAL1, localLLR, ptr, decoded2D[p] ,DYNAMIC_TWO_PHASE_TH);
		ptr += cascadeQC[p].N;
		allSuccess &= cascadeQC[p].decodeSuccess; 		
		allSwitch2Phase&= cascadeQC[p].swith2Phase;
		iterationCntLocal += cascadeQC[p].iterationCnt;
	}
	iterationCntGlobal = 0;
	cascadeQC[CASCADE_NUM_LOCAL].failCheckNum.clear();	
	// begin turbo decoding
	
	if(allSuccess == false){ 
		//cout << "allSuccess = " << allSuccess << endl;
		for(int itr0 = 0; itr0 < MAX_DECODE_ITERATION; ++ itr0){	
			
			vector<bool> decodedGlobal;
			cascadeQC[CASCADE_NUM_LOCAL].minSumDecode(MAX_DECODE_ITERATION_GLOBAL, globalLLR, 0, decodedGlobal ,0);
			iterationCntGlobal+=cascadeQC[CASCADE_NUM_LOCAL].iterationCnt;

			allSuccess = true;
			ptr = 0;
			for(int p=0; p < CASCADE_NUM_LOCAL; ++p){						  
				cascadeQC[p].minSumDecode(MAX_DECODE_ITERATION_LOCAL2, localLLR, ptr, decoded2D[p] ,0);
				ptr += cascadeQC[p].N;
				allSuccess &= cascadeQC[p].decodeSuccess; 		
				iterationCntLocal+= cascadeQC[p].iterationCnt;
			}

			if(allSuccess){				
				break;
			}
			
		}	
		
	}

}
