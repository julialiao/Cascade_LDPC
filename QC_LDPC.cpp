#include"QC_LDPC.h"
QC_LDPC::QC_LDPC(){}
QC_LDPC::~QC_LDPC(){}
void QC_LDPC::baseMatrixMasking() {
	// randomly mask one element per column
	for (int j = 0; j < nColBlockH; ++j) {
		int idx = rand() % nRowBlockH;
		baseH[idx][j] = -1;
	}
}

void QC_LDPC::maskingH(){
	

	rowMask.resize(M,false);
	colMask.resize(N,false);


	

#ifdef N_MASK_ROW
	#ifdef RANDOM_ROW_MASK

	std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 generator(rd()); 
	//std::default_random_engine generator;
  	std::uniform_int_distribution<int> distribution(ROW_MASK_BEGIN,ROW_MASK_END);
	int nRow = 0;	
	while(nRow <N_MASK_ROW ){
		int t = distribution(generator);
		if(!rowMask[t]){
			rowMask[t] = true;
			++nRow;
		}		
		
	}
	cout << "randomly drop " <<  accumulate(rowMask.begin(),rowMask.end(),0 ) << " rows" << endl;
	
	#else
#ifdef EVEN_LOCAL_ROW_MASK  // evenly drop the first rows of each local parts

	int nDropLocals = (int)(float)N_MASK_ROW/(float)BASE_MATRIX_P;
	int iOffset = 0;
	int nDrops = 0;
	int baseHd1[BASE_MATRIX_P] = BASE_MATRIX_D1;
	for(int p=0; p <BASE_MATRIX_P-1; ++p ){
		for(int i=0; i < nDropLocals; ++i){
			
			rowMask[i+iOffset] = true;
			++nDrops;
		}
		iOffset += (baseHd1[p]*BASE_MATRIX_PxL);
	}
	int i = iOffset;
	while(nDrops < N_MASK_ROW){
		rowMask[i] = true;
		++nDrops;
		++i;
	}
	//cout << nDropLocals <<endl;
	cout << "evenly dropping first (locals)" <<  accumulate(rowMask.begin(),rowMask.end(),0 ) << " rows" << endl;		
#else
	for(int i=0; i < N_MASK_ROW; i++){
		rowMask[i] = true;
		 
	}

	cout << "dropping first " <<  accumulate(rowMask.begin(),rowMask.end(),0 ) << " rows" << endl;		

#endif
	
	#endif
#endif

#ifdef N_MASK_COL

	#ifdef RANDOM_COL_MASK
	std::random_device rd2;  //Will be used to obtain a seed for the random number engine
    std::mt19937 generator2(rd2()); 
	//std::default_random_engine generator;
  	std::uniform_int_distribution<int> distribution2(0,N-1);
	int nCol = 0;	
	while(nCol <N_MASK_COL ){
		int t = distribution2(generator2);
		if(!colMask[t]){
			colMask[t] = true;
			++nCol;
		}		
		
	}
	cout << "randomly drop " <<  accumulate(colMask.begin(),colMask.end(),0 ) << " columns" << endl;
	
	#else
#ifdef EVEN_LOCAL_COL_MASK
	int nDropLocal = (int)(float)N_MASK_COL/(float)BASE_MATRIX_P;
	int jOffset = 0;
	int nDrop = 0;
	int baseJOffset[BASE_MATRIX_P] = BASE_MATRIX_COL_OFFSET;
	for(int p=0; p <BASE_MATRIX_P-1; ++p ){
		for(int j=0; j < nDropLocal; ++j){
			
			colMask[j+jOffset] = true;
			++nDrop;
		}
		jOffset += BASE_MATRIX_LAMBDA*BASE_MATRIX_PxL;
	}
	int j =  0;
	while(nDrop < N_MASK_COL){
		colMask[N-1-j] = true;
		++nDrop;
		++j;
	}

	//cout << nDropLocal << endl;
	cout << "evenly drop  " <<  accumulate(colMask.begin(),colMask.end(),0 ) << " columns" << endl;
	
#else
	for(int i=0; i < N_MASK_COL; ++i){
		colMask[N-1-i] = true;
	}
	cout << "drop last " <<  accumulate(colMask.begin(),colMask.end(),0 ) << " columns" << endl;
#endif
	
	#endif
#endif


}

void QC_LDPC::maskRowH(int nMaskRow){
	

	rowMask.resize(M,false);
	
#ifdef RANDOM_ROW_MASK

	std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 generator(rd()); 
	//std::default_random_engine generator;
  	std::uniform_int_distribution<int> distribution(ROW_MASK_BEGIN,ROW_MASK_END);
	int nRow = 0;	
	while(nRow <nMaskRow ){
		int t = distribution(generator);
		if(!rowMask[t]){
			rowMask[t] = true;
			++nRow;
		}		
		
	}
	cout << "randomly drop " <<  accumulate(rowMask.begin(),rowMask.end(),0 ) << " rows" << endl;
	
#else

		for(int i=0; i < nMaskRow; i++){
			rowMask[M-1-i] = true;
		 
	}

	cout << "dropping last " <<  accumulate(rowMask.begin(),rowMask.end(),0 ) << " rows" << endl;		


	
#endif

}

void QC_LDPC::maskColH(int nMaskCol){
	colMask.resize(N,false);
#ifdef RANDOM_COL_MASK
	std::random_device rd2;  //Will be used to obtain a seed for the random number engine
    std::mt19937 generator2(rd2()); 
	//std::default_random_engine generator;
  	std::uniform_int_distribution<int> distribution2(0,N-1);
	int nCol = 0;	
	while(nCol <nMaskCol ){
		int t = distribution2(generator2);
		if(!colMask[t]){
			colMask[t] = true;
			++nCol;
		}		
		
	}
	cout << "randomly drop " <<  accumulate(colMask.begin(),colMask.end(),0 ) << " columns" << endl;
	
#else
	for(int i=0; i < nMaskCol; ++i){
		colMask[N-1-i] = true;
	}
	cout << "drop last " <<  accumulate(colMask.begin(),colMask.end(),0 ) << " columns" << endl;
#endif



}

void QC_LDPC::B2H(){
	H.resize(M);
	for (int i = 0; i< M; ++i) {
		H[i].resize(N, 0);
	}
	int idx0 = 0;
	
	for(int i=0; i < baseH.size(); ++i){
		int idx1 = 0;
		for(int j=0; j < baseH[0].size(); ++j){

			if(baseH[i][j] >=0){
				int idxRow = idx0;
				int idxClm = idx1 + baseH[i][j];

				for(int k=0; k < CPM_SIZE; ++k){
//cout << idxRow << "," << idxClm << endl;
					if(!colMask[idxClm] && !rowMask[idxRow]){
						H[idxRow][idxClm] = 1;
					}
					idxClm = idx1 +(idxClm+1)%CPM_SIZE;
					++idxRow;
				}
			}

			idx1 += CPM_SIZE;
		}
		idx0 += CPM_SIZE;
	}
}

void QC_LDPC::generatorMatrixGen() {

	// row / column operation for systematic encoding

	 
	vector<vector<bool>> H0;
	H0 = H;

	columnPermIdx.resize(N);
	for (int j = 0; j < N; ++j) {
		columnPermIdx[j] = j;
	}
	for (int i =0; i <M; ++i) {

		bool flg = 0;
		if (H0[i][i] == 0) {
			for (int j = i + 1; j <M; ++j) {
				if (H0[j][i] == 1) {
					vector<bool> tmp = H0[i];
					H0[i] = H0[j];
					H0[j] = tmp;
					flg = 1;
					break;
				}
			}
		}
		else {
			flg = 1;
		}

		if (flg == 0) {
			bool fnd = 0;
			for (int j = i + 1; j < N; ++j) {
				for (int k = i; k < M; ++k) {
					if (H0[k][j] == 1) {
						// swap columns
						 
						for (int ii = 0; ii < M; ++ii) {
							bool tmp = H0[ii][i];
							H0[ii][i] = H0[ii][j];
							H0[ii][j] = tmp;
						}

						// swap rows
						vector<bool> tmp = H0[i];
						H0[i] = H0[k];
						H0[k] = tmp;

						//swap permutation index
						int tmpIdxP = columnPermIdx[i];
						columnPermIdx[i] = j;
						columnPermIdx[j] = tmpIdxP;
						fnd = 1;
						break;
					}
					
				}
				if (fnd) {
					break;
				}

			}
		}
		for (int j = i - 1; j >= 0; --j) {
			if (H0[j][i] == 1) {
				for (int k = 0; k < N; ++k) {
					H0[j][k] = H0[j][k]^ H0[i][k];
				}
			}

		}
		for (int j = i + 1; j<M; ++j) {
			if (H0[j][i] == 1) {
				for (int k = 0; k < N; ++k) {
					H0[j][k] = H0[j][k] ^ H0[i][k];
				}
			}
		}

	}

	// remove redundant rows and column permutation
	vector<int> idxZ;
	nzIdxG.resize(M);

	int ss = 0;
	for (int i = 0; i < M; ++i) {
		int s1 = accumulate(H0[i].begin(), H0[i].end(), 0);
		if(s1 < 2) cout << i << ": " << s1 << endl;
		if (s1 == 0) {
			idxZ.push_back(i);
		//	cout << i << ",";
		}
		ss +=  accumulate(H0[i].begin(), H0[i].end(), 0)==0;
	}

	if (idxZ.size() > 0) {
	
		int kk = 0;
		int m = 0;
		int pLen = M - idxZ.size();
		for (int i = 0; i < M; ++i) {

			if (i != idxZ[kk]) {
				for (int j = M - idxZ.size(); j < N; ++j) {
					if (H0[i][j] == 1) {

						nzIdxG[m].push_back(j - pLen);
					}
				}

				++m;

			}
			else {
				kk += (kk < idxZ.size() - 1);
			}
		}
	}
	else {
		for (int i = 0; i < M; ++i) {
			for (int j = M ; j < N; ++j) {
				if (H0[i][j] == 1) {
				//cout << j << ",";
					nzIdxG[i].push_back(j - M);
				}
				
			}
			//cout << endl;
		}
		
	}
	



	//columnPermIdx = idxZ;

 
	dataLength = N - M + idxZ.size();
	K = dataLength;
	#ifdef N_MASK_COL
	//dataLength -= N_MASK_COL;
	#endif

	maskBits.resize(dataLength,false);
	int startIdx = N-dataLength;
	for(int i=0; i < dataLength; ++i){
		 maskBits[i] = !colMask[columnPermIdx[startIdx+i]];
		 
		//cout << i << "," << columnPermIdx[startIdx+i] << "," << maskBits[i] <<"," << colMask[startIdx+i] << endl;
	}

}


void QC_LDPC::encode(vector<bool> & data, vector<bool> & codeword) {
	if (data.size() != dataLength) {
		cout << "information length miss match, whitch should be " << dataLength << "rather than" << data.size() << "\n";
		exit(0);
	}

	codeword.resize(N,0);
	vector<bool> P(N - dataLength);
 
	for (int i = 0; i < P.size(); ++i) {
		bool tmp = 0;
	//	cout<< i << ":" << nzIdxG[i].size() << ",";
		for (int j = 0; j < nzIdxG[i].size(); ++j) {
		//	cout << data[nzIdxG[i][j]];
			tmp = tmp ^ data[nzIdxG[i][j]];
		}
		P[i] = tmp;
//cout << P[i] <<  endl ;
	}

	//reorder the parity bits and the information bits
	for (int i = 0; i < N - dataLength; ++i) {
		 
		codeword[columnPermIdx[i]] = P[i];
	}
 
	int k = 0;
	for (int i = N - dataLength; i <N ; ++i) {
 
		codeword[columnPermIdx[i]] =  data[k];
		++k;
	}

	if (checkSumH(codeword)) {
		cout << "encoding check sum fails" << endl;
	}
	 



}

bool QC_LDPC::checkSumH(vector<bool> & cword) {


	int rowIdx = 0;

	for (int i = 0; i < baseH.size(); ++i) {
		
		for (int k = 0; k < CPM_SIZE; ++k) {
			bool p = 0;
			if(!rowMask[rowIdx]){
				for (int j = 0; j < baseH[i].size(); j++) {
					if (baseH[i][j] >= 0 ) {
						int idx = j*CPM_SIZE + (baseH[i][j] + k) % CPM_SIZE;
						if(!colMask[idx]){					
							p ^= cword[idx];
						}
						//else{
							//cout << "idx = " << idx << endl;
						//}
					}					
				}

				if (p) {
					return true;
				}
				
			}
			rowIdx++;
		}
	}

	return false;
}
void QC_LDPC::initDecoder(){
	
#ifndef FIX_POINT_EN
	
	memC2V.resize(M);
	int cnIdx = 0;
	for(int i=0; i < baseH.size(); ++i){
		int degree = cnDegree[i];
		for(int k=0; k < CPM_SIZE; ++k){
			memC2V[cnIdx].resize(degree);
			for(int d=0; d < degree; ++d){
				memC2V[cnIdx][d] = 0.;
			}
			cnIdx++;
		}
	}
	

	memSumLLRV2C.resize(N);
	for(int j=0; j < N; ++j){
		memSumLLRV2C[j] = 0.;
	}
 
#else

	memC2V_i.resize(M);
	int cnIdx = 0;
	for(int i=0; i < baseH.size(); ++i){
		int degree = cnDegree[i];
		for(int k=0; k < CPM_SIZE; ++k){
			memC2V_i[cnIdx].resize(degree);
			for(int d=0; d < degree; ++d){
				memC2V_i[cnIdx][d] = 0;
			}
			cnIdx++;
		}
	}
	

	memSumLLRV2C_i.resize(N);
	for(int j=0; j < N; ++j){
		memSumLLRV2C_i[j] = 0;
	}
#endif

}



void QC_LDPC::minSumDecode(vector<float> & llr, vector<bool> &decoded){

	// initialize

	initDecoder();
	memSumLLRV2C = llr;
	vector<bool> cword(N);
	//iterative decoding

	iterationCnt = MAX_DECODE_ITERATION;
cout << 517 << endl;	
	for(int itr=0; itr < MAX_DECODE_ITERATION; ++itr){

		//local check node update	
		int cnIdx = 0;		
		for(int i=0; i < nLocalRowBlk; ++i){
			for(int k=0; k <CPM_SIZE; ++k){				
				int vOffset = 0;
								 
				vector<float *> sumLLRV2C;
				vector<bool> bitMask;
				for(int j=0; j < baseH[i].size();++j){

					if(baseH[i][j]>=0){
						int vnIdx = vOffset + (baseH[i][j] + k)%CPM_SIZE;
						 
						sumLLRV2C.push_back(&memSumLLRV2C[vnIdx]);
						bitMask.push_back(colMask[vnIdx]);
						 

					}	
					vOffset += CPM_SIZE;			
				}

				if(!rowMask[cnIdx]){
					minSumCNU(sumLLRV2C, memC2V[cnIdx],bitMask);
				}


				++cnIdx;
			}
			
		}// if i

		// global check node update
		for(int i=nLocalRowBlk; i <baseH.size(); ++i ){
			for(int k=0; k <CPM_SIZE; ++k){	
				int vOffset = 0;
				
				vector<float *> sumLLRV2C;//(BASE_MATRIX_PxL);
				vector<bool> bitMask;
				for(int j=0; j < baseH[i].size(); ++j){
					if(baseH[i][j]>=0){
						int vnIdx = vOffset + (baseH[i][j] + k)%CPM_SIZE;
						
						sumLLRV2C.push_back(&memSumLLRV2C[vnIdx]);
						bitMask.push_back(colMask[vnIdx]);
						 
						 
					}
					vOffset += CPM_SIZE;	
				}
				 
				if(!rowMask[cnIdx]){
					minSumCNU(sumLLRV2C, memC2V[cnIdx],bitMask);
				}
				++cnIdx;
			}
			
		}
		
	
		//var. node update (accumulation) and decode
		int vnIdx = 0;
		
		for(int j=0; j <nColBlockH; ++j ){
			for(int k = 0; k < CPM_SIZE; ++k){
#ifndef ROW_SHUFFLE
			if(!colMask[vnIdx])
				VNU(vnIdx,llr[vnIdx],j);	
#endif
 	
			if(!colMask[vnIdx]){
				cword[vnIdx] = (memSumLLRV2C[vnIdx]< 0);	
			}
			else{
				cword[vnIdx] = 0;
			}
 			++vnIdx;
			}
		}
	/*	
		if (itr > 10) {
			cout << 10;
			for (int i = 0; i < N; ++i) {
				if (llr[i] * memSumLLRV2C[i] < 0) {
					cout << i << "," << llr[i] << "," << memSumLLRV2C[i] << endl;
				}
			}
		}*/
		// parity check
		if( checkSumH(cword)==false){
			//cout << "early termination @ " << itr<< endl;
			iterationCnt = itr+1;
			break;
		}

	}

	

	// rearrange the data
	decoded.resize(dataLength);
	int startIdx = N-dataLength;
	for(int j=0; j < dataLength; ++j){
		decoded[j] = cword[columnPermIdx[j+startIdx]];
	}


}
void QC_LDPC::minSumCNU(vector<float * >& sumLLRV2C, vector<float> &mC2V, vector<bool> bitMsk){

	
	// abs and sign
	bool sgn = 0;
	unsigned cnuDegree = sumLLRV2C.size();
	vector <float> absLLR(cnuDegree);
	vector <bool> sgnLLR(cnuDegree);
	for(int i=0; i < cnuDegree; ++i){
		if(!bitMsk[i]){
			absLLR[i] = fabs(*sumLLRV2C[i] - mC2V[i]);
			absLLR[i] > CNU_MAX_LLR_ABS ? CNU_MAX_LLR_ABS : absLLR[i];
			sgnLLR[i] = (*sumLLRV2C[i] - mC2V[i]) < 0;
			sgn ^= sgnLLR[i];
		}
		else{
			absLLR[i] = CNU_MAX_LLR_ABS;
			sgnLLR[i] = 0;
		}
		
	}

	//assume degree > 2
	float min1 = absLLR[0];
	float min2 = absLLR[1];
	unsigned idxMin1 = 0;
	

	if(min2 < min1){
		float tmp = min1;
		min1 = min2;
		min2 = tmp;
		idxMin1 = 1;
		 
		
	}

	for(int i=2; i < cnuDegree; ++i){
		if(absLLR[i] < min2){
			min2 = absLLR[i];
			 
			if(absLLR[i]<min1){
				float tmp = min1;
				min1 = absLLR[i];
				min2 = tmp;				 
				idxMin1 = i;				 
			}
		}
	}

	

	if (cnuDegree > CNU_DEG_TH) {
		min1 *= CNU_SCALE1;
		min2 *= CNU_SCALE1;
	}
	else {
		min1 *= CNU_SCALE2;
		min2 *= CNU_SCALE2;
	}
	

	 
	for(int i=0; i < cnuDegree; ++i){
		if(!bitMsk[i]){
	#ifdef ROW_SHUFFLE 
			if(i!=idxMin1){
				float tmp = sgnLLR[i]^sgn ? -min1 : min1;
				*sumLLRV2C[i] = *sumLLRV2C[i] - mC2V[i] + tmp;
				mC2V[i] = tmp;
			}
			else{
				float tmp = sgnLLR[i]^sgn ? -min2 : min2;
				*sumLLRV2C[i] = *sumLLRV2C[i] - mC2V[i] + tmp;
				mC2V[i] = tmp;
			}
			*sumLLRV2C[i] = *sumLLRV2C[i] > CNU_MAX_LLR_ABS ? CNU_MAX_LLR_ABS : *sumLLRV2C[i];
			*sumLLRV2C[i] = *sumLLRV2C[i] < -CNU_MAX_LLR_ABS ?  -CNU_MAX_LLR_ABS : *sumLLRV2C[i];
	#else
			if(i!=idxMin1){
				mC2V[i] = sgnLLR[i]^sgn ? -min1 : min1;
			}
			else{
				mC2V[i] = sgnLLR[i]^sgn ? -min2 : min2;
			}
	#endif
		
		}
	}

}


void QC_LDPC::VNU(int vIndex, float llr, int j){  

	float sumLLRV2C = llr;
	for(int i=0; i < baseH.size(); ++i){
		if(baseH[i][j]>=0){
			int cnIdx = i*CPM_SIZE +(CPM_SIZE+vIndex-baseH[i][j])%CPM_SIZE;		
			int j2 = baseHOffsetJ[i][j];			
			sumLLRV2C+= memC2V[cnIdx][j2];
		}					 
	 
	}
	
//	sumLLRV2C = sumLLRV2C > MAX_LLR_ABS ?  MAX_LLR_ABS : sumLLRV2C;
//	sumLLRV2C = sumLLRV2C < -MAX_LLR_ABS ?  -MAX_LLR_ABS : sumLLRV2C;
	memSumLLRV2C[vIndex] = sumLLRV2C;
	

}



void QC_LDPC::minSumDecode(vector<int> & llr, vector<bool> &decoded){

	// initialize

	initDecoder();
	memSumLLRV2C_i = llr;
	vector<bool> cword(N);
	//iterative decoding

	iterationCnt = MAX_DECODE_ITERATION;
	
	for(int itr=0; itr < MAX_DECODE_ITERATION; ++itr){

		//local check node update	
		int cnIdx = 0;		
		for(int i=0; i < nLocalRowBlk; ++i){
			for(int k=0; k <CPM_SIZE; ++k){				
				int vOffset = 0;
								 
				vector<int *> sumLLRV2C;
				vector<bool> bitMask;
				for(int j=0; j < baseH[i].size();++j){

					if(baseH[i][j]>=0){
						int vnIdx = vOffset + (baseH[i][j] + k)%CPM_SIZE;
						 
						sumLLRV2C.push_back(&memSumLLRV2C_i[vnIdx]);
						bitMask.push_back(colMask[vnIdx]);
						 

					}	
					vOffset += CPM_SIZE;			
				}

				if(!rowMask[cnIdx]){
					minSumCNU(sumLLRV2C, memC2V_i[cnIdx],bitMask);
				}


				++cnIdx;
			}
			
		}// if i

		// global check node update
		for(int i=nLocalRowBlk; i <baseH.size(); ++i ){
			for(int k=0; k <CPM_SIZE; ++k){	
				int vOffset = 0;
				
				vector<int *> sumLLRV2C;//(BASE_MATRIX_PxL);
				vector<bool> bitMask;
				for(int j=0; j < baseH[i].size(); ++j){
					if(baseH[i][j]>=0){
						int vnIdx = vOffset + (baseH[i][j] + k)%CPM_SIZE;
						
						sumLLRV2C.push_back(&memSumLLRV2C_i[vnIdx]);
						bitMask.push_back(colMask[vnIdx]);
						 
						 
					}
					vOffset += CPM_SIZE;	
				}
				 
				if(!rowMask[cnIdx]){
					minSumCNU(sumLLRV2C, memC2V_i[cnIdx],bitMask);
				}
				++cnIdx;
			}
			
		}
		
	
		//var. node update (accumulation) and decode
		int vnIdx = 0;
		
		for(int j=0; j <nColBlockH; ++j ){
			for(int k = 0; k < CPM_SIZE; ++k){
#ifndef ROW_SHUFFLE
			if(!colMask[vnIdx])
				VNU(vnIdx,llr[vnIdx],j);	
#endif
 	
			if(!colMask[vnIdx]){
				cword[vnIdx] = (memSumLLRV2C_i[vnIdx]< 0);	
			}
			else{
				cword[vnIdx] = 0;
			}

 			++vnIdx;
			}
		}
	/*	
		if (itr > 10) {
			cout << 10;
			for (int i = 0; i < N; ++i) {
				if (llr[i] * memSumLLRV2C[i] < 0) {
					cout << i << "," << llr[i] << "," << memSumLLRV2C[i] << endl;
				}
			}
		}*/
		// parity check
		if( checkSumH(cword)==false){
			//cout << "early termination @ " << itr<< endl;
			iterationCnt = itr+1;
			break;
		}

	}

	

	// rearrange the data
	decoded.resize(dataLength);
	int startIdx = N-dataLength;
	for(int j=0; j < dataLength; ++j){
		decoded[j] = cword[columnPermIdx[j+startIdx]];
	}


}
void QC_LDPC::minSumCNU(vector<int * >& sumLLRV2C, vector<int> &mC2V, vector<bool> bitMsk){

	
	// abs and sign
	bool sgn = 0;
	unsigned cnuDegree = sumLLRV2C.size();
	vector <int> absLLR(cnuDegree);
	vector <bool> sgnLLR(cnuDegree);
	for(int i=0; i < cnuDegree; ++i){
		if(!bitMsk[i]){
			absLLR[i] = fabs(*sumLLRV2C[i] - mC2V[i]);
			absLLR[i] > CNU_MAX_LLR_ABS ? CNU_MAX_LLR_ABS : absLLR[i];
			sgnLLR[i] = (*sumLLRV2C[i] - mC2V[i]) < 0;
			sgn ^= sgnLLR[i];
		}
		else{
			absLLR[i] = CNU_MAX_LLR_ABS;
			sgnLLR[i] = 0;
		}
		
	}

	//assume degree > 2
	int min1 = absLLR[0];
	int min2 = absLLR[1];
	unsigned idxMin1 = 0;
	
	if(min2 < min1){
		int tmp = min1;
		min1 = min2;
		min2 = tmp;
		idxMin1 = 1;
		 
		
	}

	for(int i=2; i < cnuDegree; ++i){
		if(absLLR[i] < min2){
			min2 = absLLR[i];
			 
			if(absLLR[i]<min1){
				int tmp = min1;
				min1 = absLLR[i];
				min2 = tmp;				 
				idxMin1 = i;				 
			}
		}
	}
	// scaling factor = 0.5
	min1 >> 2;
	min2 >> 1;

/*
	if (cnuDegree > CNU_DEG_TH) {
		min1 *= CNU_SCALE1;
		min2 *= CNU_SCALE1;
	}
	else {
		min1 *= CNU_SCALE2;
		min2 *= CNU_SCALE2;
	}
*/	

	 
	for(int i=0; i < cnuDegree; ++i){
		if(!bitMsk[i]){
	#ifdef ROW_SHUFFLE 
			if(i!=idxMin1){
				int tmp = sgnLLR[i]^sgn ? -min1 : min1;
				*sumLLRV2C[i] = *sumLLRV2C[i] - mC2V[i] + tmp;
				mC2V[i] = tmp;
			}
			else{
				int tmp = sgnLLR[i]^sgn ? -min2 : min2;
				*sumLLRV2C[i] = *sumLLRV2C[i] - mC2V[i] + tmp;
				mC2V[i] = tmp;
			}
			*sumLLRV2C[i] = *sumLLRV2C[i] > CNU_MAX_LLR_ABS ? CNU_MAX_LLR_ABS : *sumLLRV2C[i];
			*sumLLRV2C[i] = *sumLLRV2C[i] < -CNU_MAX_LLR_ABS ?  -CNU_MAX_LLR_ABS : *sumLLRV2C[i];
	#else
			if(i!=idxMin1){
				mC2V[i] = sgnLLR[i]^sgn ? -min1 : min1;
			}
			else{
				mC2V[i] = sgnLLR[i]^sgn ? -min2 : min2;
			}
	#endif
		
		}
	}

}


void QC_LDPC::VNU(int vIndex, int llr, int j){  

	int sumLLRV2C = llr;
	for(int i=0; i < baseH.size(); ++i){
		if(baseH[i][j]>=0){
			int cnIdx = i*CPM_SIZE +(CPM_SIZE+vIndex-baseH[i][j])%CPM_SIZE;		
			int j2 = baseHOffsetJ[i][j];			
			sumLLRV2C+= memC2V[cnIdx][j2];
		}					 
	 
	}
	
//	sumLLRV2C = sumLLRV2C > MAX_LLR_ABS ?  MAX_LLR_ABS : sumLLRV2C;
//	sumLLRV2C = sumLLRV2C < -MAX_LLR_ABS ?  -MAX_LLR_ABS : sumLLRV2C;
	memSumLLRV2C[vIndex] = sumLLRV2C;
	

}

void QC_LDPC::minSumDecode(int maxIteration,  vector<float*>llrPtrs, int llrStartIdx,  vector<bool> & decoded, int swtichThreshold ){
		// initialize

	initDecoder();
	swith2Phase = false;
	//memSumLLRV2C = llr;
	for(int i=0; i <N; ++i){
		memSumLLRV2C[i] = *llrPtrs[i+ llrStartIdx]; //*llr[i+llrStartIdx];
	}
	
	vector<bool> cword(N);
	//iterative decoding
	for(int i=0; i < N; ++i){	
		if(!colMask[i]){
			cword[i] = (memSumLLRV2C[i]< 0);	
		}
		else{
			cword[i] = 0;
		}
		++i;
		
	}
	if( checkSumH(cword)==false){
		iterationCnt = 0;
		decodeSuccess = true;
		
	}
	else{
		iterationCnt = maxIteration;
		decodeSuccess = false;
	
		for(int itr=0; itr < maxIteration; ++itr){

			
			int cnIdx = 0;		
			for(int i=0; i < nRowBlockH; ++i){
				for(int k=0; k <CPM_SIZE; ++k){				
					int vOffset = 0;
									 
					vector<float *> sumLLRV2C;
					vector<bool> bitMask;
					for(int j=0; j < baseH[i].size();++j){

						if(baseH[i][j]>=0){
							int vnIdx = vOffset + (baseH[i][j] + k)%CPM_SIZE;
							 
							sumLLRV2C.push_back(&memSumLLRV2C[vnIdx]);
							bitMask.push_back(colMask[vnIdx]);
							 

						}	
						vOffset += CPM_SIZE;			
					}

					if(!rowMask[cnIdx]){
						minSumCNU(sumLLRV2C, memC2V[cnIdx],bitMask);
					}


					++cnIdx;
				}
			
			}// if i

		
	
			//var. node update (accumulation) and decode
			int vnIdx = 0;
		
			for(int j=0; j <nColBlockH; ++j ){
				for(int k = 0; k < CPM_SIZE; ++k){
	#ifndef ROW_SHUFFLE
				if(!colMask[vnIdx])
					VNU(vnIdx,*llr[vnIdx],j);	
	#endif
	 	
				if(!colMask[vnIdx]){
					cword[vnIdx] = (memSumLLRV2C[vnIdx]< 0);	
				}
				else{
					cword[vnIdx] = 0;
				}
	 			++vnIdx;
				}
			}
		
			// parity check
 			
			
 			if( checkSumH(cword)==false){
				//cout << "early termination @ " << itr<< endl;
				iterationCnt = itr+1;
				decodeSuccess = true;
				break;
			}

			int nFailChecks =countFailChecks(cword,0,nRowBlockH);
#ifdef PRINT_FAIL_CHECK_NUM
			 failCheckNum.push_back(nFailChecks);
#endif
			if(nFailChecks < nRowBlockH){
				swith2Phase = true;
				break;
			}		
		}
	}
	

	// rearrange the data
	decoded.resize(dataLength);
	int startIdx = N-dataLength;
	for(int j=0; j < dataLength; ++j){
		decoded[j] = cword[columnPermIdx[j+startIdx]];
	}

	// output the LLR 
	for(int i=0; i < N; ++i){
		*llrPtrs[i + llrStartIdx] = memSumLLRV2C[i];
	}

}

int QC_LDPC::countFailChecks(vector<bool> & cword, int beginRowBlk, int endRowBlk) {


	int rowIdx = beginRowBlk*CPM_SIZE;

	int failCheckCount = 0;

	for (int i = beginRowBlk; i < endRowBlk; ++i) {

		for (int k = 0; k < CPM_SIZE; ++k) {
			bool p = 0;
			if (!rowMask[rowIdx]) {
				for (int j = 0; j < baseH[i].size(); j++) {
					if (baseH[i][j] >= 0) {
						int idx = j*CPM_SIZE + (baseH[i][j] + k) % CPM_SIZE;
						if (!colMask[idx]) {
							p ^= cword[idx];
						}
						//else{
						//cout << "idx = " << idx << endl;
						//}
					}
				}

				if (p) {
					++failCheckCount;
				}

			}
			rowIdx++;
		}
	}

	return failCheckCount;
}
