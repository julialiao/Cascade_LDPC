#ifndef CASDADE_GC_LDPC_H
#define CASDADE_GC_LDPC_H


#include "QC_LDPC.h"

class CASCADE_GC_LDPC{
public:
CASCADE_GC_LDPC();
~CASCADE_GC_LDPC();



// parameters
int cascadeN;
vector<int> localN;
vector<int> localM;
vector<int> localK;
int globalN;
int globalK;
int globalM;

//encoding functions
void baseMatrixGen();
void generatorMatrixGen();
void cascadeEncoding(vector< vector <bool>> &data2D, vector<bool> &codeword);


//decoding
void cascadeDecoding( vector <float>&llr, vector<vector<bool> > &decoded2D);
int iterationCntLocal;
int iterationCntGlobal;


vector< QC_LDPC >  cascadeQC;
private:
int baseB_localN[CASCADE_NUM_LOCAL] = CASCADE_LOCAL_N;
int baseB_localD[CASCADE_NUM_LOCAL] = CASCADE_LOCAL_D;
#ifdef N_MASK_ROW_LOCAL
int nMaskRowLocal[CASCADE_NUM_LOCAL] = N_MASK_ROW_LOCAL;
#endif
#ifdef N_MASK_COL_LOCAL
int nMaskColLocal[CASCADE_NUM_LOCAL] = N_MASK_COL_LOCAL;
#endif


};
#endif
