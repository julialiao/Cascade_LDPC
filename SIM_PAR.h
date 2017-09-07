//////////////// base matrix related
// CASCADE-LDPC parameters
#define CASCADE_NUM_LOCAL 1
#define CASCADE_LOCAL_N   {47}
#define CASCADE_GLOBAL_N  49
#define CASCADE_LOCAL_D {4}
#define CASCADE_GLOBAL_D 2
#define MAX_DECODE_ITERATION_LOCAL1	10
#define MAX_DECODE_ITERATION_LOCAL2	3
#define MAX_DECODE_ITERATION_GLOBAL	2
//#define SQUARE_MASKING


//#define N_MASK_ROW_LOCAL {0}
//#define N_MASK_COL_LOCAL {0}
//#define N_MASK_ROW_GLOBAL 0
//#define N_MASK_COL_GLOBAL 0


// RS-GC-LDPC parameters

#define BASE_MATRIX_P_LOCAL  1 //base matrix lement is (beta^p)^i
#define BASE_MTRIX_L	151// 73 
#define BASE_MATRIX_P_GLOBAL	5
#define BASE_MATRIX_PXL 151 //219 // (QC_BASE_MTRIX_L*QC_BASE_MATRIX_P)
#define CPM_SIZE	 151//219// (QC_BASE_MTRIX_L*QC_BASE_MATRIX_P)	


//#define BASE_CSHIFT	1

#define BASE_MATRIX_GLOBAL_DOWNSAMPLE	2

#define N_COL_BLK	CASCADE_GLOBAL_N
//#define N_MASK_ROW	108// 64 //10
//#define N_MASK_COL	50
//#define RANDOM_ROW_MASK
//#define RANDOM_COL_MASK
//#define EVEN_LOCAL_ROW_MASK	
//#define EVEN_LOCAL_COL_MASK
//#define ROW_MASK_BEGIN 0
//#define ROW_MASK_END  CPM_SIZE*2-1



////////// decoding related parmaters
#define MAX_DECODE_ITERATION 2
#define ROW_SHUFFLE
#define INTERLACED_MINSUM
//#define PRINT_FAIL_CHECK_NUM
#define DYNAMIC_TWO_PHASE_TH 0

//#define HARD_DECODE
//#define TWO_BIT_SOFT
//#define FIX_POINT_EN
//#define Q_BIT 2
//#define SOFT_Q_BIT


#define CNU_MAX_LLR_ABS 127
#define MAX_LLR_ABS	CNU_MAX_LLR_ABS*16*MAX_DECODE_ITERATION	//assume column degree <, or = 16

#define CNU_SCALE1	0.5
#define CNU_SCALE2  0.5
#define CNU_DEG_TH	60

//#define DEBUG_DATA3
