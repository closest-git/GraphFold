#pragma once
#include <memory.h>

#include "GORD_util.h"


#define MLGP_CONNECTED 0x10000
#define MLGP_CONTRACT_ENSURE 0x20000

#define MLGP_DUMP 0x10000000

// MultiLevel Graph Partition
class MLGP;
typedef MLGP *hMLGP;

class MLGP {
  // The adjacency structure of the graph is stored using the compressed storage format (CSR).
  MLGP *parent=nullptr;
  G_INT_64 *project=nullptr;
  G_INT_64 *ptr=nullptr, *adj=nullptr;

  float *gain=nullptr, *wPart=nullptr, wNode_max=0, wBorder=0;
  float *wAbsorb=nullptr, *wCluster=nullptr;  //ÿ���ڵ������ֵ
  double gain_limit=0, gain_unit_edge=0;
  double wCutEdge=0;


 protected:
  G_INT_64 PickConnectComponent(G_INT_64 seed,G_INT_64 C_no,G_INT_64 *C_id,G_INT_64 T_c,G_INT_64 flag);
  G_INT_64 Verify_Symmetry(G_INT_64 flag);
  G_INT_64 _Init_AAt(G_INT_64 nCol, G_INT_64 *ptr, G_INT_64 *ind, G_INT_64 *weight, G_INT_64 nzAAt, G_INT_64 *deg, G_INT_64 flag);
  G_INT_64 _Init_A(G_INT_64 nCol, G_INT_64 *ptr, G_INT_64 *ind, G_INT_64 *weight, G_INT_64 flag);
  G_INT_64 Init_Rank_List(G_INT_64 flag);
  G_INT_64 GetBlock(G_INT_64 B_dim, G_INT_64 **B_ptr, G_INT_64 **B_adj, G_INT_64 *set, G_INT_64, G_INT_64 flag);
  G_INT_64 GetXBlock_withcut(G_INT_64 B_dim, G_INT_64 *B_w, G_INT_64 **B_ptr, G_INT_64 **B_adj, G_INT_64 *set, G_INT_64, G_INT_64 flag);
  void Separator_ISO(G_INT_64, G_INT_64, G_INT_64 flag);
  void DumpInfo(G_INT_64 flag);
  G_INT_64 UpdateGain_V(G_INT_64 no, G_INT_64 *x, float base, G_INT_64 flag);
  G_INT_64 UpdateGain_E(G_INT_64 no, G_INT_64 *x, float base, float *obj, G_INT_64 flag);
  G_INT_64 MoveSeparator(G_INT_64 flag);

  G_INT_64 Kernel_Refine(G_INT_64 flag);

  G_INT_64 VFM_Refine(float *, G_INT_64 flag);
  G_INT_64 Flow_Refine(G_INT_64 flag);
  G_INT_64 Flow_Refine_2(G_INT_64 flag);
  G_INT_64 Apart_Refine(G_INT_64 flag);
  G_INT_64 Tabu_SelectCand(SORT_BUCKET *moves, SORT_BUCKET **move_cand, G_INT_64 *modify, G_INT_64 *tenure, float *, float eval_0,
                      G_INT_64 stp, G_INT_64 flag);
  G_INT_64 Tabu_Refine_0(float *obj, G_INT_64 flag);
  G_INT_64 VFM_Perturb(G_INT_64 flag);
  G_INT_64 UpdateBestPass(float eval, G_INT_64 *nP, G_INT_64 *P_tags, float *P_eval, G_INT_64 flag);

  bool _HEM_match(G_INT_64 cur, G_INT_64 *m_1, G_INT_64 *m_2, G_INT_64 *match, G_INT_64 flag);
  G_INT_64 Match_(G_INT_64 *match, G_INT_64 *, G_INT_64 loop_max, G_INT_64 flag);
  G_INT_64 Match_N(G_INT_64 *match, G_INT_64 *, G_INT_64 loop_max, G_INT_64 flag);
  G_INT_64 Match_4(G_INT_64 *match, G_INT_64 *, G_INT_64 loop_max, float z_thrsh, G_INT_64 flag);

  G_INT_64 Pertub_2(G_INT_64 type, G_INT_64 *temp, G_INT_64 flag);

 public:
  static G_INT_64 THRSH_COARSEST, MDO, MD_SWITCH, SEPRATOR, DOWN, DUMP, PARTI, MEASURE;
  static G_INT_64 nBisection, nDown, nMDBlock, nVFM, nTaRf, nEval;
  static G_INT_64 TABU_OUTER, TABU_INNER, CAND_TABU;
  static float T_COARSEST_RELAX, STAR, STAR_0, rMove, rEval, TABU_rMove;
  static float tPartition, tBiSection, tSub, tSplit, tX, tFlowRefine, tRefine, tEval;
  static G_INT_64 DWON_LOOP_MIN;
  static G_INT_64 DumpParam(G_INT_64 flag);
  static G_INT_64 InitParam(G_INT_64 flag);

  float *wNode=nullptr, *wEdge=nullptr;
  G_INT_64 nNode=0, nEdge=0, nMaxPart=0, type=0, sepa=0, iter=0, *moveto=nullptr, *kind=nullptr, *tag=nullptr, *stmp=nullptr;
  float wNode_p=0, wEdge_p=0;  //,T_density;
  G_INT_64 nLU=0;
  bool isGreedyCutEdge = true;
  enum {
    TABU_ONLY_CUT = 0x100,
    TABU_CAND_EVAL = 0x200,
    TABU_TENTURE_LONG = 0x400,
    INNER_DYNAMIC = 0x800,
    REFINE_APART_LOOP = 0x1000,
    VFM_NO_ABSORB = 0x2000,

  };
  enum { VFM_REFINE, TABU_REFINE, KMEANS_REFINE, UPDATE_WPART = 0x2000 };
  enum { SEPARATOR_UNKNOWN = 0x0, SEPARATOR_E, SEPARATOR_V };
  enum { NORMAL_CUT = 0x1, NORMAL_ASSOCIATE, RATIO_CUT, RATIO_ASSOCIATE, MIN_CUT };
  MLGP() { }
  MLGP(G_INT_64 nCol, G_INT_64 *ptr, G_INT_64 *ind, G_INT_64 *weight, G_INT_64 t, G_INT_64 sep, G_INT_64 flag);
  MLGP(G_INT_64 nCol, G_INT_64 *ptr, G_INT_64 *ind, G_INT_64 *weight, G_INT_64 *comp, G_INT_64 t, G_INT_64 sep, G_INT_64 flag);
  MLGP(MLGP *graph, G_INT_64 nC, G_INT_64 *project, G_INT_64 t, G_INT_64 flag);
  ~MLGP(void);

  MLGP *Parent() { return parent; }

  MLGP *ProjectDown(G_INT_64 *code, G_INT_64 flag);
  MLGP *ProjectDown_2(G_INT_64 *code, G_INT_64 loop_max, G_INT_64 flag);
  //MLGP *ProjectDown_slim(G_INT_64 *code, G_INT_64 loop_max, G_INT_64 flag);
  MLGP *Sub(G_INT_64 no, G_INT_64 *map, G_INT_64 *code, G_INT_64 flag);
  float Eval(G_INT_64 *part, G_INT_64 flag,bool  isRecalc = true, float w1 = -1, float w0 = -1);
  float Objectives(G_INT_64 *part, G_INT_64 type, G_INT_64 cand, G_INT_64 flag);
  G_INT_64 Valid(G_INT_64 flag);
  G_INT_64 isSymmetry(G_INT_64 flag);

  //	G_INT_64 Compress_N( G_INT_64 *comp,G_INT_64 alg,G_INT_64 flag );
  G_INT_64 Partition_3(G_INT_64 flag);
  G_INT_64 Partition_2(G_INT_64 *nP, G_INT_64 *P_tags, float *P_eval, G_INT_64 flag);
  G_INT_64 P3_GG(G_INT_64 flag);
  G_INT_64 P2_RG(G_INT_64 flag);
  G_INT_64 P2_GreedyGrow(G_INT_64 *nP, G_INT_64 *P_tags, float *P_eval, G_INT_64 flag);
  G_INT_64 P3_GreedyGrow(G_INT_64 *nP, G_INT_64 *P_tags, float *P_eval, G_INT_64 flag);
  G_INT_64 P3_Spectral(G_INT_64 flag);
  G_INT_64 P3_Karger(G_INT_64 flag);
  //	G_INT_64 VS4( G_INT_64 flag );
  G_INT_64 Refine_3(G_INT_64 type, G_INT_64 *temp, G_INT_64 flag);
  G_INT_64 Refine_2(G_INT_64 type, G_INT_64 flag);
  G_INT_64 CheckPartition(G_INT_64 flag);
  //	bool isCoarsest( G_INT_64 flag )	{	return nNode<=THRSH_COARSEST;	}
  MLGP *GetParent() { return parent; }
  G_INT_64 SplitComponent(hMLGP **arrGraph, G_INT_64 *map, G_INT_64 flag);

  G_INT_64 ProjectUp(G_INT_64 flag);

  G_INT_64 BiSection(G_INT_64 *nS, G_INT_64 *nP_1, G_INT_64 *nP_2, G_INT_64 alg_, G_INT_64 flag);
  G_INT_64 BiSection_2(G_INT_64 *nS, G_INT_64 *nP_1, G_INT_64 *nP_2, G_INT_64 alg_, G_INT_64 level, G_INT_64 flag);
  G_INT_64 MdOrder(G_INT_64 algorithm,G_INT_64 flag);
  G_INT_64 NestSectionOrder(G_INT_64 *pCol, G_INT_64 base, G_INT_64 alg_, G_INT_64 level, G_INT_64 flag);

  G_INT_64 GetEdgeTable(G_INT_64 **e_n1, G_INT_64 **e_n2, float **w, G_INT_64 flag);
};

#ifdef __cplusplus
extern "C" {
#endif
G_INT_64 qmd_order(G_INT_64 nRow, G_INT_64 nCol, G_INT_64 nnz, G_INT_64 *rowind, G_INT_64 *colptr, G_INT_64 *pCol);
#ifdef __cplusplus
}
#endif