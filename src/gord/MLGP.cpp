

#include "MLGP.h"

#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "Flow.h"
#include "GMDO.h"
#include "gord.h"

float MLGP::T_COARSEST_RELAX = 1.2f;
float MLGP::STAR = 0.0;
float MLGP::STAR_0 = 0.0;
int MLGP::CAND_TABU = 2;    // 0: Min-Eval cand	1: First tenture-valid		2:	First
int MLGP::MD_SWITCH = 100;  //ԽСԽ��
int MLGP::THRSH_COARSEST = 50;
int MLGP::SEPRATOR = 1;  // 1 Edge	2 Node
int MLGP::DOWN = 2;      // 2 HEM match		3 RANDOM match	 4 Fast Random match	��������
int MLGP::DWON_LOOP_MIN = 1;
int MLGP::TABU_OUTER = 5, MLGP::TABU_INNER = 1;
float MLGP::TABU_rMove = 0.1f;
int MLGP::MDO = 3;  // minimum degree order
int MLGP::PARTI = 2;
int MLGP::MEASURE = MLGP::MIN_CUT;
int MLGP::DUMP = 1;
int MLGP::nBisection = 0;
int MLGP::nVFM = 0;
int MLGP::nTaRf = 0;
int MLGP::nMDBlock = 0;
int MLGP::nDown = 0, MLGP::nEval = 0;
float MLGP::rMove = 0;
float MLGP::rEval = 0;
float MLGP::tPartition = 0, MLGP::tBiSection = 0, MLGP::tSub = 0, MLGP::tSplit = 0;
float MLGP::tX = 0, MLGP::tFlowRefine = 0, MLGP::tRefine = 0, MLGP::tEval = 0;
clock_t g_tX = 0;

/*
        Copyright 2008-present, Grusoft.
        v0.1	cys
                9/25/2009
*/
void _quick_sort_2(float *arr, void *arr_2, int low, int high, int flag) {
  int i, scanUp, scanDown, SEP, *no = (int *)arr_2, t, no_pivot;
  float pivot, a;

  //	ASSERT( flag==1 );
  if (high - low <= 0)
    return;
  else if (high - low == 1) {
    if (arr[high] < arr[low]) {
      SWAP(arr[high], arr[low], a);
      SWAP(no[high], no[low], t);
    }
    return;
  }
  //	SEP=(low+high)/2;
  SEP = -1;
  for (i = low + 1; i <= high; i++) {
    if (arr[i] < arr[i - 1]) {
      SEP = i;
      break;
    }
  }
  if (SEP == -1) {  //������
    return;
  }

  pivot = arr[SEP];
  no_pivot = no[SEP];
  SWAP(arr[SEP], arr[low], a);
  SWAP(no[SEP], no[low], t);
  scanUp = low + 1;
  scanDown = high;
  do {
    while (scanUp <= scanDown && arr[scanUp] <= pivot) scanUp++;
    while (pivot < arr[scanDown]) scanDown--;
    if (scanUp < scanDown) {
      SWAP(arr[scanUp], arr[scanDown], a);
      SWAP(no[scanUp], no[scanDown], t);
    }

  } while (scanUp < scanDown);
  arr[low] = arr[scanDown];
  no[low] = no[scanDown];
  arr[scanDown] = pivot;
  no[scanDown] = no_pivot;

  _quick_sort_2(arr, arr_2, low, scanDown - 1, flag);
  _quick_sort_2(arr, arr_2, scanDown + 1, high, flag);
}

/*
        get nz in A+A' (The column ind must be sorted )
        len [j]: length of column j of A+A', except the diagonal
        temp[nCol]

        v0.1	cys
                7/5/2004 CCS2AAT_len( )
*/
int CCS2AAT_len(int nCol, int *ptr, int *ind, int *len, int *temp) {
  int p1, p2, p, i, j, pj, pj2, k, nzdiag, nzboth, nz, nzaat, isNull = temp == NULL;

  if (isNull) {
    temp = (int *)malloc(sizeof(int) * nCol);
  }
  memset(len, 0x0, sizeof(int) * nCol);
  //	cblas_iset( nCol,0,len,1 );

  nzdiag = 0;
  nzboth = 0;
  nz = ptr[nCol];

  for (k = 0; k < nCol; k++) {
    p1 = ptr[k];
    p2 = ptr[k + 1];
    /* construct A+A' */
    for (p = p1; p < p2;) {
      /* scan the upper triangular part of A */
      j = ind[p];
      if (j < k) {
        len[j]++;
        len[k]++;
        p++;
      } else if (j == k) { /* skip the diagonal */
        p++;
        nzdiag++;
        break;
      } else /* j > k */
      {
        break;
      }

      //			ASSERT (temp [j] != EMPTY) ;
      //			ASSERT (ptr [j] <= temp [j] && temp [j] <= ptr
      //[j+1]) ;
      pj2 = ptr[j + 1];
      for (pj = temp[j]; pj < pj2;) {
        i = ind[pj];
        if (i < k) {
          /* A (i,j) is only in the lower part, not in upper.
           * add both A (i,j) and A (j,i) to the matrix A+A' */
          len[i]++;
          len[j]++;
          pj++;
        } else if (i == k) {
          /* entry A (k,j) in lower part and A (j,k) in upper */
          pj++;
          nzboth++;
          break;
        } else /* i > k */
        {
          break;
        }
      }
      temp[j] = pj;
    }
    /* temp [k] points to the entry just below the diagonal in column k */
    temp[k] = p;
  }

  /* clean up, for remaining mismatched entries */
  for (j = 0; j < nCol; j++) {
    for (pj = temp[j]; pj < ptr[j + 1]; pj++) {
      i = ind[pj];
      len[i]++;
      len[j]++;
    }
  }

  nzaat = 0;
  for (k = 0; k < nCol; k++) {
    nzaat += len[k];
  }
  // nCol���ܲ�����nzDiag!!!
  nzaat += nCol;

  if (isNull) {
    free(temp);
  }
  return (nzaat);
}

/*
        v0.1	cys
                9/27/2011
*/
int MLGP::_Init_A(int nCol, int *A_ptr, int *A_ind, int *weight, int flag) {
  int nMaxEdge = A_ptr[nCol], i, j, r, ret = GORD_OK;

  nNode = nCol;
  ptr = new int[nCol + 1];
  adj = new int[nMaxEdge];
  ptr[0] = 0;
  nEdge = 0;
  for (i = 0; i < nCol; i++) {
    for (j = A_ptr[i]; j < A_ptr[i + 1]; j++) {
      r = A_ind[j];
      if (r == i) continue;
      adj[nEdge++] = r;
    }
    ptr[i + 1] = nEdge;
  }
  ASSERT(nEdge <= nMaxEdge);

  return ret;
}

/*
        v0.1	cys
                6/8/2012
*/
int MLGP::_Init_AAt(int nCol, int *A_ptr, int *A_ind, int *weight, int nzAAt, int *deg, int flag) {
  int i, j, nnz = 0, ret = GORD_OK;

  ptr = new int[nCol + 1];
  int *len = new int[nCol];
  adj = new int[nzAAt];
  if (adj == NULL) return -1;
  for (i = 0; i < nzAAt; i++) adj[i] = -1;

  ptr[0] = 0;
  for (i = 0; i < nCol; i++) {
    len[i] = deg[i];
    ptr[i + 1] = ptr[i] + deg[i];
  }
  int p1, p2, p, pj, pj2, k, *cur = new int[nCol], *temp = new int[nCol];
  for (k = 0; k < nCol; k++) cur[k] = ptr[k];

  for (k = 0; k < nCol; k++) {
    p1 = A_ptr[k];
    p2 = A_ptr[k + 1];
    for (p = p1; p < p2;) {
      j = A_ind[p];
      if (j < k) { /* scan the upper triangular part of A */
        adj[cur[j]++] = k;
        adj[cur[k]++] = j;
        p++;
      } else if (j == k) { /* skip the diagonal */
        p++;
        break;
      } else /* j > k */
      {
        break;
      }

      pj2 = A_ptr[j + 1];
      for (pj = temp[j]; pj < pj2;) {
        i = A_ind[pj];
        if (i < k) { /*only in the lower part*/
          adj[cur[i]++] = j;
          adj[cur[j]++] = i;
          pj++;
        } else if (i == k) {
          pj++;
          break;
        } else /* i > k */
        {
          break;
        }
      }
      temp[j] = pj;
    }
    temp[k] = p;
  }
  /* clean up, for remaining mismatched entries */
  for (j = 0; j < nCol; j++) {
    for (pj = temp[j]; pj < A_ptr[j + 1]; pj++) {
      i = A_ind[pj];
      adj[cur[i]++] = j;
      adj[cur[j]++] = i;
    }
  }
  delete[] len;
  delete[] cur;
  delete[] temp;

  /*	for( i = 0; i < nCol; i++ )		temp[i]=0;
          for( i = 0; i < nCol; i++ )	{
                  for( j = A_ptr[i]; j < A_ptr[i+1]; j++ )	{
                          r = A_ind[j];
                          if( r==i )
                                  continue;
                          temp[r]++;		temp[i]++;
                          nnz+=2;
                  }
          }
          ptr=new int[nCol+1];
          ptr[nCol] = nnz;
          for( i=nCol-1; i >=0; i-- )		{
                  ptr[i] = ptr[i+1]-temp[i];		temp[i]=ptr[i];
          }
          ASSERT( temp[0]==0 );

          adj=new int[nnz];
          for( i = 0; i < nCol; i++ )	{
                  for( j = A_ptr[i]; j < A_ptr[i+1]; j++ )	{
                          r = A_ind[j];
                          if( r==i )
                                  continue;
                          adj[temp[i]++] = r;
                          adj[temp[r]++] = i;
                  }
          }
          for( i=nCol-1; i >=0; i-- )		{
                  ASSERT( temp[i]==ptr[i+1] );
          }*/

  nNode = nCol;
  nEdge = ptr[nCol];
  ASSERT(nEdge <= nzAAt);
  //	ASSERT( isSymmetry(0x0)==1 );

  return ret;
}

/*
        Init from CCS formant of MATRIX
        v0.1	cys
                9/27/2011
*/
MLGP::MLGP(int nCol, int *ccs_ptr, int *ccs_ind, int *weight, int t, int sep, int flag) {
  //memset(this, 0x0, sizeof(MLGP));    isGreedyCutEdge = true;

  int *deg = new int[nCol], nzAAT = CCS2AAT_len(nCol, ccs_ptr, ccs_ind, deg, NULL);
  if (nzAAT != ccs_ptr[nCol]) {
    _Init_AAt(nCol, ccs_ptr, ccs_ind, weight, nzAAT, deg, flag);
  } else {
    _Init_A(nCol, ccs_ptr, ccs_ind, weight, flag);
  }
  delete[] deg;
  ASSERT(isSymmetry(0x0) == 1);

  sepa = sep;
  type = t;
  nMaxPart = 2;
  wPart = new float[nMaxPart + 1];
  wNode_p = nNode;
  wEdge_p = nEdge;

  wNode_max = T_COARSEST_RELAX * (nNode) / THRSH_COARSEST;

  int i, j;
  wEdge = new float[nEdge];
  if (weight != NULL) {
    for (i = 0; i < nEdge; i++) wEdge[i] = weight[i];
  } else {
    for (i = 0; i < nEdge; i++) wEdge[i] = 1.0;
  }

  wNode = new float[nNode];
  tag = new int[nNode];
  wAbsorb = new float[nNode];
  gain = new float[nNode];
  stmp = new int[nNode];
  moveto = new int[nNode];
  kind = new int[nNode];

  for (i = 0; i < nNode; i++) {
    wNode[i] = 1.0;
    tag[i] = -1;
    wAbsorb[i] = 0.0;
    gain[i] = 0.0;
    moveto[i] = i;
    kind[i] = 0;
  }
  //	Init_Rank_List( 0x0 );
  gain_limit = 0.0;
  gain_unit_edge = FLT_MAX;
  double a = 0;
  for (i = 0; i < nNode; i++) {
    a = 0;
    for (j = ptr[i]; j < ptr[i + 1]; j++) {
      a += wEdge[j];
      gain_unit_edge = MIN(gain_unit_edge, wEdge[j]);
    }
    gain_limit = MAX(gain_limit, a);
  }
  gain_unit_edge = MAX(gain_unit_edge, 1.0);
}

/*
        Init from project of GRAPH
        v0.1	cys
                9/28/2011
*/
MLGP::MLGP(MLGP *hP, int nC, int *proj, int t, int flag) {
  //memset(this, 0x0, sizeof(MLGP));    isGreedyCutEdge = true;

  sepa = hP->sepa;
  type = t;

  nMaxPart = hP->nMaxPart;
  wPart = new float[nMaxPart + 1];
  int i, j, k, *head = NULL, *next = NULL, *adj_pos = NULL, stp = 0, no, cur;
  nNode = nC;
  nEdge = 0;
  int nMaxEdge = hP->nEdge;
  adj = new int[nMaxEdge];
  wEdge = new float[nMaxEdge];
  for (i = 0; i < nMaxEdge; i++) wEdge[i] = 0.0;

  wAbsorb = new float[nNode];
  gain = new float[nNode];
  wNode = new float[nNode];
  tag = new int[nNode];
  head = new int[nNode];
  adj_pos = new int[nNode];
  stmp = new int[nNode];
  moveto = new int[nNode];
  kind = new int[nNode];
  for (i = 0; i < nNode; i++) {
    wNode[i] = 0.0;
    tag[i] = -1;
    head[i] = -1;
    stmp[i] = 0;
    wAbsorb[i] = 0.0;
    gain[i] = 0.0;
    moveto[i] = i;
    kind[i] = 0;
  }
  next = new int[hP->nNode];
  for (i = 0; i < hP->nNode; i++) {
    next[i] = -1;
  }
  for (i = 0; i < hP->nNode; i++) {
    no = proj[i];
    if (no == -1) continue;
    ASSERT(no >= 0 && no < nC);
    if (head[no] == -1)
      head[no] = i;
    else {
      cur = head[no];
      next[i] = cur;
      head[no] = i;
    }
  }
double t_2 = G_TIC();
  wNode_p = 0.0;
  wEdge_p = -1.0;
  ptr = new int[nNode + 1];
  ptr[0] = 0;
  for (i = 0; i < nNode; i++) {
    stp++;
    cur = head[i];
    while (cur != -1) {
      ASSERT(proj[cur] == i);
      wNode_p += hP->wNode[cur];
      wNode[i] += hP->wNode[cur];
      wAbsorb[i] += hP->wAbsorb[cur];
      for (j = hP->ptr[cur]; j < hP->ptr[cur + 1]; j++) {
        no = proj[hP->adj[j]];
        if (no == -1) continue;
        if (no == i) {
          wAbsorb[i] += hP->wEdge[j];
          continue;
        }
        if (stmp[no] == stp) {  //�Ѵ���
          k = adj_pos[no];
          ASSERT(adj[k] == no);
          wEdge[k] += hP->wEdge[j];
          /*	int isFind=0;
                  for( k = ptr[i]; k < nEdge; k++ )	{
                          if( adj[k]==no )
                          {	wEdge[k] += hP->wEdge[j];	isFind=1;
             break;	}
                  }
                  ASSERT( isFind==1 );*/
        } else {
          stmp[no] = stp;
          adj_pos[no] = nEdge;
          adj[nEdge] = no;
          wEdge[nEdge] += hP->wEdge[j];
          nEdge++;
        }
      }
      cur = next[cur];
    }
    ptr[i + 1] = nEdge;
  }
  MLGP::tX += G_TIC() - t_2;

  double a = 0;
  for (i = 0; i < nNode; i++) {
    a = MAX(a, wNode[i]);
  }
  wNode_max = (float)(T_COARSEST_RELAX * MAX(a, (wNode_p) / THRSH_COARSEST));

  // wNode_max = hP->wNode_max;
  parent = hP;
  project = proj;

  delete[] head;
  delete[] next;
  delete[] adj_pos;
  //	float C_n=wNode_p*1.0/nNode,C_e = wEdge_p*1.0/nNode;
  //	T_density = hP->T_density;

  //	Init_Rank_List( 0x0 );
  gain_limit = 0.0;
  gain_unit_edge = FLT_MAX;
  a = 0;
  for (i = 0; i < nNode; i++) {
    a = 0;
    for (j = ptr[i]; j < ptr[i + 1]; j++) {
      //		if( wEdge[j]==36 )
      //			wEdge[j] = 36;
      a += wEdge[j];
      gain_unit_edge = MIN(gain_unit_edge, wEdge[j]);
    }
    gain_limit = MAX(gain_limit, a);
  }
  gain_unit_edge = MAX(gain_unit_edge, 1.0);
  
}

int MLGP::isSymmetry(int flag) {
  int isSymm = 1, i, j, k, no, isMatch;
  for (i = 0; i < nNode; i++) {
    for (j = ptr[i]; j < ptr[i + 1]; j++) {
      no = adj[j];
      isMatch = 0;
      for (k = ptr[no]; k < ptr[no + 1]; k++) {
        if (adj[k] == i) {
          isMatch = 1;
          break;
        }
      }
      ASSERT(isMatch == 1);
      if (isMatch != 1) {
        isSymm = -1;
        goto END;
      }
    }
  }
END:
  return isSymm;
}
/*
        Init from CCS formant of MATRIX
        v0.1	cys
                5/28/2012
*/
MLGP::MLGP(int nCol, int *ccs_ptr, int *ccs_ind, int *weight, int *proj, int t, int sep, int flag) {
  //memset(this, 0x0, sizeof(MLGP));    isGreedyCutEdge = true;

  sepa = sep;
  type = t;
  nMaxPart = 2;
  wPart = new float[nMaxPart + 1];

  int *next = new int[nCol], *head = new int[nCol], i, j, k, no, stp = 0, cur;
  int nMaxEdge = 0;
  for (i = 0; i < nCol; i++) {
    next[i] = -1;
    head[i] = -1;
  }
  for (i = 0; i < nCol; i++) {
    no = proj[i];
    if (no == -1) continue;
    ASSERT(no >= 0 && no < nCol);
    if (head[no] == -1) {
      head[no] = i;
      nNode++;
      nMaxEdge += ccs_ptr[i + 1] - ccs_ptr[i];
    } else {
      cur = head[no];
      next[i] = cur;
      head[no] = i;
    }
  }

  wNode = new float[nNode];
  tag = new int[nNode];
  wAbsorb = new float[nNode];
  gain = new float[nNode];
  stmp = new int[nNode];
  moveto = new int[nNode];
  kind = new int[nNode];
  for (i = 0; i < nNode; i++) {
    wNode[i] = 0.0;
    tag[i] = -1;
    wAbsorb[i] = 0.0;
    gain[i] = 0.0;
    moveto[i] = i;
    kind[i] = 0;
    stmp[i] = stp;
  }

  ptr = new int[nNode + 1];
  wEdge = new float[nMaxEdge];
  adj = new int[nMaxEdge];
  ptr[0] = 0;
  for (i = 0; i < nNode; i++) {
    stp++;
    cur = head[i];
    while (cur != -1) {
      ASSERT(proj[cur] == i);
      wNode[i] += 1;
      //	wAbsorb[i] += hP->wAbsorb[cur];
      for (j = ccs_ptr[cur]; j < ccs_ptr[cur + 1]; j++) {
        no = proj[ccs_ind[j]];
        if (no == -1) continue;
        if (no == i) {
          wAbsorb[i] += 1;  // hP->wEdge[j];
          continue;
        }
        if (stmp[no] == stp) {  //�Ѵ���
          int isFind = 0;
          for (k = ptr[i]; k < nEdge; k++) {
            if (adj[k] == no) {
              wEdge[k] += 1;
              isFind = 1;
              break;
            }
          }
          ASSERT(isFind == 1);
        } else {
          stmp[no] = stp;
          adj[nEdge] = no;
          wEdge[nEdge] = 1;  // hP->wEdge[j];
          nEdge++;
        }
      }
      cur = next[cur];
    }
    ptr[i + 1] = nEdge;
  }
  ASSERT(nEdge < nMaxEdge);
  delete[] head;
  delete[] next;
  wNode_p = nCol;
  wEdge_p = ccs_ptr[nCol];
  wNode_max = T_COARSEST_RELAX * (wNode_p) / THRSH_COARSEST;

  if (nEdge < nMaxEdge / 2) {
  }
  double a = 0.0, b = 0;
#ifdef _DEBUG
  ASSERT(isSymmetry(0x0));
  for (i = 0; i < nNode; i++) {
    a += wNode[i];
    b += wAbsorb[i];
  }
  ASSERT(a == nCol);
  a = 0.0;
  for (i = 0; i < nEdge; i++) {
    a += wEdge[i];
  }
  ASSERT(a + b == ccs_ptr[nCol] * 1.0);
#endif

  //	Init_Rank_List( 0x0 );
  gain_limit = 0.0;
  gain_unit_edge = FLT_MAX;
  for (i = 0; i < nNode; i++) {
    a = 0;
    for (j = ptr[i]; j < ptr[i + 1]; j++) {
      a += wEdge[j];
      gain_unit_edge = MIN(gain_unit_edge, wEdge[j]);
    }
    gain_limit = MAX(gain_limit, a);
  }
  gain_unit_edge = MAX(gain_unit_edge, 1.0);
}

/*
        v0.1	cys
                10/5/2011
*/
int MLGP::Valid(int flag) {
  int i, no, ret = GORD_OK;
  for (i = 0; i < nEdge; i++) {
    no = adj[i];
    ASSERT(no >= 0 && no < nNode);
  }

  return ret;
}

/*
        Init from CCS formant of MATRIX
        v0.1	cys
                9/27/2011
*/
MLGP::~MLGP(void) {
  if (ptr != NULL) delete[] ptr;
  if (adj != NULL) delete[] adj;
  if (tag != NULL) delete[] tag;
  if (stmp != NULL) delete[] stmp;
  if (moveto != NULL) delete[] moveto;
  if (kind != NULL) delete[] kind;

  if (wPart != NULL) delete[] wPart;
  if (wNode != NULL) delete[] wNode;
  if (wAbsorb != NULL) delete[] wAbsorb;
  if (wCluster != NULL) delete[] wCluster;
  if (wEdge != NULL) delete[] wEdge;
  if (gain != NULL) delete[] gain;

  if (project != NULL) delete[] project;

  //memset(this, 0x0, sizeof(MLGP));
}

/*
        Init from subgrahph of GRAPH
        v0.1	cys
                5/4/2012
*/
int MLGP::GetEdgeTable(int **e_n1, int **e_n2, float **e_w, int flag) {
  int nzE = nEdge / 2, nz = 0, i, j;
  int *e_1 = new int[nzE], *e_2 = new int[nzE];
  float *w = new float[nzE];
  for (i = 0; i < nNode; i++) {
    for (j = ptr[i]; j < ptr[i + 1]; j++) {
      if (adj[j] < i) continue;
      w[nz] = wEdge[j];
      e_1[nz] = i;
      e_2[nz] = adj[j];
      nz++;
    }
  }
  ASSERT(nz == nzE);
  *e_n1 = e_1;
  *e_n2 = e_2;
  *e_w = w;

  return nzE;
}

/*
        Init from subgrahph of GRAPH
        v0.1	cys
                9/29/2011
*/
MLGP *MLGP::Sub(int no, int *map, int *code, int flag) {
  double t_0 = G_TIC();
  *code = GORD_OK;
  MLGP *sub = NULL;

  int i, nC = 0, *project = new int[nNode];

  for (i = 0; i < nNode; i++) project[i] = -1;
  for (i = 0; i < nNode; i++) {
    if (map[i] == no) {
      project[i] = nC;
      nC++;
    }
  }
  sub = new MLGP(this, nC, project, MLGP_CONNECTED, 0x0);
  sub->Valid(0x0);
  MLGP::tSub += G_TIC() - t_0;
  return sub;
}
/*
        [REF]:	Comparison of coarsening schemes for multilevel graph
   partitioning

        v0.1	cys
                9/27/2011
        v0.2	cys
                5/21/2011
*/
MLGP *MLGP::ProjectDown(int *code, int flag) {
  *code = GORD_OK;
  MLGP *child = NULL, *comp = NULL;
  int i, j, no, *project = NULL, nIso = 0, cur = -1, nC = -1;
  int loop, l_max = 10, *l_project = NULL, l_nC, nDense = 0;
  float w_0 = 0.0, w_sum = 0.0, s, cls_1 = 0, cls_0 = 10000.0f, c_thrsh = 0.9f, eval, eval_0 = FLT_MAX;
  //	ASSERT( isSymmetry( 0x0 )==1 );
  project = new int[nNode];
  l_project = new int[nNode];
  for (i = 0; i < nNode; i++) {
    project[i] = -1;
    w_sum += wNode[i];
    if (ptr[i + 1] - ptr[i] > nNode / 10) {
      nDense++;
    }
  }
  s = nEdge * 1.0f / nNode;
  w_0 = 1.5f * (w_sum) / THRSH_COARSEST;
  //	ASSERT( wNode_max == w_0 );

  for (loop = 0; loop < l_max; loop++) {
    l_nC = 0;
    eval = 0;
    for (i = 0; i < nNode; i++) {
      l_project[i] = -1;
    }
    // HEM matching
    for (i = 0; i < nNode; i++) {
      //	cur = loop==0 ? list[i] : G_RAND( )%nNode;
      cur = G_RAND() % nNode;
      if (l_project[cur] != -1) continue;
      if (wNode[cur] > wNode_max) continue;
      no = -1;
      w_0 = 0.0;
      for (j = ptr[cur]; j < ptr[cur + 1]; j++) {
        if (l_project[adj[j]] != -1) continue;
        if (wNode[cur] + wNode[adj[j]] > wNode_max) continue;
        //		s = wCluster[j];
        s = wEdge[j];
        if (s > w_0) {
          w_0 = s;
          no = j;
        }
      }
      if (no != -1) {
        //			cls_1=MAX( cls_1,wCluster[no] );
        //			cls_0=MIN( cls_0,wCluster[no] );
        l_project[cur] = l_project[adj[no]] = l_nC;
        l_nC++;
      }
    }

    for (i = 0; i < nNode; i++) {
      if (l_project[i] != -1) continue;
      nIso++;
      if (0 && ptr[i + 1] - ptr[i] > 0) {
        no = -1;
        w_0 = 0;
        for (j = ptr[i]; j < ptr[i + 1]; j++) {
          if (wEdge[j] > w_0) {
            w_0 = wEdge[j];
            no = j;
          }
          //	if( wNode[adj[j]]<w_0 )	{
          //		w_0 = wNode[adj[j]];		no=j;
          //	}
        }
        ASSERT(no != -1);
        cur = adj[no];
        //	if( l_project[cur]==-1 )
        //		l_project[cur] = l_nC++;
        l_project[i] = l_project[cur];
      }
      if (l_project[i] == -1) {
        l_project[i] = l_nC++;
      }
    }
    for (i = 0; i < nNode; i++) {
      ASSERT(l_project[i] >= 0 && l_project[i] < l_nC);
      for (j = ptr[i]; j < ptr[i + 1]; j++) {
        if (l_project[i] != l_project[adj[j]]) eval += wEdge[j];
      }
    }
    if (eval < eval_0) {
      eval_0 = eval;
      nC = l_nC;
      memcpy(project, l_project, sizeof(int) * nNode);
    }
  }
  delete[] l_project;

  if (STAR_0 == -1.0) {
    STAR_0 = nC * 1.0f / nNode;
  }
  STAR = MAX(STAR, nC * 1.0f / nNode);
#ifdef _DEBUG
  //	printf( "====G=(%d,%d),ProjectDown
  // cls_1=%g%%,cls_0=%g%%\r\n",nNode,nEdge,cls_1,cls_0 );
#endif  //	if( nIso<nNode )
  if (nC < nNode * c_thrsh)
    child = new MLGP(this, nC, project, type, 0x0);
  else {
    //	type = 0x0;
    *code = GORD_DOWN_FAIL;
    child = NULL;
  }

  if (child != NULL) {
    if (1) {  // verify absorb down
      double w_sum = 0.0, w_0 = 0.0;
      for (i = 0; i < child->nNode; i++) w_sum += child->wNode[i];
      ASSERT(w_sum == wNode_p);
      w_sum = 0.0;
      w_0 = 0.0;
      for (i = 0; i < nNode; i++) w_sum += wAbsorb[i];
      for (i = 0; i < nEdge; i++) w_sum += wEdge[i];
      for (i = 0; i < child->nNode; i++) w_0 += child->wAbsorb[i];
      for (i = 0; i < child->nEdge; i++) w_0 += child->wEdge[i];
      ASSERT(w_sum == w_0 /*&& w_sum==wEdge_p*/);
    }
    //		if( 1 )
    //			child->Verify_Symmetry( 0x0 );
  }
  return child;
}

/*
        v0.1	cys
                11/16/2011
*/
int MLGP::Verify_Symmetry(int flag) {
  int i, j, k, sym, cur;
  for (i = 0; i < nNode; i++) {
    for (j = ptr[i]; j < ptr[i + 1]; j++) {
      sym = -1;
      cur = adj[j];
      for (k = ptr[cur]; k < ptr[cur + 1]; k++) {
        if (adj[k] == i) {
          sym = k;
          break;
        }
      }
      ASSERT(sym != -1);
      ASSERT(wEdge[sym] == wEdge[j]);
    }
  }
  return 0x0;
}
/*
        v0.1	cys
                9/28/2011
*/
void MLGP::DumpInfo(int flag) {
  float w_N = 0.0, w_E = 0.0, w_N_1 = 0.0, w_N_0 = FLT_MAX, w_E_1 = 0.0, w_E_0 = FLT_MAX, w;
  int i;
  for (i = 0; i < nNode; i++) {
    w = wNode[i];
    w_N += w;
    w_N_0 = MIN(w, w_N_0);
    w_N_1 = MAX(w, w_N_1);
  }
  for (i = 0; i < nEdge; i++) {
    w = wEdge[i];
    w_E += w;
    w_E_0 = MIN(w, w_E_0);
    w_E_1 = MAX(w, w_E_1);
  }
  printf("(%d,%d) wN=(%6.1f,%.0f,%.0f) wE=(%6.1f,%.0f,%.0f)\r\n", nNode, nEdge, w_N, w_N_0, w_N_1, w_E, w_E_0, w_E_1);
}

/*
        v0.1	cys
                9/28/2011
*/
int MLGP::Init_Rank_List(int flag) {
  int ret = GORD_OK;

  if (nNode > 1 && BIT_TEST(type, MLGP_CONNECTED)) {
  } else
    goto END;

  //	wN_max = 1.5*(wNode_p)/THRSH_COARSEST;
  //	wCluster = new float[nEdge];
  int i, j, *temp = new int[nNode], nz_1 = 0, nz_0 = nNode + 1;

  float w = 0.0, s;
  double gain;
  gain_limit = 0.0;
  for (i = 0; i < nNode; i++) {
    //	if( wNode[i]>wN_max )
    //		continue;
    w = 0.0;
    gain = 0;
    for (j = ptr[i]; j < ptr[i + 1]; j++) {
      //	if( wNode[i]+wNode[adj[j]]>wN_max )
      //		continue;
      s = wEdge[j];
      gain += wEdge[j];
      /*	c_n = wNode[i]+wNode[adj[j]];		ASSERT( c_n>=2 );
              c_e = wAbsorb[i]+wAbsorb[adj[j]]+wEdge[j]*2;
              s = c_e*100.0/(c_n*(c_n-1));*/
      if (wCluster != NULL) {
        wCluster[j] = s;
        w = MAX(w, wCluster[j]);
      }
    }
    //		w = wNode[i];
    ASSERT(w >= 0);
    temp[i] = (int)(w + 1);
    nz_1 = MAX(nz_1, temp[i]);
    nz_0 = MIN(nz_0, temp[i]);
  }
  /*
          list.Init( nz_1-nz_0+1,nNode,0x0 );
          for( i = 0; i < nNode; i++ )	{
                  pos = temp[i]-nz_0;
                  list.Insert( pos,i,0x0 );
          }
          list.Sequence( 0x0 );

          for( i = 0; i < nNode; i++ )	{
                  cur = list[i];
                  ASSERT( temp[cur]>0 );
                  temp[cur]=-1;
          }*/
  delete[] temp;
END:
  return ret;
}

/*
        v0.1	cys
                9/29/2011
*/
extern "C" int _ccs_component(int nCol, int *ptr, int *ind, int *map, int flag) {
  int i, j, seed, nCom = 0, cur, top = 0, *stack, no;

  stack = new int[nCol];
  for (i = 0; i < nCol; i++) {
    map[i] = -1;
  }
  for (seed = 0; seed < nCol; seed++) {
    if (map[seed] != -1) continue;
    nCom++;
    ASSERT(top == 0);
    stack[top++] = seed;
    while (top > 0) {
      cur = stack[--top];
      map[cur] = nCom;
      for (j = ptr[cur]; j < ptr[cur + 1]; j++) {
        no = ind[j];
        if (map[no] != -1) {
          ASSERT(map[no] == nCom);
        } else {
          map[no] = nCom;
          stack[top++] = no;
        }
      }
    }
  }
  delete[] stack;

  return nCom;
}

/*
        v0.1	cys
                9/29/2011
*/
int MLGP::SplitComponent(hMLGP **arrGraph, int *map, int flag) {
  double t_0 = G_TIC();
  int nCom = 0, i, code;
  nCom = _ccs_component(nNode, ptr, adj, map, 0x0);
  *arrGraph = new hMLGP[nCom];
  if (nCom == 1) {
    (*arrGraph)[0] = this;
  } else {
    for (i = 0; i < nCom; i++) {
      (*arrGraph)[i] = this->Sub(i + 1, map, &code, MLGP_CONNECTED);
    }
  }
  MLGP::tSplit += G_TIC() - t_0;
  return nCom;
}

/*
        2 components + Veterx Separator
        v0.1	cys
                11/18/2011
*/
int MLGP::Partition_3(int flag) {
  int ret = GORD_VSEP_FAIL, type = 0;
  switch (type) {
    case 0:
      ret = P3_GG(flag);
      break;
    case 1:
      ret = P3_Spectral(flag);
      break;
    case 2:
      ret = P3_Karger(flag);
      break;
    case 3:
      ret = P3_GreedyGrow(NULL, NULL, NULL, 0x0);
      break;
    default:
      ASSERT(0);
      break;
  }
  return ret;
}

/*
        Graph Growing
        v0.1	cys
                9/29/2011
*/
int MLGP::P3_GG(int flag) {
  if (nNode < 2) return GORD_TOO_SMALL;
  int ret = GORD_VSEP_FAIL;
  int i = 0, j, no, *tag_0, T_s, T_1, T_2, seed, cur, nRand, l, *stack, top;
  float eval, eval_0 = FLT_MAX, s_0, s_1 = 0.0;

  stack = new int[nNode];
  tag_0 = new int[nNode];
  s_0 = nEdge * 100.0f / (nNode * (nNode - 1));
  if (1) {
    for (i = 0; i < nNode; i++) {
      if (ptr[i + 1] - ptr[i] <= 2) s_1++;
      /*		stack[nNode-1-i] = list[i];
                      if( wNode[i]==1 )
                              continue;
                      s = wAbsorb[i]*100.0/(wNode[i]*(wNode[i]-1));
                      if( s_1<s )
                              s_1=s;*/
    }
  }
  //	srand( (unsigned)time( NULL ) );
  nRand = MIN(20, nNode / 2);
  for (l = 0; l < nRand; l++) {
    seed = G_RAND() % nNode;
    //	seed = list[l];
    T_1 = 0;
    for (i = 0; i < nNode; i++) tag[i] = 2;
    top = 0;
    stack[top++] = seed;
    tag[seed] = 1;
    while (top > 0) {
      cur = stack[--top];
      for (j = ptr[cur]; j < ptr[cur + 1]; j++) {
        no = adj[j];
        if (tag[no] == 2) {
          tag[no] = 1;
          T_1++;
          stack[top++] = no;
        }
        if (T_1 > nNode / 2) goto NEXT;
      }
    }
  NEXT:
    T_s = 0;
    T_1 = 0;
    T_2 = 0;
    for (i = 0; i < nNode; i++) {
      if (tag[i] == 0) continue;
      for (j = ptr[i]; j < ptr[i + 1]; j++) {
        no = adj[j];
        if (tag[i] != tag[no] && (tag[no] != 0)) {
          tag[i] = 0;
          T_s++;
          break;
        }
      }
      if (tag[i] == 1)
        T_1++;
      else if (tag[i] == 2) {
        T_2++;
      }
    }
    if (/*T_s==0 ||*/ T_1 == 0 || T_2 == 0) continue;

    ASSERT(T_s + T_1 + T_2 == nNode);
    if (0) {
      Refine_3(MLGP::VFM_REFINE, NULL, 0x0);
    }
    eval = Eval(tag, 0x0);
    if (eval < eval_0) {
      eval_0 = eval;
      memcpy(tag_0, tag, sizeof(int) * nNode);
      ret = GORD_OK;
    }
  }
  if (ret == GORD_OK) {
    memcpy(tag, tag_0, sizeof(int) * nNode);
    eval = Eval(NULL, 0x0);
    ASSERT(eval == eval_0);
#ifdef _DEBUG
    //	printf( "====Partition
    // G=(%d,%d),(%d,%d,%d),r=%g%%\r\n",nNode,nEdge,T_s,T_1,T_2,nEdge*100.0/(nNode*(nNode-1))
    //);
#endif
    for (i = 0; i < nMaxPart + 1; i++) wPart[i] = 0.0;
    for (i = 0; i < nNode; i++) {
      no = tag[i];
      ASSERT(no >= 0 && no <= nMaxPart);
      wPart[no] += wNode[i];
    }
    for (i = 1; i < nMaxPart + 1; i++) ASSERT(wPart[i] > 0.0);
  }
  delete[] tag_0;
  delete[] stack;

  return ret;
}

int MLGP::P3_Karger(int flag) {
  int ret = GORD_OK, i, nC = 0, nG = 2;
  float a = log(nNode * 1.0f), cut_0 = FLT_MAX;
  int nLoop = int(1.1 * a * a + 1.0);
  int *cut = NULL, *cut_set = new int[nNode];
  MLGP *hBest = NULL, *hG = NULL, *hG_0;
  int nCom = _ccs_component(nNode, ptr, adj, cut_set, 0x0);
  ASSERT(nCom == 1);
  delete[] cut_set;

  for (i = 0; i < nLoop; i++) {
    //	Karger_MinCut( this,&cut_0,&hBest,0x0 );
  }
  if (hBest != NULL) {
    hG = hBest;
    ASSERT(hG != this);
    while (hG != this) {
      hG_0 = hG->parent;
      hG->ProjectUp(0x0);
      delete hG;
      hG = hG_0;
    }
  }

  return ret;
}

/*
        Multilevel Algorithms for Partitioning Power-Law Graphs
        v0.1	cys
                11/15/2011
*/
int MLGP::P3_Spectral(int flag) {
  int ret = GORD_VSEP_FAIL;
  /*	if( nNode<2 )
                  return GORD_TOO_SMALL;
          int
  ret=GORD_VSEP_FAIL,i,j,no,T_s=0,T_1=0,T_2=0,info,lwork=nNode*nNode,loop,check;
          EIGEN_DOUBLE *mat=new EIGEN_DOUBLE[nNode*nNode],e_sum=0.0,*w=new
  EIGEN_DOUBLE[nNode],s,w_pre; EIGEN_DOUBLE *work=new EIGEN_DOUBLE[lwork];
          EIGEN_DOUBLE *D=new EIGEN_DOUBLE[nNode];
          EIGEN_DOUBLE fiedler,f_thrsh=0.0,c_n,c_e;

          char jobz='V',uplo='L';//,matrix_order=LAPACK_ROW_MAJOR;
          for( i = nNode*nNode-1; i >= 0; i-- )	mat[i] = 0.0;
          for( i = 0; i < nNode; i++ )	{		//����wCluster
                  e_sum = 0.0;
                  for( j = ptr[i]; j < ptr[i+1]; j++ )	{
                          s = wEdge[j];
                          no = adj[j];
                          c_n = wNode[i]+wNode[no];		ASSERT( c_n>=2
  ); c_e = wAbsorb[i]+wAbsorb[no]+wEdge[j]*2; s = c_e*100.0/(c_n*(c_n-1));
                          wCluster[j] = s;
                          e_sum += wCluster[j];
                  }
                  D[i]=sqrt(e_sum);
          }
  //	Verify_Symmetry( 0x0 );
          for( i = 0; i < nNode; i++ )	{
                  e_sum = 0.0;
                  for( j = ptr[i]; j < ptr[i+1]; j++ )	{
                          no = adj[j];
                          mat[i*nNode+no] = -wCluster[j]/D[i]/D[no];
                  }
                  mat[i*nNode+i]=1.0;
  //		mat[i*nNode+i]=D[i]*D[i];
          }
          if( 1 )	{
                  for( i = 0; i < nNode; i++ )	{
                          e_sum = 0.0;
                          for( j = 0; j < nNode; j++ )	{
                                  e_sum += mat[i*nNode+j];
                                  s = mat[i*nNode+j]-mat[j*nNode+i];
                                  ASSERT( fabs(s)<DBL_EPSILON );
                          }
                  //	ASSERT( fabs(e_sum)<DBL_EPSILON );
                  }
          }
  //	lapack_int code = LAPACKE_ssyev( matrix_order, jobz, uplo, nNode, mat,
  nNode, w ); dsyev( &jobz, &uplo, &nNode, mat,&nNode, w,work, &lwork,&info );
          ASSERT( info==0 );
          delete[] work;
          float eval_0=FLT_MAX,eval;
          int *tag_0 = new int[nNode];
          w_pre = w[0];
          for( loop = 1; loop < nNode; loop++ )	{
                  s = w[loop]-w_pre;
                  w_pre = w[loop];
                  if( fabs(s)<DBL_EPSILON )
                          continue;
                  for( i = 0; i < nNode; i++ )	{
                          fiedler = mat[loop+i*nNode]/D[i];
                          tag[i]= fiedler<f_thrsh ? 1 : 2;
                  }
                  T_s=0,	T_1=0,	T_2=0;
                  for( i = 0; i < nNode; i++ )	{
                          if( tag[i]==0 )
                                  continue;
                          for( j = ptr[i]; j < ptr[i+1]; j++ )	{
                                  no = adj[j];
                                  if( tag[i]!=tag[no] && (tag[no]!=0) )
                                  {	tag[i]=0;		T_s++;
  break;	}
                          }
                          if( tag[i]==1 )
                                  T_1++;
                          else if( tag[i]==2 )	{
                                  T_2++;
                          }
                  }
                  eval = Eval( tag,0x0 );
                  if( eval<eval_0 )	{
                          eval_0 = eval;		memcpy(
  tag_0,tag,sizeof(int)*nNode ); ret = GORD_OK;
                  }
          }
          memcpy( tag,tag_0,sizeof(int)*nNode );
          delete[] tag_0;

          delete[] mat;	delete[] w;				delete[] D;*/

  return ret;
}

/*
        [REF]:	Linear time heuristic for improving network partitions
        v0.1	cys
                9/30/2011
*/
int MLGP::Refine_3(int type, int *temp, int flag) {
  int ret = GORD_OK, loop = 0, l_max = 10, *cut = temp, nCut = 0, i, j;
  float w_0, eval;
  ASSERT(nMaxPart == 2);
  if (0) {
    for (i = 0; i < nNode; i++) {
      for (j = ptr[i]; j < ptr[i + 1]; j++) {
        if (tag[adj[j]] != tag[i]) {
          cut[nCut++] = i;
          break;
        }
      }
    }
    for (i = 0; i < nCut; i++) tag[cut[i]] = 0;
  }
  switch (type) {
    case VFM_REFINE:
      VFM_Refine(&eval, 0x0);
      break;
    default:
      break;
  };
  if (1) {
    w_0 = wPart[0];
    while (loop++ < l_max) {
      if (Flow_Refine(0x0) != GORD_OK) break;
      if (wPart[0] >= w_0) break;
      w_0 = wPart[0];
    }
  }
  //	if( !BIT_TEST( flag,REFINE_APART_LOOP ) )
  //		Apart_Refine( 0x0 );
  return ret;
}

/*
        v0.1	cys
                9/30/2011
*/
float MLGP::Eval(int *part_0, int flag,bool  isRecalc = true, float w1, float w0) {
  double t_1 = G_TIC();
  if(isRecalc) wCutEdge=-1;
  //	ASSERT( sepa==SEPARATOR_V );
  float eval = 0.0, alpha = 0.75, beta = 0.25, V_avg = 0.0, nV[3] = {0.0, 0.0, 0.0};
  int i, j, S_n = 0, no, w_sum = 0, e_type = 1;
  bool isDump = BIT_TEST(flag, MLGP_DUMP);
  int *part = part_0 == NULL ? tag : part_0;

  if (flag == 1) {
    for (i = 0; i < nMaxPart + 1; i++) w_sum += (int)wPart[i];
#ifdef _DEBUG
      /*	float *T_w = new float[nMaxPart+1];
              for( i = 0; i < nMaxPart+1; i++ )	T_w[i]=0;
              for( i = 0; i < nNode; i++ )	{
                      no = part[i];
                      T_w[no]+=wNode[i];
              }
              for( i = 0; i < nMaxPart+1; i++ )
                      ASSERT( T_w[i]==wPart[i] );
              delete[] T_w;*/
#endif
  } else {
    for (i = 0; i < nMaxPart + 1; i++) {
      wPart[i] = 0;
    }
    for (i = 0; i < nNode; i++) {
      no = part[i];
      w_sum += (int)(wNode[i]);
      switch (no) {
        case 0:
          //			S_n+=wNode[i];
          //			break;
        case 1:
        case 2:
          wPart[no] += wNode[i];
          nV[no]++;
          break;
        default:
          ASSERT(0x0);
          break;
      }
    }
  }

  if (e_type == 0) {
    S_n = (int)(wPart[0]);
    V_avg = (w_sum - S_n) * 1.0f / nMaxPart;
    eval = alpha * S_n * S_n;
    //	eval = alpha*nV[0]*nV[0];
    for (i = 0; i < nMaxPart; i++) {
      eval += beta * (wPart[i + 1] - V_avg) * (wPart[i + 1] - V_avg);
    }
  } else {
    double w_1 = 0.0, w_0 = DBL_MAX, w_2 = 1.0;
    if (isRecalc) {  //summary of CutEdge
      for (i = 0; i < nNode; i++) {
        if (tag[i] != 0) continue;
        /*	w_2 += ptr[i+1]-ptr[i];*/
        for (j = ptr[i]; j < ptr[i + 1]; j++) {
          if (tag[adj[j]] != 0) w_2 += wEdge[j];
        }
      }
      if (wCutEdge > 0 && wCutEdge != w_2 ) {
        ASSERT(wCutEdge == w_2);
      }
      wCutEdge = w_2;
      // w_2 = 1.0f + w_2 / nEdge * 5;
    } else {
    }
    

    w_2 = 1.0f + wCutEdge / nEdge * 5;
    for (i = 0; i < nMaxPart; i++) {
      w_1 = MAX(w_1, wPart[i + 1]);
      w_0 = MIN(w_0, wPart[i + 1]);
    }
    //		if( w_0 == 0.0 )
    if (w_1 > w_sum * 0.999)
      eval = FLT_MAX;
    else if (w_0 == 0) {
      eval = FLT_MAX;
    } else
      eval = (float)((wPart[0] + 1.0) * (1.0 + 0.05 * w_1 / w_0) / w_2);
  }

  if (isDump) {
    printf("G=(%d,%d),eval=%.1f,sep=%d,iter=%d", nNode, nEdge, eval, S_n, iter);
    printf("\r\n");
  }
  if(isRecalc){
    MLGP::tEval += G_TIC() - t_1;     MLGP::nEval++;      
  }
    
  return eval;
}

/*
        v0.1	cys
                9/29/2011
*/
int MLGP::CheckPartition(int flag) {
  if (sepa == SEPARATOR_E) return GORD_OK;

  int i, j, no, ret = GORD_OK, cur, nzE = 0;

  for (i = 0; i < nNode; i++) {
    cur = tag[i];
    if (cur == 0) continue;
    for (j = ptr[i]; j < ptr[i + 1]; j++) {
      no = adj[j];
      if (cur != tag[no]) {
        ASSERT(tag[no] == 0);
        nzE++;
      }
    }
  }
  //	printf( " nzE=%d:",nNode,nEdge,nzE );
  //	float eval = Eval( NULL, MLGP_DUMP );

  return ret;
}

/*
        gain=|s_0|-|s_1|

        v0.1	cys
                10/3/2011
*/
int MLGP::UpdateGain_V(int no, int *x, float base, int flag) {
  int j, ret = GORD_OK, neibor;
  float g_1 = wNode[no], g_2 = wNode[no];
  for (j = ptr[no]; j < ptr[no + 1]; j++) {
    neibor = adj[j];
    if (tag[neibor] == 1) g_2 -= wNode[neibor];
    if (tag[neibor] == 2) g_1 -= wNode[neibor];
  }

  if (wPart[1] < wPart[2] / 3 && g_1 > 0) {
    *x = 1;
    gain[no] = g_1;
  } else if (wPart[1] > wPart[2] * 3 && g_2 > 0) {
    *x = 2;
    gain[no] = g_2;
  } else {  // P_1,P_2�����൱
    if (g_1 > g_2) {
      *x = 1;
      gain[no] = g_1;
    } else if (g_1 < g_2) {
      *x = 2;
      gain[no] = g_2;
    } else {
      if (wPart[1] <= wPart[2]) {
        *x = 1;
        gain[no] = g_1;
      } else {
        *x = 2;
        gain[no] = g_2;
      }
    }
  }

  gain[no] = gain[no] + base;
  return ret;
}

/*
        v0.1	cys
                9/30/2011
*/
int MLGP::Apart_Refine(int flag) {
  int ret = GORD_OK, i, j, eval_0, eval_1;
  int *cut = new int[nNode], *map = new int[nNode], nCut_0 = 0, nCut_1 = 0, nApart = 0, code, no;
  int *tag_0 = new int[nNode];
  float *wcut = new float[nNode], *u_adj = new float[nNode], u_sum = 0.0, u_avg = 0.0;

  memcpy(tag_0, tag, sizeof(int) * nNode);
  eval_0 = Eval(NULL, 0x0);
  for (i = 0; i < nNode; i++) {
    map[i] = 1;
    if (tag[i] != 0) continue;
    u_adj[nCut_0] = 0.0;
    for (j = ptr[i]; j < ptr[i + 1]; j++) {
      if (tag[adj[j]] == 0) continue;
      u_adj[nCut_0] += wNode[adj[j]];
    }
    u_sum += u_adj[nCut_0];
    //	ASSERT( u_adj[nCut]>0.0 );
    wcut[nCut_0] = wNode[i];
    map[i] = -1;
    cut[nCut_0++] = i;
  }
  u_avg = u_sum / nCut_0;
  _quick_sort_2(u_adj, cut, 0, nCut_0 - 1, 0x0);
  for (i = 0; i < nCut_0; i++) {
    if (u_adj[i] > u_avg * 2.0) break;
    map[cut[i]] = 1;
  }
  nApart = nCut_0 - i - 1;

  MLGP *hSub = Sub(1, map, &code, MLGP_CONNECTED);
  hSub->BiSection(NULL, NULL, NULL, 0x0, 0x0);
  nCut_1 = 0;
  for (i = 0; i < nNode; i++) {
    no = hSub->project[i];
    if (no == -1) {
      tag[i] = 0;
    } else {
      tag[i] = hSub->tag[no];
    }
    if (tag[i] == 0) {
      nCut_1++;
    }
  }
  delete hSub;
  eval_1 = Eval(NULL, 0x0);
  if (eval_1 >= eval_0) {
    memcpy(tag, tag_0, sizeof(int) * nNode);
  } else {
    int i = 0;
  }

  delete[] cut;
  delete[] map;
  delete[] wcut, delete[] u_adj;
  delete[] tag_0;
  return ret;
}

/*
        v0.1	cys
                9/30/2011
*/
int MLGP::MoveSeparator(int flag) {
  return -1;
  int ret = GORD_OK, nMove = 2, l, i, j, no, *T_tag = new int[nNode], x;
  float e_0 = Eval(NULL, 0x0), e_1;

  memcpy(T_tag, tag, sizeof(int) * nNode);
  for (l = 0; l < nMove; l++) {
    x = wPart[1] < wPart[2] ? 2 : 1;
    for (i = 0; i < nNode; i++) {  // absort isolated partion to separator
      if (tag[i] != 0) continue;
      for (j = ptr[i]; j < ptr[i + 1]; j++) {
        no = adj[j];
        if (tag[no] == x) {
          T_tag[no] = 0;
        }
      }
      T_tag[i] = x == 1 ? 2 : 1;
    }
    memcpy(tag, T_tag, sizeof(int) * nNode);
    e_1 = Eval(NULL, 0x0);
  }
  return ret;
}

/*
        ������Separator��صĹ����ڵ�

        v0.1	cys
                11/18/2011
*/
void MLGP::Separator_ISO(int, int, int flag) {
  int i, j, nIso = 0;
  bool isIso, isChange = false;

  for (i = 0; i < nNode; i++) {
    if (tag[i] == 0) {
      isIso = true;
      for (j = ptr[i]; j < ptr[i + 1]; j++) {
        if (tag[adj[j]] != 0) {
          isIso = false;
          break;
        }
      }
      if (isIso) {
        //	tag[i] = wPart[1]<wPart[2] ? 1 : 2;
        //	isChange=true;
      }
    }
  }
  if (isChange) Eval(NULL, 0x0);
  for (i = 0; i < nNode; i++) {
    if (tag[i] == 0) continue;
    isIso = true;
    for (j = ptr[i]; j < ptr[i + 1]; j++) {
      if (tag[adj[j]] != 0) {
        isIso = false;
        break;
      }
    }
    if (isIso) {
      //	tag[i] = 0;
      tag[i] = wPart[1] < wPart[2] ? 1 : 2;
      nIso++;
    }
  }
  if (nIso > 0) {
    Eval(NULL, 0x0);
  }
}

/*
        v0.1	cys
                4/24/2012
*/
int MLGP::Flow_Refine(int flag) {
  MaxFlow_push flow;
  int ret = GORD_FLOW_REFINE_FAIL, i, j, k, no, cur, next, part = -1, nz, f_0 = 0, f_1 = 0, f_2 = 0, nS = 0, nX = 0;
  int *B_s = new int[nNode], *B_x = new int[nNode], *cut = new int[nNode + 2], *map = new int[nNode];

  for (i = 0; i < nNode; i++) {
    map[i] = -1;
  }
  part = wPart[1] < wPart[2] ? 2 : 1;
  nz = 0;
  for (i = 0; i < nNode; i++) {
    no = tag[i];
    if (no != 0) continue;
    map[i] = nS;
    B_s[nS++] = i;
    for (j = ptr[i]; j < ptr[i + 1]; j++) {
      cur = adj[j];
      if (tag[cur] == part && map[cur] == -1) {
        B_x[nX++] = cur;
        map[cur] = 1;
      }
      for (k = ptr[cur]; k < ptr[cur + 1]; k++) {
        if (tag[adj[k]] == part && map[adj[k]] == -1) {
          B_x[nX++] = adj[k];
          map[adj[k]] = 1;
        }
      }
      nz += ptr[cur + 1] - ptr[cur];
    }
    nz += ptr[i + 1] - ptr[i];
  }
  for (i = 0; i < nX; i++) {
    map[B_x[i]] = i + nS;
    B_s[nS + i] = B_x[i];
  }
  //	ASSERT( nS+nX==nV[0] );
  int F_nV = nS + nX + 2;
  if (F_nV < 200) goto END;
  nz = MIN(nz * 2 + 2 * (nS + nX), ptr[nNode] * 2 + 2 * (nS + nX));
  int *F_ptr = new int[F_nV + 1], *F_adj = new int[nz], *F_cap = new int[nz], F_nz = 0;
  F_ptr[0] = 0;
  for (i = 0; i < nS; i++) {  // source node
    F_adj[F_nz] = i + 1;
    F_cap[F_nz] = wNode[B_s[i]];
    F_nz++;
    f_0 += wNode[B_s[i]];
  }
  F_ptr[1] = F_nz;
  for (i = 0; i < nS + nX; i++) {
    cur = B_s[i];
    if (i < nS) {  // separator
      F_adj[F_nz] = 0;
      F_cap[F_nz] = 0;
    } else {  // X part
      F_adj[F_nz] = nS + nX + 1;
      F_cap[F_nz] = wNode[cur];
    }
    F_nz++;
    for (j = ptr[cur]; j < ptr[cur + 1]; j++) {
      next = map[adj[j]];
      if (next == -1) continue;
      if ((i < nS && next < nS) || (i >= nS && next >= nS)) continue;
      F_adj[F_nz] = next + 1;
      F_cap[F_nz] = i < nS ? CAPCITY_INFINITE : 0;
      F_nz++;
    }
    F_ptr[i + 2] = F_nz;
  }
  for (i = nS; i < nS + nX; i++) {  // sink node
    F_adj[F_nz] = i + 1;
    F_cap[F_nz] = 0;
    F_nz++;
  }
  F_ptr[nS + nX + 2] = F_nz;
  ASSERT(F_nz <= nz);
  flow.Init(F_nV, F_ptr, F_adj, F_cap, 0x0);
  f_1 = flow.Push_Relabel(0x0, 0x0);
  if (f_1 < f_0) {
    flow.GetCut(cut, NULL, 0x0);
    for (i = 0; i < nS + nX; i++) {
      cur = B_s[i];
      if (cut[i + 1] == 1) {
        tag[cur] = 0;
        f_2 += wNode[cur];
      } else {
        if (tag[cur] == 0) {
          tag[cur] = part == 2 ? 1 : 2;
        } else {
        }
      }
    }
    ASSERT(f_2 == f_1);
  }
  CheckPartition(0x0);
  Eval(NULL, 0x0);
  delete[] F_ptr, delete[] F_adj, delete[] F_cap;
  ret = GORD_OK;
END:
  delete[] map;
  delete[] cut;
  delete[] B_s;
  delete[] B_x;
  return ret;
}

/*
        [REF]:	1 Linear time heuristic for improving network partitions

        v0.1	cys
                9/30/2011
*/
int MLGP::VFM_Refine(float *vfm_eval, int flag) {
  double t_1 = G_TIC(), t_0 = G_TIC(), t_gain = 0, t_2;
  int sepa_old = sepa;
  if (sepa != SEPARATOR_V) sepa = SEPARATOR_V;
  int ret = GORD_OK, i, j, k, no, cand = -1, nMove = 0, nMove_0;
  int *set = NULL, *elig = NULL, *T_tag = new int[nNode], *modify = new int[nNode], stp = 1, nM = 0, nIso = 0;
  int *pass = new int[nNode];
  float eval_0 = 0, eval, G_thrsh = 0.0, w_max = 0, base = w_max / 2;
  bool isBetter = false, isIso;

  nVFM++;
  //	if( nVFM==278 )
  //		nVFM=278;
  elig = new int[nNode];  // elig=(0-lock,1-move to V1,2-move to V2)
  SORT_BUCKET moves;
  for (i = 0; i < nNode; i++) {
    //	w_max+=wNode[i];
    eval = 0.0;
    for (j = ptr[i]; j < ptr[i + 1]; j++) {
      eval += wNode[adj[j]];
    }
    w_max = MAX(eval, w_max); /**/
                              //	w_max =	MAX( w_max,ptr[i+1]-ptr[i] );
    stmp[i] = stp;
    pass[i] = -1;
  }
  memcpy(T_tag, tag, sizeof(int) * nNode);
  //	base=w_max/5;
  //	base=w_max<10 ? w_max/2 : w_max<100 ? w_max/4 : w_max/10;
  base = w_max;
  moves.Init((int)(base + w_max + 1), nNode, 0x0);
  eval_0 = Eval(NULL, 0x0);
  iter = 0;
  t_2 = G_TIC();
  do {
    iter++;
    isBetter = false;
    nMove = 0, nMove_0 = 0;
    //		t_ = clock( );
    if (!BIT_TEST(flag,VFM_NO_ABSORB)) {  // absort isolated partion to separator
      for (i = 0; i < nNode; i++) {
        if (tag[i] == 0) {
          isIso = true;
          for (j = ptr[i]; j < ptr[i + 1]; j++) {
            if (tag[adj[j]] != 0) {
              isIso = false;
              break;
            }
          }
          if (isIso) {
            //	tag[i] = wPart[1]<wPart[2] ? 1 : 2;
            tag[i] = i % 2 + 1;
          }
        } else {
          isIso = true;
          for (j = ptr[i]; j < ptr[i + 1]; j++) {
            if (tag[adj[j]] != 0) {
              isIso = false;
              break;
            }
          }
          if (isIso) {
            tag[i] = 0;
            nIso++;
          }
        }
      }
    }
    Eval(NULL, 0x0);    
    for (i = 0; i < nNode; i++) {
      elig[i] = -1;
      if (tag[i] != 0) continue;
      UpdateGain_V(i, elig + i, base, 0x0);
      if (gain[i] > G_thrsh) {
        moves.Insert((int)(gain[i]), i, 0x0);
        nMove_0++;
      }      
    }
    //		t_gain += clock( )-t_;
    while (!moves.isEmpty(0x0)) {
      stp++;
      nMove++;
      nM = moves.FirstRankList(modify, 0x0);
      ASSERT(nM > 0);
      cand = nM > 1 ? modify[G_RAND() % nM] : modify[0];
      if (nM > 1) {  //�����������ƽ����ԣ�����֣�
        /*	for( i = 0; i < nM; i++ )	{
                        int x = wPart[1]<wPart[2] ? 2 : 1;
                        if( elig[modify[i]]==x )
                        {	cand = modify[i];		break;	}
                }*/
      }
      ASSERT(pass[cand] < iter);
      pass[cand] = iter;
      ASSERT(tag[cand] == 0);
      ASSERT(gain[cand] >= G_thrsh);
      ASSERT(elig[cand] == 1 || elig[cand] == 2);
      tag[cand] = elig[cand] == 1 ? 1 : 2;
      wPart[0] -= wNode[cand];
      wPart[tag[cand]] += wNode[cand];
      for (j = ptr[cand]; j < ptr[cand + 1]; j++) {  // update partition
        no = adj[j];      
        if (tag[no] != tag[cand]) {
          wPart[tag[no]] -= wNode[no];
          if (elig[no] != 0) {
            tag[no] = 0;
          } else {  //�������
            tag[no] = 0;
            elig[no] = -1;
          }
          wPart[tag[no]] += wNode[no];
        }
        if(isGreedyCutEdge){
          if (tag[no] == 0)
            wCutEdge += wEdge[j];
          else {
            wCutEdge -= wEdge[j];
          }
        }
      }
      
      elig[cand] = 0;
      moves.Remove((int)(gain[cand]), cand, 0x0);
      eval = Eval(NULL, 1,!isGreedyCutEdge, -1, -1);      
      if (eval < eval_0) {
        if(isGreedyCutEdge)
          eval = Eval( NULL,1 );      //w_CutEdge is only approximation
        if (eval < eval_0) {
          eval_0 = eval;
          memcpy(T_tag, tag, sizeof(int) * nNode);
          isBetter = true;
        }
      }
      nM = 0;
      for (j = ptr[cand]; j < ptr[cand + 1]; j++) {
        no = adj[j];
        if (tag[no] != 0) continue;
        if (stmp[no] == stp) continue;
        stmp[no] = stp;
        modify[nM++] = no;
        if (1) {
          for (k = ptr[no]; k < ptr[no + 1]; k++) {
            if (tag[adj[k]] != 0) continue;
            if (stmp[adj[k]] == stp) continue;
            stmp[adj[k]] = stp;
            modify[nM++] = adj[k];
          }
        }
      }
      for (j = 0; j < nM; j++) {
        no = modify[j];
        ASSERT(tag[no] == 0);
        if (pass[no] == iter) continue;
        if (gain[no] > G_thrsh && elig[no] != -1) moves.Remove((int)(gain[no]), no, 0x0);
        UpdateGain_V(no, elig + no, base, 0x0);
        if (gain[no] > G_thrsh)
          moves.Insert((int)(gain[no]), no, 0x0);
        else {
        }
      }
    }
    memcpy(tag, T_tag, sizeof(int) * nNode);
    eval = Eval(NULL, 0x0);
  } while (isBetter && iter < 10);
  
  //	Separator_ISO( 0x0,0x0,0x0 );

  delete[] T_tag;
  delete[] elig;
  delete[] modify;
  delete[] pass;
  *vfm_eval = Eval(NULL, 0x0);
  ASSERT(*vfm_eval == eval_0);
  sepa = sepa_old;
  // printf("\tRefine:	time=%.3g,gain=%.3g,move=%d\r\n", G_TIC() - t_0, t_gain, nMove);
  MLGP::tRefine += G_TIC() - t_1;
  return ret;
}

/*
        v0.1	cys
                11/4/2011
*/
int MLGP::VFM_Perturb(int flag) {
  int ret = GORD_OK, i, j;
  int *cut = new int[nNode], nCut_0 = 0, nCut_1 = 0, nApart = 0, no;
  int *tag_0 = new int[nNode];
  float *wcut = new float[nNode], *u_adj = new float[nNode], u_sum = 0.0, u_avg = 0.0;

  memcpy(tag_0, tag, sizeof(int) * nNode);
  for (i = 0; i < nNode; i++) {
    if (tag[i] != 0) continue;
    u_adj[nCut_0] = 0.0;
    for (j = ptr[i]; j < ptr[i + 1]; j++) {
      if (tag[adj[j]] == 0) continue;
      u_adj[nCut_0] += wNode[adj[j]];
    }
    u_sum += u_adj[nCut_0];
    //	ASSERT( u_adj[nCut]>0.0 );
    cut[nCut_0++] = i;
  }
  u_avg = u_sum / nCut_0;
  _quick_sort_2(u_adj, cut, 0, nCut_0 - 1, 0x0);
  Eval(NULL, 0x0);
  int move_to = wPart[1] < wPart[2] ? 1 : 2, cur;
  for (i = 0; i < nCut_0; i++) {
    if (u_adj[i] > u_avg * 1.5) break;
    cur = cut[i];
    tag[cur] = move_to;
    for (j = ptr[cur]; j < ptr[cur + 1]; j++) {
      no = tag[adj[j]];
      if (no == 0 || no == move_to) continue;
      tag[adj[j]] = 0;
    }
  }
  Eval(NULL, 0x0);
  if (wPart[1] == 0 || wPart[2] == 0) {
    memcpy(tag, tag_0, sizeof(int) * nNode);
    ret = GORD_PERTUB_FAIL;
    goto END;
  }
  CheckPartition(0x0);
END:
  delete[] cut;
  delete[] tag_0;
  delete[] wcut;
  delete[] u_adj;
  return ret;
}

/*
        v0.1	cys
                9/29/2011
*/
int MLGP::ProjectUp(int flag) {
  ASSERT(parent != NULL);
  int i, no, ret = GORD_OK;
  ASSERT(parent->nMaxPart == nMaxPart);

  parent->nMaxPart = nMaxPart;
  for (i = 0; i < nMaxPart + 1; i++) parent->wPart[i] = 0.0;
  for (i = 0; i < parent->nNode; i++) {
    ASSERT(project[i] != -1);
    no = tag[project[i]];
    ASSERT(no >= 0 && no <= nMaxPart);
    parent->tag[i] = no;
    parent->wPart[no] += parent->wNode[i];
  }

  return ret;
}

/*
        v0.1	cys
                10/4/2011
*/
int MLGP::BiSection(int *nS_, int *nP_1_, int *nP_2_, int alg_, int flag) {
  int ret = GORD_OK, d_ret, i, j;
  MLGP *child = NULL, *parent = NULL, *cur = this, *g_base;
  float s_0, eval, eval_0;
  bool isBetter = false, isPartition = false;
  int inner_refine[10] = {TABU_ONLY_CUT, TABU_ONLY_CUT, TABU_ONLY_CUT, TABU_ONLY_CUT, 0x0,
                          TABU_ONLY_CUT, 0x0,           TABU_ONLY_CUT, 0x0,           TABU_ONLY_CUT};
  ASSERT(TABU_INNER <= 10);

  nBisection++;
  s_0 = nEdge * 1.0f / nNode;
  // down
  while (cur->nNode > THRSH_COARSEST) {
    child = cur->ProjectDown(&d_ret, 0x0);
    nDown++;

    if (d_ret != GORD_OK) {
      ASSERT(child == NULL);
      break;
    }
    //	s = child->nEdge*1.0/child->nNode;
    //	if( s > s_0*1.5 )
    //		break;
    cur = child;
  }
  // Init Partition
  g_base = cur;
  int P_dim = cur->nNode, *P_tag = new int[P_dim * TABU_OUTER], nPass = TABU_OUTER, o_best, i_best;
  int *tag_0 = new int[nNode];
  float *P_eval = new float[TABU_OUTER], obj_0, obj;
  for (i = 0; i < TABU_OUTER; i++) P_eval[i] = FLT_MAX;
  if (MLGP::SEPRATOR == MLGP::SEPARATOR_E) {
    isPartition = cur->Partition_2(&nPass, P_tag, P_eval, 0x0) == GORD_OK;
    ASSERT(isPartition);
  } else
    isPartition = cur->Partition_3(0x0) == GORD_OK;
  if (!isPartition) {
    delete[] P_tag;
    P_tag = NULL;
    delete[] P_eval;
    P_eval = NULL;
    nPass = 1;
  }
  eval_0 = FLT_MAX;  // cur->Objectives( NULL, MLGP::MEASURE,-1,0x0 );
  isBetter = false;
  for (i = 0; i < nPass; i++) {  // up
    //	cur = g_base;
    if (i >= 1) {
      for (j = 0; j < nNode; j++) {
        if (tag[j] == 0) tag[j] = G_RAND() % 2 + 1;  // wPart[2]>wPart[1] ? 1 : 2;
      }
    }
    if (0 && P_eval != NULL) {
      if (P_eval[i] == FLT_MAX) continue;
      memcpy(cur->tag, P_tag + i * P_dim, sizeof(int) * P_dim);
      eval = cur->Objectives(NULL, MLGP::MEASURE, -1, 0x0);
      ASSERT(eval == P_eval[i]);
    }
    while (true) {
      if (!isPartition) {
        if (MLGP::SEPRATOR == MLGP::SEPARATOR_E) {
          isPartition = cur->Partition_2(NULL, NULL, NULL, 0x0) == GORD_OK;
        } else {
          isPartition = cur->Partition_3(0x0) == GORD_OK;
        }
      }
      if (isPartition) {
        obj_0 = FLT_MAX;
        i_best = 0;
        for (j = 0; j < TABU_INNER; j++) {
          if (MLGP::SEPRATOR == MLGP::SEPARATOR_E) {
            cur->Refine_2(MLGP::TABU_REFINE, inner_refine[j]);
            obj = cur->Objectives(NULL, MLGP::MEASURE, -1, 0x0);
            if (obj < obj_0) {
              obj_0 = obj;
              i_best = j;
            }
            if (j > i_best) {
              break;
            }
          } else {
            cur->Refine_3(MLGP::VFM_REFINE, NULL, inner_refine[j]);
          }
        }
      }

      if (cur != this) {
        if (isPartition) cur->ProjectUp(0x0);
        parent = cur->parent;
        parent->CheckPartition(0x0);
        cur = parent;
      } else {
        //			eval = cur->Objectives( NULL,
        // MLGP::MEASURE,-1,0x0
        //);
        if (MLGP::SEPRATOR == MLGP::SEPARATOR_E) {
          Flow_Refine_2(0x0);
          eval = Objectives(NULL, MLGP::MEASURE, -1, 0x0);
        }
        //	for( j = 0; j < TABU_INNER; j++ )	{
        //					Refine_3( MLGP::VFM_REFINE,0x0
        //);
        VFM_Refine(&eval, 0x0);
        if (eval < eval_0) {
          memcpy(tag_0, tag, sizeof(int) * nNode);
          eval_0 = eval;
          isBetter = true;
          o_best = i;
        }
        if (i > o_best + 1) {
          goto OUTER_NEXT;
        }
        //	}
        break;
      }
    }
  }
OUTER_NEXT:
  //	Refine_3( MLGP::VFM_REFINE,0x0 );
  if (isBetter == true) {
    memcpy(tag, tag_0, sizeof(int) * nNode);
    eval = Eval(NULL, 0x0);
    //		ASSERT( eval==eval_0 );
  }

  cur = g_base;
  while (cur != NULL && cur != this) {
    parent = cur->parent;
    delete cur;
    cur = parent;
  }
  if (P_eval != NULL) delete[] P_eval;
  if (P_tag != NULL) delete[] P_tag;
  if (!isPartition) {
    ret = GORD_BISEC_FAIL;
  } else {
  }

  int no, nS = 0, nP_1 = 0, nP_2 = 0;
  for (i = 0; i < nNode; i++) {
    no = cur->tag[i];
    if (no == 0) {
      nS++;  // s_0[nS++] = i;
    } else if (no == 1) {
      nP_1++;  // s_1[nP_1++] = i;
    } else {
      ASSERT(no == 2);
      nP_2++;  // s_2[nP_2++] = i;
    }
  }
  if (nS_ != NULL) *nS_ = nS;
  if (nP_1_ != NULL) *nP_1_ = nP_1;
  if (nP_2_ != NULL) *nP_2_ = nP_2;

  if (nP_1 == 0 || nP_2 == 0) {
    ret = GORD_BISEC_FAIL;
  }

  return ret;
}

extern "C" int ccs_symamd( int nRow,int nCol,int nnz,int *rowind,int* colptr,int* pCol );
/*
        v0.1	cys
                9/27/2011
*/
int MLGP::MdOrder(int algorithm,int flag) {
  clock_t t_0 = G_TIC();
  int set_MIN = 3,set_BIG=MLGP::MD_SWITCH/10,nBigSet=0;    //1
  int ret = GORD_OK, i, j, base = 0, stp = 0;
  int nCut = 0, nLeft = 0, no, cur, *cut = stmp, blk = -1, nSet = 0, nzMax = 0, nz, nMD = 0;
  ASSERT(VERIFY_PERMUTE_VECTOR(nNode, moveto, stmp) == 0);

  int *perm = new int[nNode], *set = new int[nNode], *temp = new int[nNode];
  int B_w = 0, *B_ptr = NULL, *B_adj = NULL;
  for (i = 0; i < nNode; i++) {
    perm[moveto[i]] = i;
    stmp[i] = stp;
  }
  i = 0;
  while (i < nNode) {
    no = perm[i];
    ASSERT(moveto[no] == i);
    blk = kind[no];
    nSet = 0;
    nzMax = 0;
    while (kind[no] == blk && i < nNode) {
      set[nSet++] = no;
      nzMax += ptr[no + 1] - ptr[no];
      if (++i < nNode)
        no = perm[i];
      else
        break;
    }
    ASSERT(blk <= nMDBlock);
    if (blk == 0 && nMDBlock > 1) {  // cut
      base += nSet;
      continue;
    }
    nMD++;
    if (nSet > set_MIN) {
      if(nSet>set_BIG)    nBigSet++;
      stp++;
      B_w = 0, B_ptr = NULL, B_adj = NULL;

      if(algorithm==62){
        nz = GetBlock(nSet, &B_ptr, &B_adj, set, stp,0x0);
        ASSERT(nz <= nzMax);        
        ret = ccs_symamd( nSet,nSet,nz, B_adj,B_ptr, temp );
        for(int k=0;k<nSet;k++){
          B_ptr[temp[k]] = k;
        }
        memcpy(temp,B_ptr,sizeof(int)*nSet);
        //printf( "\tMLGP::62_%d: A=[%d,%d] pCol=[%d,%d,%d]\r\n\n",nMD,nSet,nz,temp[0],temp[nSet/2],temp[nSet-1] );        
      }else{
        nz = GetXBlock_withcut(nSet, &B_w, &B_ptr, &B_adj, set, stp,0x0);
        ASSERT(nz <= nzMax);
        GMDO gmdo(nSet, B_w, B_ptr, B_adj, 0x0);
        ret = gmdo.MD_0(temp, 0x0);
        ASSERT(ret == GORD_OK);
        nLU += gmdo.nLU;
      }
     
      delete[] B_ptr, delete[] B_adj;
      
    } else {
      for (j = 0; j < nSet; j++) temp[j] = j;
      nLU += nzMax;
    }
    for (j = 0; j < nSet; j++) {
      //		no = set[temp[j]];
      //		moveto[no] = ord++;
      moveto[set[j]] = base + temp[j];
    }
    base += nSet;
  }
  ASSERT(base == nNode);
  // ASSERT( nMD==nMDBlock );
  if (0) {
    for (i = 0; i < nNode; i++) {
      cur = perm[i];
      if (kind[cur] == 1) {
        cut[nCut++] = cur;
      } else {
        moveto[cur] = nLeft++;
      }
    }
    for (i = 0; i < nCut; i++) {
      cur = cut[i];
      moveto[cur] = nLeft++;
    }
    ASSERT(nLeft == nNode);
  }
  ASSERT(VERIFY_PERMUTE_VECTOR(nNode, moveto, temp) == 0);
  delete[] perm;
  delete[] set;
  delete[] temp;

  GMDO::t_MD += G_TIC() - t_0;
  printf( "\tMLGP::MdOrder nBigSet=%d(%d) MD_SWITCH=%d T=%.3g\r\n",nBigSet,nMD,MD_SWITCH,GMDO::t_MD );    
  return ret;
}

int CCS2AAT_struc(int dim, int nnz, int *a_ptr, int *a_ind, int *ptr, int *ind, int *temp);
/*
        Init from CCS formant of MATRIX
        v0.1	cys
                9/27/2011

2020/10/10	cys
	1) Greedy minimum degree order is important and need careful redesign. Maybe rewrite syam with high efficiency code
	2) BiSection_2!!!	It's really hard to optimize and the nzLU is so sensitive to this function

*/
extern "C" int ccs_MLGP_(int nCol, int *ptr_0, int *ind_0, int *weight, int *perm,int algorithm, int flag) {
  int ret = GORD_UNKNOWN_ERROR, nCom = 0, i, top = 0, base = 0, nnz = ptr_0[nCol], nzAAT, C_dim;
  int *ptr = NULL, *ind = NULL;
  int sepa = MLGP::SEPRATOR, *comp = new int[nCol], *temp = new int[nCol];
  float e = 0.0;

  MLGP::InitParam(0x0);
  MLGP::MD_SWITCH = 10000;
  // MLGP::THRSH_COARSEST = MLGP::MD_SWITCH;
  if (MLGP::SEPRATOR == MLGP::SEPARATOR_V) {
    MLGP::TABU_OUTER = 1, MLGP::TABU_INNER = 1;
  }

  double t_0 = G_TIC(), t_r;
  g_tX = 0;
  GMDO::ClearStat(0x0);

  bool isComp = false, isExpand = false;
  MLGP *graph_0 = NULL;
  if (!BIT_TEST(flag, GORD_MAT_SYMMETRY)) {
    int *deg = comp;
    nzAAT = CCS2AAT_len(nCol, ptr_0, ind_0, deg, temp);
    if (nzAAT > nnz) {
      ptr = new int[nCol + 1];
      ind = new int[nzAAT];
      ptr[nCol] = nzAAT;
      for (i = nCol - 1; i >= 0; i--) ptr[i] = ptr[i + 1] - (deg[i] + 1);
      ASSERT(ptr[0] == 0);
      CCS2AAT_struc(nCol, nnz, ptr_0, ind_0, ptr, ind, temp);
      isExpand = true;
    }
  }
  if (isExpand) {
  } else {
    ptr = ptr_0;
    ind = ind_0;
  }
  C_dim = CCS_Compress_(nCol, ptr, ind, weight, comp, 0x0, 0x0);
  isComp = C_dim < nCol;
  
  if (isComp) {
    graph_0 = new MLGP(nCol, ptr, ind, weight, comp, MLGP_CONNECTED, sepa, flag);
    if (MLGP::DUMP > 0) {
      printf("CCS_Compress_(%d-%d)=>(%d,%d)\r\n", nCol, ptr[nCol], graph_0->nNode, graph_0->nEdge);
    }
  } else {
    graph_0 = new MLGP(nCol, ptr, ind, weight, MLGP_CONNECTED, sepa, flag);
  }
  ASSERT(graph_0->isSymmetry(0x0));

  if (isExpand) {
    delete[] ptr;
    delete[] ind;
  }
  t_0 = G_TIC();
  graph_0->NestSectionOrder(NULL, -1, 0x0, 1, 0x0);
  t_r = G_TIC() - t_0;
  graph_0->MdOrder(algorithm,0x0);
  MLGP::DumpParam(0x0);

  int *moveto = graph_0->moveto, nNode = graph_0->nNode, nCut = 0, nLeft = 0, no, cur, *cut = graph_0->stmp;
  ASSERT(VERIFY_PERMUTE_VECTOR(nNode, moveto, NULL) == 0);
  if (isComp) {
    int *next = new int[nCol], *head = new int[nCol], i, nz = 0;
    for (i = 0; i < nCol; i++) {
      next[i] = -1;
      head[i] = -1;
    }
    for (i = 0; i < nCol; i++) {
      no = comp[i];
      if (no == -1) continue;
      ASSERT(no >= 0 && no < nCol);
      if (head[no] == -1) {
        head[no] = i;
      } else {
        cur = head[no];
        next[i] = cur;
        head[no] = i;
      }
    }
    nz = 0;
    for (i = 0; i < nNode; i++) {
      perm[moveto[i]] = i;
    }
    for (i = 0; i < nNode; i++) {
      no = perm[i];
      cur = head[no];
      while (cur != -1) {
        temp[cur] = nz++;
        cur = next[cur];
      }
    }
    for (i = 0; i < nCol; i++) {
      perm[temp[i]] = i;
    }
    delete[] next;
    delete[] head;
  } else {
    ASSERT(nCol == nNode);
    for (i = 0; i < nCol; i++) perm[moveto[i]] = i;
  }
  GMDO::DumpStat(0x0);
  printf("\tccs_MLGP_:Algorithm=%d GreedyCutEdge=%d\tnLU=%d,time=(%.3g,NS=%.3g)\r\n", algorithm,(int)(graph_0->isGreedyCutEdge),graph_0->nLU, G_TIC() - t_0, t_r);
  printf(
      "\tMLGP::tPartition=%.3g,tBiSection=%.3g(%d),tSub=%.3g,tSplit=%."
      "3g\n\ttRefine=%.3g,%.3g,\ttEval=%.3g(%d)\tMLGP::tX=%.3g\n\n",
      MLGP::tPartition, MLGP::tBiSection, MLGP::nBisection, MLGP::tSub, MLGP::tSplit, MLGP::tFlowRefine, MLGP::tRefine,
      MLGP::tEval, MLGP::nEval, MLGP::tX);
  delete graph_0;

  ASSERT(VERIFY_PERMUTE_VECTOR(nCol, perm, NULL) == 0);
  ret = GORD_OK;

  delete[] temp;
  delete[] comp;
#ifdef _DEBUG
//	_CrtDumpMemoryLeaks( );
#endif
  return ret;
}

/*
        v0.1	cys
                9/23/2020
*/
int qmd_order(int nRow, int nCol, int nnz, int *rowind, int *colptr, int *pCol) {
  int flag = 0x0, iRet = 0x0;
  // size_t nnz = colptr[nCol];
  MLGP mlgp(nCol, colptr, rowind, NULL, 0x0, 0x0, flag);
  iRet = mlgp.MdOrder(6,flag);
  return iRet;
}
