#pragma once

class SORT_BUCKET;
class GMDO_NODE{
public:
	enum{
		FIX=0x1000,
	};
//	G_INT_64 id,
	G_INT_64 nTie,super,key,ord,flag,rank_0,extra;
	float outlay,degree;
	static G_INT_64 OUTLAY_CEIL;
	G_INT_64 *tie,nSize;

	float Degree( );
	bool isFix( )	{	return BIT_TEST(flag,FIX);		}
	G_INT_64 Init( G_INT_64 no,G_INT_64 len,G_INT_64 *a,G_INT_64 flag );
	G_INT_64 RankList( SORT_BUCKET* hList,G_INT_64 flag );
//	G_INT_64 GetExtra( GMDO_NODE *hPivot,G_INT_64 *stmp,G_INT_64 stp,G_INT_64 flag );
	G_INT_64 SuperAbsorb( G_INT_64 nNode,GMDO_NODE *nodes,G_INT_64 pivot,G_INT_64 flag );
	G_INT_64 Eliminate( GMDO_NODE *hPivot,G_INT_64 *stmp,G_INT_64 stp,G_INT_64 flag );
	G_INT_64 DetectMatch( GMDO_NODE *hNext,G_INT_64 *stmp,G_INT_64 stp,G_INT_64 flag );

	~GMDO_NODE( );
	GMDO_NODE( ){	memset( this,0x0,sizeof(GMDO_NODE) );	}	
};
/*
*/	

class GMDO	{
	G_INT_64 height,width,*temp,*hash,*stmp,stp,nFix;
//	G_INT_64 *ptr,*adj,*elim;
	GMDO_NODE *nodes;	

	G_INT_64 SelectCand( SORT_BUCKET *list,G_INT_64 *modify,G_INT_64 flag );
	float UpdateOutlay( G_INT_64 no,float best,G_INT_64 flag );
//	float NodeOutlay( G_INT_64 no,G_INT_64 flag );
	G_INT_64 SuperDetect( SORT_BUCKET *list,G_INT_64 pivot,G_INT_64 flag );
	G_INT_64 isSymmetry( G_INT_64 flag );
	bool isExtraBorder( G_INT_64 no )	{	ASSERT( no>=0&& no<width);	return no>=height; }
public:
	G_INT_64 nLU;
	enum{
		UPDATE_OUTLAY=0x10000,
	};
	static G_INT_64 nMD,CAND,MD_LU;
	static float rZERO,rSuper,rDim,t_MD,t_X;

	GMDO(G_INT_64 h,G_INT_64 w,G_INT_64 *ptr,G_INT_64 *adj,G_INT_64 flag);
	~GMDO( );

	G_INT_64 MD_0( G_INT_64 *permut,G_INT_64 flag );
	static G_INT_64 DumpStat( G_INT_64 flag );
	static G_INT_64 ClearStat( G_INT_64 flag );
};