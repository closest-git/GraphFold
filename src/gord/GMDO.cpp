#include <memory.h> 
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "gord.h"
#include "GORD_util.h"
#include "GMDO.h"
#include "MLGP.h"
#include "../util/grus_util.h"

/*
	1 node����LIST��˳��ֱ��Ӱ�������������selectcand?
	2 Ϊ��DEGREE�Ĳ�ͬѡ�������ô��Ĳ���?
	3 SelectCand ��������Ƿ���Ч��
	4 SelectCand Ϊ��outlay��Ч��

*/
float GMDO::t_MD=0,GMDO::t_X=0;
int GMDO::nMD=0;
int GMDO::CAND=0;
int GMDO::MD_LU=0;
float GMDO::rZERO=0;
float GMDO::rSuper=0;
float GMDO::rDim=0;
int GMDO_NODE::OUTLAY_CEIL=0;

GMDO_NODE::~GMDO_NODE( )	{	
	if( tie!=NULL )		
		::free( tie );	
}
/*
	v0.1	cys
		5/29/2012
*/
int GMDO_NODE::Init( int no,int len,int *a,int flag )	{
	int i=0,ret = GORD_OK;
//	id = no;
	nTie = len+1;		
	nSize = MAX( nTie*2,64 );
	tie = (int*)::malloc( sizeof(int)*nSize );
	tie[0] = no;		super=1;		//basic supervariable 
	key = 0;	
	for( i = 0; i < len; i++ )	{
		tie[i+1] = a[i];
		key += tie[i+1];
		ASSERT( tie[i+1]!=no );
	}
	degree = len;				extra = 0;
	outlay = nTie-1;
	ord = -1;
	rank_0 = -1;

	return ret;
}

/*
	v0.1	cys
		6/1/2012

float GMDO_NODE::Degree( )	{	
//	return degree;	
	return degree+super;	
//	return nTie;	
//	return nTie-super;	
}
*/
/*
	v0.1	cys
		5/30/2012
*/
int GMDO_NODE::RankList( SORT_BUCKET* hList,int flag )	{	
//	int rank = (int)degree;
	int rank = (int)(degree+super);	
//	int rank = nTie;	
//	int rank = nTie-super;	
	if( rank>OUTLAY_CEIL )
		rank = OUTLAY_CEIL;
	if( rank!=rank_0 )	{
		int id = tie[0];
		if( rank_0>=0 )	{
			hList->Remove( rank_0,id,0x0 );
		}
		if( rank>=0 )	{
			hList->Insert( rank,id,0x0 );
		}

		rank_0 = rank;
	}
	return rank;	
}
/*
	v0.1	cys
		5/29/2012
*/
int GMDO_NODE::Eliminate( GMDO_NODE *hPivot,int *stmp,int stp,int flag )	{
	int i=0,no,nz=0,nExtra=0,id=tie[0];

	key = 0;			extra=0;
	for( i = 0; i < nTie; i++ )	{
		no = tie[i];		
		if( stmp[no]==-stp )	{	//supervariable of hPivot
		}else	{
			tie[nz++] = no;			
			key+= no;			
			stmp[no]*=-1;
		}
	}

	stmp[id] = 0;
	for( i = hPivot->super; i < hPivot->nTie; i++ )	{
		no = hPivot->tie[i];	
		if( stmp[no]>0 )	{
			nExtra++;
		}
	}
	if( nz+nExtra>nSize )	{
		int nSize_1 = MAX( nSize*2,nz+nExtra*2 );
		tie = (int*)::realloc( tie,nSize_1*sizeof(int) );
		nSize = nSize_1;
	}
	for( i = hPivot->super; i < hPivot->nTie; i++ )	{
		no = hPivot->tie[i];	//ASSERT( no>=0 );
		if( stmp[no]<=0 )	{
			stmp[no] *= -1;
		}else	{
			tie[nz++] = no;			key += no;
		}
	}
	stmp[id]=stp;

	nTie = nz;
	ASSERT( nTie<=nSize );

	if( 0 )	{		//�����ظ�
		for( i = 0; i < nTie; i++ )	{
			no = tie[i];	
			for( int j = i+1; j < nTie; j++ )
				ASSERT( tie[j]!=no );			
		}
	}

	return nExtra;
}

/*
	v0.1	cys
		5/29/2012
*/
int GMDO_NODE::SuperAbsorb( int nNode,GMDO_NODE *nodes,int pivot,int flag )	{
	int nz = 0,i,no,nAborb=0;
	ASSERT( super>0 );
	for( i = 0; i < super; i++ )	{
		no = tie[i];
		nz++;
		ASSERT( nodes[no].super>0 || -nodes[no].super==tie[0] );
	}
	for( i = super; i < nTie; i++ )	{
		no = tie[i];
		if( no<nNode && nodes[no].super==-pivot )
			continue;
		tie[nz++] = no;
	}
	nAborb = nTie-nz;
	nTie = nz;

	return nAborb;
}

/*
	v0.1	cys
		5/30/2012
*/
int GMDO_NODE::DetectMatch( GMDO_NODE *hNext,int *stmp,int stp,int flag )	{
	int i=0,no,isMatch=1,pos=-1,id=tie[0];
	ASSERT( nTie==hNext->nTie && hNext->super>0 );
	for( i = 0; i < nTie; i++ )	{
		no = hNext->tie[i];
		if( stmp[no]!=stp )
		{	isMatch=0;	break;		}
		stmp[no]*=-1;
	}
	for( i = 0; i < nTie; i++ )	{
		no = tie[i];	
		if( stmp[no]<=0 )	{
			stmp[no] *= -1;
		}else	{
			isMatch=0;	
		}
	}
	if( isMatch == 1	)	{
		for( i = 0; i < super; i++ )	{
			no = tie[i];
			stmp[no]*=-1;
		}
		for( i = 0; i < hNext->super; i++ )	{
			no = hNext->tie[i];
			if( stmp[no]<0 )	{
				ASSERT( 0 );
			}else	{
				tie[super++] = no;
				stmp[no]*=-1;
			}
		}
		int nz = super;
		for( i = hNext->super; i < nTie; i++ )	{
			no = hNext->tie[i];
			if( stmp[no]<=0 )	{
				continue;
			}else	{
				tie[nz++] = no;
			}
		}
		for( i = 0; i < super; i++ )	{
			no = tie[i];
			stmp[no]*=-1;
		}
		ASSERT( nz==nTie );
		hNext->super = -id;
	}

	return isMatch;
}

/*
	v0.1	cys
		5/29/2012
*/
GMDO::GMDO(int h,int w,int *p,int *a,int flag)	{
	GMDO_NODE *hNode = NULL;
	nodes = NULL;
	int i;

	nFix = 0;
	height = h;				width = w;
	ASSERT( width>=height );
	stmp = new int[width];		//elim = new int[dim];
	hash = new int[width];
	temp = new int[width];
	nodes = new	GMDO_NODE[height];
	for( i = 0; i < height; i++ )	{
		hNode = nodes+i;
		hNode->Init( i,p[i+1]-p[i],a+p[i],0x0 );
	}
	for( i = height-nFix; i < height; i++ )	{
		hNode = nodes+i;
		hNode->super=-1;
	}
}
GMDO::~GMDO( )	{
	if( stmp!=NULL )	delete[] stmp;
	//delete[] elim;
	if( nodes!=NULL )	delete[] nodes;
	if( hash!=NULL )	delete[] hash;
	if( temp!=NULL )	delete[] temp;
}

/*
	v0.1	cys
		5/29/2012
*/
int GMDO::SelectCand( SORT_BUCKET *list,int *modify,int flag )	{
	int cand=-1,flr,nM,j;
	float w1,w2,w1_best,w2_best,w3_best;
	GMDO_NODE *hNode = NULL,*hAdj=NULL;

	flr = list->Floor_Bin( );					ASSERT( flr>=0 );
	nM = list->Bin_List( flr,modify,0x0 );		ASSERT( nM>0 );
	if( CAND==0 )	{
		cand = modify[0];
		return cand;
	}
	if( nM>1 )	{
		cand = modify[rand( )%nM];
		w1_best = FLT_MAX,				w2_best=FLT_MAX;			w3_best=FLT_MAX;
		for( j = 0; j < nM; j++ )	{
			hNode = nodes+modify[j];
//			w1 = hNode->outlay;			w2=hNode->Degree( );
			ASSERT( hNode->rank_0==flr );
			w1 = hNode->rank_0;		//hNode->Degree( );	
			if( hNode->outlay>=0.0 )
				w2 = hNode->outlay;
			else
				w2 = UpdateOutlay( modify[j],w2_best,UPDATE_OUTLAY );		
//			w3 = stp-hNode->history;
			if( w1<w1_best )	{
				w1_best = w1;			w2_best=w2;		//deg_1 = hNode->nTie;
				cand = modify[j];
			}else if( w1==w1_best )	{
				if( w2<w2_best )	{
					w2_best = w2;		/*w3_best = w3;*/
					cand = modify[j];
				}/*else if( w2==w2_best )	{
					if( w3<w3_best )	{
						w3_best = w3;		
						cand = modify[j];
					}
				}*/
			}
		}		
	}else
		cand = modify[0];

	return cand;
}

/*
	v0.1	cys
		5/29/2012
*/
int GMDO::MD_0( int *moveto,int flag )	{
	int ret = GORD_OK,i,j,k,cand,no,flr,nM=0,nExtra,nElim=0,nSuper=0;
	int nZero=0,stp_0,nNeedElim=height-nFix,deg_0,deg_base;
	int *modify = new int[width];
	clock_t t_1;
	SORT_BUCKET list;
	GMDO_NODE *hNode = NULL,*hAdj;

	for( i = 0; i < width; i++ )	stmp[i]=-1;
	stp=0;
	GMDO_NODE::OUTLAY_CEIL = width+1;
	list.Init( GMDO_NODE::OUTLAY_CEIL+1,height,0x0 );
	for( i = 0; i < height; i++ )	{
		hNode = nodes+i;
		if( hNode->isFix( ) )
			continue;
		UpdateOutlay( i,FLT_MAX,UPDATE_OUTLAY );
		hNode->RankList( &list,0x0 );
		stmp[i] = stp;			hash[i]=stp;	
	}
	nElim=0;					nLU=0;
	while( nElim<nNeedElim )	{
		stp++;		nSuper++;	
		flr = list.Floor_Bin( );		
		if( flr==-1 )	{
			for( j = 0; j < height; j++ )	{
				hNode = nodes+j;
				ASSERT( hNode->ord>=0 );
			}			
		}
		if( flr==0 )		nZero++;
		ASSERT( flr>=0 );
		t_1=clock( );
		cand = SelectCand( &list,modify,0x0 );			ASSERT( cand!=-1 );	
		t_X += clock( )-t_1;
		hNode = nodes+cand;
		list.Remove( flr,cand,0x0 );
//		ASSERT( isSymmetry( 0x0 ) );
		stp++;
		deg_base = 0;
		for( j = 0; j < hNode->super; j++ )			{
			no = hNode->tie[j];				ASSERT( nodes[no].ord==-1 );
			ASSERT( -nodes[no].super==cand || nodes[no].super>0  );		
			deg_base += nodes[no].super;
			nodes[no].ord = nElim++;				
			stmp[no] = -stp;
		}
		nLU += (2*hNode->degree+hNode->super)*hNode->super;
		for( j = hNode->super; j < hNode->nTie; j++ )			{	
			no = hNode->tie[j];		
			if( !isExtraBorder(no) )
			{	ASSERT( nodes[no].super>0 && nodes[no].ord==-1);}
			stmp[no] = stp;
		}
		for( j = hNode->super; j < hNode->nTie; j++ )	{
			no = hNode->tie[j];		
	//		if( cand==54  )
	//			cand = 54;
			if( isExtraBorder(no) )
				continue;
			hAdj = nodes+no;		ASSERT( hAdj->super>=0 );
	//		rank_0 = hAdj->Rank( 0x0 );
			nExtra = hAdj->Eliminate( hNode,stmp,stp,0x0 );
			hAdj->key = (hAdj->key)%width;	
		//	if( hNode->outlay==0.0 )
		//		ASSERT( nExtra==0 );
		//	ASSERT( hAdj->nAdj>0 );
		/*	if( hAdj->nAdj==0 )	{	//absorb
				elim[no]=1;
				hNode = nodes+no;
				hNode->ord = ++i;
			}else*/	
		}
		nM = 0;			stp_0 = stp;
		SuperDetect( &list,cand,0x0 );
		for( j = hNode->super; j < hNode->nTie; j++ )	{
			hAdj = nodes+hNode->tie[j];
			if( isExtraBorder( hNode->tie[j] ) )
				continue;
			if( hAdj->super>0 )	{
				modify[nM++] = hNode->tie[j];
				deg_0=hAdj->degree;			hAdj->degree = 0;
				for( k = hAdj->super;k < hAdj->nTie; k++ )	{
					no = hAdj->tie[k];
					hAdj->degree+=isExtraBorder(no) ? 1 : nodes[no].super;
				}
			//	ASSERT( hAdj->degree+deg_base-deg_0>=0 );
			//	hAdj->extra += hAdj->degree+deg_base-deg_0;
			}else	{
				list.Remove( hAdj->rank_0,hNode->tie[j],0x0 );
			}
		/*	
			for( k = 0; k < hAdj->nAdj; k++ )	{
				no = hAdj->tie[k];
				if( stmp[no]==stp )		
					continue;
				modify[nM++] = no;
				stmp[no]=stp;
			}*/
		}
		for( j = 0; j < nM; j++ )	{	
			nodes[modify[j]].outlay = -1.0;
		//	UpdateOutlay( modify[j],0x0 );
			nodes[modify[j]].RankList( &list,0x0 );
		}
	}

	for( i = 0; i < height; i++ )	{
		no = nodes[i].ord;
		moveto[i] = no;
		if( i>=nElim )
		{	ASSERT( nodes[i].isFix( ) );}
	//	moveto[no] = i;
	}
	ASSERT( VERIFY_PERMUTE_VECTOR( height,moveto,temp ) == 0 );	
	delete[] modify;		

	MD_LU += nLU;
	rDim += height;
	rSuper += nSuper*1.0/height;
	rZERO += nZero*1.0/height;
	
	nMD++;

	return ret;		
}

/*
	v0.1	cys
		5/30/2012
*/
int GMDO::SuperDetect( SORT_BUCKET *list,int pivot,int flag )	{
	GMDO_NODE *hPivot=nodes+pivot,*hNode=NULL,*hAdj=NULL;
	int i,j,nz,no,adj;
	int stp_0=stp,nSuper=0,nMatch=0;

	for( i = hPivot->super; i < hPivot->nTie; i++ )	{
		no = hPivot->tie[i];	
		if( isExtraBorder(no) )
			continue;
		hNode = nodes+no;
		if( hNode->super<=0 )
			continue;
		if( hash[hNode->key]!=stp_0 )	{	//supervariable detection
			hash[hNode->key]=stp_0;			continue;
		}
		stp++;
		for( j = 0; j < hNode->nTie; j++ )	{
			stmp[hNode->tie[j]] = stp;
		}
		nMatch = 0;
//		for( j = hNode->super; j < hNode->nTie; j++ )	{
//			adj = hNode->tie[j];
		for( j = i+1; j < hPivot->nTie; j++ )	{
			adj = hPivot->tie[j];
			if( isExtraBorder(adj) )
				continue;
			hAdj = nodes+adj;
			if( hAdj->super<=0 )
				continue;
			if( hNode->key==hAdj->key && hNode->nTie==hAdj->nTie )	{
				if( hNode->DetectMatch( hAdj,stmp,stp,0x0 ) )	{
					nMatch++;
				}
			}
		};
		if( nMatch>0 )	
			temp[nSuper++]=no;
	}
	
	for( i = 0; i < nSuper; i++ )	{
		hNode = nodes+temp[i];
		ASSERT( hNode->super>0 );
		for( j = hNode->super; j < hNode->nTie; j++ )	{
			adj = hNode->tie[j];
			if( isExtraBorder(adj) )				
				continue;
			hAdj = nodes+adj;
			if( hAdj->super<=0 )
				continue;
			nz = hAdj->SuperAbsorb( height,nodes,temp[i],0x0 );
			ASSERT( nz>0 );		
		}
	}
	return 0x0;
}
/*
	v0.1	cys
		5/30/2012
*/
int GMDO::DumpStat( int flag )	{
	printf( "\tGMDO(%d,%g):CAND=%d,SUPER=%f,ZERO=%f,Time=%d(%d),LU=%d\r\n",nMD,rDim/nMD,CAND,rSuper/nMD,rZERO/nMD,t_MD,t_X,MD_LU );

	return 0x0;
}
int GMDO::ClearStat( int flag )	{
	MD_LU =0;
	nMD = 0;			t_MD=0;					t_X=0;
	rZERO=0.0,			rSuper=0.0;				rDim=0.0;

	return 0x0;
}

/*
	v0.1	cys
		5/30/2012
*/
float GMDO::UpdateOutlay( int var,float w_best,int flag )	{
	GMDO_NODE *hNode=nodes+var,*hAdj=NULL;
	ASSERT( hNode->super>0 );
	int no,i,j,no_j;
	float w = 0,w_0=0,a;		

//	rank_0 = hNode->Rank( 0x0 );
	stp++;
	
	for( i = hNode->super; i < hNode->nTie; i++ )	{
		no = hNode->tie[i];		
		temp[no] = no>=height? 1 : nodes[no].super;			
		ASSERT( temp[no]>0 );
		w += temp[no];
		stmp[no] = stp;	//i < hNode->super ? -stp : stp;
	}
	if( !BIT_TEST( flag,UPDATE_OUTLAY ) )	{
		hNode->outlay = -1.0;		goto END;
	}
	ASSERT( hNode->degree==w );

	if( 1 )	{
		w_0 = w*w;		w=w_0;
		//stmp[var] = -stp;
		for( i = hNode->super; i < hNode->nTie; i++ )	{
			no = hNode->tie[i];
			if( isExtraBorder(no) )
				continue;
			hAdj = nodes+no;
			a=0;
			for( j = 0; j < hAdj->nTie; j++ )	{
				no_j = hAdj->tie[j];
				if( stmp[no_j]==stp )
					a += temp[no_j];
			}
			w -= hAdj->super*a;
	//		w += hAdj->GetExtra( hNode,stmp,stp,0x0 );	
		}
		ASSERT( w>=0 );		
	}else	{
		w = 0;
		for( i = hNode->super; i < hNode->nTie; i++ )	{
			no = hNode->tie[i];
			hAdj = nodes+no;
			a = (hAdj->degree+hNode->degree)/2;
			w += a;
		}
	}
	hNode->outlay = w;

END:
	return w;
}

int GMDO::isSymmetry( int flag )	{
	int isSymm=1,i,j,k,no,isMatch;
	GMDO_NODE *hNode=NULL,*hAdj=NULL;
	for( i = 0; i < height; i++ )	{
		hNode = nodes+i;
		if( hNode->super<=0 )
			continue;
		if( hNode->ord>=0 )
			continue;
		for( j = hNode->super; j < hNode->nTie; j++ )	{
			no = hNode->tie[j];
			if( no>=height )
				continue;
			hAdj = nodes+no;
			isMatch = 0;
			for( k = hAdj->super; k < hAdj->nTie; k++ )	{
				if( hAdj->tie[k]==i )	{
					isMatch=1;		break;
				}
			}
			ASSERT( isMatch==1 );
		}
	}
	return isSymm;
}