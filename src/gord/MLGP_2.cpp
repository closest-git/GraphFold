#include <memory.h> 
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "gord.h"
#include "GORD_util.h"
#include "MLGP.h"
#include "Flow.h"
#include "GMDO.h"
#include "../util/grus_util.h"
extern clock_t g_tX;

int MLGP::DumpParam( int flag )	{
	printf( "\tMLGP:THRSH=(%d,%g),CAND=(%d),MD=%d(%d),DOWN=%d,MEASURE=%d,PARTI=%d,TABU=(%d,%d,%g),SEPRATOR=%d\t",THRSH_COARSEST,
		T_COARSEST_RELAX,CAND_TABU,MDO,MD_SWITCH,DOWN,MEASURE,PARTI,
		TABU_OUTER,TABU_INNER,TABU_rMove,SEPRATOR );
	printf( "\r\n\tSTAR=(%g,%g),Bisec=%d,Down=%d,tX=%d,rMove=%g(%g),rEval=%g(%g)\r\n",STAR,STAR_0,
		nBisection,nDown,g_tX,rMove,(rMove/nDown),rEval,(rEval/nDown) );
	
	return 0x0;
}
int MLGP::InitParam( int flag )	{
	MLGP::STAR_0=-1.0;;

	MLGP::nBisection=0;			MLGP::nVFM=0;				MLGP::nTaRf=0;
	MLGP::nMDBlock = 0;
	MLGP::nDown=0;
	MLGP::rMove=0;
	MLGP::rEval=0;

	MLGP::tPartition=0,MLGP::tBiSection=0,MLGP::tSub=0,MLGP::tSplit=0;
	MLGP::tX = 0;		MLGP::tFlowRefine=0;		MLGP::tRefine=0;
	MLGP::tEval=0;

	return 0x0;
}
/*
	Weighted Graph Cuts without Eigenvectors

	v0.1	cys
		5/9/2012
*/
int Kernel_Kmeans( int nC,int nP,int *cluster,float *K,float *W,int flag )	{
	int ret=GORD_OK,l=0,LOOP_MAX=10,i,j,c_0,c,nMove,nMove_1=0;
	float *K_row=NULL,w=1.0,d_0,*D,*Vol,*P_w,D_sum=0.0;

	D=new float[nC];				Vol=new float[nC];
	P_w=new float[nP];		
	do	{
		D_sum = 0.0;
		nMove = 0;
		for( j = 0; j < nC; j++ )	{
			Vol[j]=0.0;			D[j]=0.0;
		}
		for( i = 0; i < nP; i++ )	{
			c_0 = cluster[i];		ASSERT( c_0>=0 && c_0<nC );
			D[c_0] += W[i];
		}
		for( i = 0; i < nP; i++ )	{
			c_0 = cluster[i];		
			P_w[i] = W[i]/D[c_0];
		}
		for( i = 0; i < nP; i++ )	{
			c_0 = cluster[i];		ASSERT( c_0>=0 && c_0<nC );
			K_row = K+i*nP;
			for( j = 0; j < nP; j++ )	{
				if( c_0 == cluster[j] )	{
					Vol[c_0] += P_w[i]*P_w[j]*K_row[j];
				}				
			}
		}
		for( i = 0; i < nP; i++ )	{
			K_row = K+i*nP;
			for( j = 0; j < nC; j++ )	{
				D[j]=0.0;	
			}
			for( j = 0; j < nP; j++ )	{
				if( K_row[j]==0.0 )
					continue;
				c = cluster[j];
				D[c] += K_row[j]*P_w[j];
			}
			d_0 = FLT_MAX;	c_0=-1;	
			for( j = 0; j < nC; j++ )	{
				D[j] = K_row[i]-2.0*D[j]+Vol[j];		ASSERT( D[j]>0.0 );
				if( D[j]<d_0 )	
				{	d_0 = D[j];		c_0 = j;	}
			}
			D_sum += D[cluster[i]];
			if( cluster[i]!=c_0 )	{
				cluster[i] = c_0;
				nMove++;
			}
		}
		nMove_1 += nMove;
	}while( l++<LOOP_MAX && nMove!=0 );
	delete[] D;				delete[] Vol;
	delete[] P_w;			

	return nMove_1;
}
/*
	Balanced Separator

	v0.1	cys
		5/9/2012
*/
int MLGP::Partition_2( int *nP,int *P_tags,float *P_eval,int flag ){
	int ret=GORD_VSEP_FAIL;
	switch( MLGP::PARTI )	{
	case 1:
		ret = P2_RG( flag );
		break;
	case 2:
		ret = P2_GreedyGrow( nP,P_tags,P_eval,flag );
		break;
	default:
		ASSERT( 0 );
		break;
	}
	return ret;
}

/*
	Region Growing
	v0.1	cys
		5/9/2012
*/
int MLGP::P2_RG( int flag ){
	if( nNode<2 )
		return GORD_TOO_SMALL;
	int ret=GORD_VSEP_FAIL;
	int i=0,j,no,*tag_0,T_1,seed,cur,nRand,l,*stack,top;
	float eval,eval_0=FLT_MAX,s_0,s_1=0.0;

	stack=new int[nNode];		tag_0=new int[nNode];
	s_0 = (float)(nEdge*100.0/(nNode*(nNode-1)));
	nRand = MIN( 200,nNode/2 );
	//	srand( (unsigned)time( NULL ) );
	for( l = 0; l < nRand; l++ )	{
		seed = G_RAND( )%nNode;
	//	seed = list[l];
		T_1 = 0;
		for( i = 0; i < nNode; i++ )		tag[i]=2;
		top=0;		stack[top++]=seed;
		tag[seed] = 1;
		while( top>0 )	{
			cur = stack[--top];
			for( j = ptr[cur]; j < ptr[cur+1]; j++ )	{
				no = adj[j];
				if( tag[no]==2)	{
					tag[no]=1;		T_1++;
					stack[top++] = no;
				}
				if( T_1>=nNode/2 )
					goto NEXT;
			}
		}
NEXT:		
		if( T_1==0 || nNode-T_1==0 )
			continue;
		eval = Objectives( tag,MLGP::MEASURE,-1,0x0 );
		if( eval<eval_0 )	{
			eval_0 = eval;		memcpy( tag_0,tag,sizeof(int)*nNode );			
			ret = GORD_OK;
		//	Kernel_Refine( 0x0 );
		}
	}
	if( ret==GORD_OK )	{
		memcpy( tag,tag_0,sizeof(int)*nNode );	
		eval = Objectives( NULL,MLGP::MEASURE,-1,0x0 );		ASSERT( eval==eval_0 );
#ifdef _DEBUG
	//	printf( "====Partition G=(%d,%d),(%d,%d,%g),r=%g%%\r\n",nNode,nEdge,T_1,nNode-T_1,eval,nEdge*100.0/(nNode*(nNode-1)) );
#endif		
	}
	delete[] tag_0;			delete[] stack;

	return ret;
}

/*
	v0.1	cys
		6/7/2012
*/
int MLGP::UpdateBestPass( float eval,int *nP,int *P_tags,float *P_eval,int flag )	{
	int ret=GORD_OK,i,pos=-1;
	float e_fit=FLT_MAX;
	for( i = 0; i <*nP; i++ )	{
		if( eval<P_eval[i] )	{
			if( pos==-1 )	{
				pos = i;	e_fit = P_eval[i];	
			}else	{
				if( P_eval[i]>e_fit )	{
					pos = i;	e_fit = P_eval[i];	
				}
			}
		}
	}
	if( pos!=-1 )	{
		P_eval[pos] = eval;
		//memcpy( tag_0,tag,sizeof(int)*nNode );	
		memcpy( P_tags+nNode*pos,tag,sizeof(int)*nNode );
	}
	return ret;
}
/*
	Greedy Region Growing
	v0.1	cys
		5/21/2012
*/
int MLGP::P2_GreedyGrow( int *nP,int *P_tags,float *P_eval,int flag ){
	if( nNode<2 )
		return GORD_TOO_SMALL;
	ASSERT( isSymmetry(0x0)==1 );
	*nP = TABU_OUTER;
	int ret=GORD_VSEP_FAIL,isBetter=0,T_pass;
	int i=0,j,no,*tag_0,T_1,seed,nRand,l,*frontier=NULL,*update=NULL,base=0,nUpdate=0,nM=0,cand=-1,stp=0;
	float eval,eval_0=FLT_MAX,eval_l,w_max=0.0;

	frontier=new int[nNode];		update=new int[nNode];
	tag_0=new int[nNode];
	nRand = MIN( 10,nNode/2 );
//	nRand = nNode;
	T_pass = nNode*0.7;
	SORT_LIST list;
	base = (int)(gain_limit/gain_unit_edge+1);
	list.Init( nNode,0x0 );
	for( l = 0; l < nRand; l++ )	{
		isBetter=0;
		list.Clear( );
		for( i = 0; i < nNode; i++ )	{
			tag[i]=2;		gain[i] = 0;		stmp[i]=stp;
		}
		eval = Objectives( tag,MLGP::MEASURE,-1,0x0 );			eval_l=eval;

		seed = G_RAND( )%nNode;
	//	seed = l;
		UpdateGain_E( seed,NULL,(float)base,NULL,0x0 );
		list.UpdateWeight( seed,gain[seed],0x0 );
		T_1 = 0;
		while( T_1< T_pass )	{
			nUpdate = 0;			stp++;
			nM = list.First( frontier,0x0 );			ASSERT( nM>0 );
		//	for( i = 0; i < nM; i++ )	{
		//		cand = frontier[i];
				cand = nM>1 ? frontier[G_RAND( )%nM] : frontier[0];
				ASSERT( ptr[cand+1]-ptr[cand]>0 );
				list.Remove( cand );
				eval = Objectives( NULL,MLGP::MEASURE,cand,UPDATE_WPART );
				tag[cand]=1;		T_1++;
				if( eval<eval_l )		eval_l = eval;
				if( eval<eval_0 )	{
				//	ASSERT(  Objectives( tag,MLGP::MEASURE,-1,0x0 )==eval );
					eval_0 = eval;		memcpy( tag_0,tag,sizeof(int)*nNode );	
					if( P_tags!=NULL )
						UpdateBestPass( eval,nP,P_tags,P_eval,flag );
					isBetter = 1;
				}
				for( j = ptr[cand]; j < ptr[cand+1]; j++ )	{
					no = adj[j];
					if( tag[no]==1 )
						continue;
					if( stmp[no]!=stp )	{
						stmp[no]=stp;
						update[nUpdate++] = no;
					}
				}
		//	}
			for( j = 0; j < nUpdate; j++ )	{
				no = update[j];
				if( tag[no]==1 )
					continue;
				UpdateGain_E( no,NULL,(float)base,NULL,0x0 );	
				list.UpdateWeight( no,gain[no],0x0 );
			}
		}

		if( T_1==0 || nNode-T_1==0 )
			continue;
		
		if( isBetter==1 )	{
			memcpy( tag,tag_0,sizeof(int)*nNode );	
			eval = Objectives( tag,MLGP::MEASURE,-1,0x0 );	
			ASSERT( eval==eval_0 );
			ret = GORD_OK;
		}
	}
	if( ret==GORD_OK )	{
		memcpy( tag,tag_0,sizeof(int)*nNode );	
		eval = Objectives( NULL,MLGP::MEASURE,-1,0x0 );		ASSERT( eval==eval_0 );
#ifdef _DEBUG
	//	s_0 = (float)(nEdge*100.0/(nNode*(nNode-1)));
	//	printf( "====Partition G=(%d,%d),(%d,%d,%g),r=%g%%\r\n",nNode,nEdge,T_1,nNode-T_1,eval,nEdge*100.0/(nNode*(nNode-1)) );
#endif		
	}
	delete[] tag_0;			delete[] update;		delete[] frontier;
	
	return ret;
}

/*
	Greedy Region Growing
	v0.1	cys
		6/12/2012
*/
int MLGP::P3_GreedyGrow( int *nP,int *P_tags,float *P_eval,int flag ){
	if( nNode<2 )
		return GORD_TOO_SMALL;
//	ASSERT( isSymmetry(0x0)==1 );
	if( nP!=NULL )		*nP = TABU_OUTER;
	int ret=GORD_VSEP_FAIL,i=0,n_1,n_0=0,nSeed,no,*tag_0,seed,nRand,l;
	float eval,eval_0=FLT_MAX;

	tag_0=new int[nNode];
	nRand = MIN( 20,nNode/2 );
	nSeed =	MAX( nNode/3,2 );
	//	srand( (unsigned)time( NULL ) );
	for( l = 0; l < nRand; l++ )	{
		for( i = 0; i < nNode; i++ )	{
			tag[i]=2;		gain[i] = 0;
		}
		n_0=0;		n_1 = 0;
		while( n_1<nSeed )	{		
			seed = G_RAND( )%nNode;			
			tag[seed]=1;		n_1++;
			for( i=ptr[seed]; i < ptr[seed+1];	i++ )	{
				no = adj[i];
				if( tag[no]==2 )
					tag[no]=0;
			}
		}
	//	if( n_1==0 || n_1+n_0==nNode )
	//		continue;
		eval = Eval( tag,0x0 );			
		VFM_Refine( &eval,VFM_NO_ABSORB );
		if( eval<eval_0 )	{
			eval_0 = eval;
			memcpy( tag_0,tag,sizeof(int)*nNode );	
			eval = Eval( tag,0x0 );	
			ASSERT( eval==eval_0 );
			ret = GORD_OK;
		}
	}
	if( ret==GORD_OK )	{
		memcpy( tag,tag_0,sizeof(int)*nNode );	
		eval = Eval( tag,0x0 );		ASSERT( eval==eval_0 );
#ifdef _DEBUG
	//	s_0 = (float)(nEdge*100.0/(nNode*(nNode-1)));
	//	printf( "====Partition G=(%d,%d),(%d,%d,%g),r=%g%%\r\n",nNode,nEdge,T_1,nNode-T_1,eval,nEdge*100.0/(nNode*(nNode-1)) );
#endif		
		for( i = 0; i < nMaxPart+1; i++ )
			wPart[i]=0.0;
		for( i = 0; i < nNode; i++ )	{
			no = tag[i];		ASSERT( no>=0 && no<=nMaxPart );
			wPart[no]+=wNode[i];
		}
		for( i = 1; i < nMaxPart+1; i++ )
			ASSERT( wPart[i]>0.0 );
	}
	delete[] tag_0;			

	return ret;
}

/*
	v0.1	cys
		5/9/2012
*/	
float MLGP::Objectives( int *part_0,int type,int cand, int flag ){
	ASSERT( sepa==SEPARATOR_E );
	float eval =FLT_MAX,cut=0.0,V_avg=0.0,w_1=0.0,cut_0,size_1,size_0;
	float w[3]={	wPart[0],wPart[1],wPart[2]	};
	int i,j,S_n=0,no;
	int *part=part_0==NULL ? tag : part_0;
	
	for( i = 0; i < nMaxPart+1; i++ )	{
		w[i] = wPart[i];				
	}	
	if( cand!=-1 )	{
		int from = tag[cand],to = from==1 ? 2 : 1;
		for( j = ptr[cand]; j < ptr[cand+1]; j++ )	{
			no = adj[j];
			if( tag[no]==from )
			{	w[0]+=2*wEdge[j];		w[from]-=2*wEdge[j];			}
			else
			{	w[0]-=2*wEdge[j];		w[to]+=2*wEdge[j];				}
		}
		cut = w[0];
		if( BIT_TEST( flag,UPDATE_WPART ) )	{
			for( i = 0; i < nMaxPart+1; i++ )	{
				wPart[i] = w[i];				
			}
		}
	}else	{
		wBorder = 0;
		for( i = 0; i < nMaxPart+1; i++ )	{
			wPart[i]=0;				
		}
		cut_0 = 0;				
		for( i = 0; i < nNode; i++ )	{
			no = part[i];			//ASSERT( no>0 && no<=nMaxPart );
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
	//			w_1 += wEdge[j];
				if( part[adj[j]]!=no )
				{	cut += wEdge[j];				}
				else
					wPart[no]+=wEdge[j];
			}
		/*	if( cut>cut_0 )	{		//is Border
				wBorder += w_1;
				cut_0 = cut;
			}*/
		}
		wPart[0] = cut;
		for( i = 0; i < nMaxPart+1; i++ )	{
			w[i] = wPart[i];				
		}		
	}
	size_1=MAX(w[1],w[2]),size_0=MIN(w[1],w[2]);
	switch( type )	{
	case NORMAL_CUT:
		if( w[1]>0 && w[2]>0 )	{
			eval = cut/w[1]+cut/w[2];
		}
		break;
	case MIN_CUT:
		if( size_0>0  /*&& size_1<size_0*4*/ )	{
			w_1 = size_1/size_0;
			eval = cut*(1.0+w_1*w_1 );
//			eval = size_1<size_0*2 ? cut : cut*(1.0+size_1/size_0);
		}
		break;
	default:
		eval =FLT_MAX;
	break;
	}
	ASSERT( eval>=0 );
	return eval;
}

/*
	v0.1	cys
		5/9/2012
*/
int MLGP::Kernel_Refine( int flag )	{
	int ret = GORD_OK,nMove=0;
	int nC=10,nP=0,*cluster,i,j,no,no_j,*P_no=new int[nNode],*map=new int[nNode];
	float *K,*W,*K_row;
	float e_0 = Objectives( tag,MLGP::MEASURE,-1,0x0 ),e_1=0;
	bool isBorder=false;
	for( i = 0; i < nNode; i++ )	{
		P_no[i] = -1;		map[i]=-1;
		no = tag[i];
		isBorder=false;
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			if( tag[adj[j]]!=no )
			{	isBorder=true;		break;	}			
		}
		if( isBorder )	{
			map[i] = nP;
			P_no[nP++] = i;
		}
	}
	cluster = new int[nP];
	K = new float[nP*nP];		W=new float[nP];
	for( i = 0; i < nP; i++ )	{
		K_row  = K + i*nP;
		for( j = 0; j < nP; j++ )	K_row[j]=0.0;
		W[i] = 0.0;
		no = P_no[i];			ASSERT( map[no]==i );
		for( j = ptr[no]; j < ptr[no+1]; j++ )	{
			no_j = map[adj[j]];
			if( no_j==-1 )
				continue;
			ASSERT( no_j>=0 && no_j<nP );
			K_row[no_j] = wEdge[j];
			W[i] += K_row[no_j];
		}
		cluster[i] = tag[no]-1;
	//	cluster[i] = G_RAND( )%nC;
		ASSERT( W[i]>0.0 );
	}
	for( i = 0; i < nP; i++ )	{
		K_row  = K + i*nP;
		for( j = 0; j < nP; j++ )	{
			K_row[j]/=(W[i]*W[j]);
		}
		K_row[i] += 50/W[i];
	}
	nMove = Kernel_Kmeans( nC,nP,cluster,K,W,0x0 );
	if( nMove>0 )	{
		for( i = 0; i < nP; i++ )	{
			no = P_no[i];
			tag[no] = cluster[i]+1;
		}
		e_1 = Objectives( tag,MLGP::MEASURE,-1,0x0 );
	//	ASSERT( e_0>e_1 );
	}
	delete[] cluster;
	delete[] K;					delete[] W;
	delete[] P_no;				delete[] map;

	return ret;
}

/*
	v0.1	cys
		5/9/2012
*/
int MLGP::Refine_2( int type,int flag )	{
	int ret = GORD_OK,loop=0,l_max=1;
//	clock_t t_0=clock( );
	ASSERT( nMaxPart==2 );
	switch( type )	{
	case TABU_REFINE:
//		Tabu_Refine_0( flag );
		break;
	case KMEANS_REFINE:
		Kernel_Refine( 0x0 );
		break;
	default:
		break;
	};	

//	if( !BIT_TEST( flag,REFINE_APART_LOOP ) )
//		Apart_Refine( 0x0 );

	return ret;
}

/*
	v0.1	cys
		5/9/2012
*/
int MLGP::Flow_Refine_2( int flag )	{
	if( wPart[0]<10 )
		return 0x0;
	MaxFlow_push flow;
	int ret = GORD_FLOW_REFINE_FAIL,i,j,no,cur,next,part=-1,nz,f_0=0,f_1=0,f_2=0,nS=0,nX=0;
	int *B_s = new int[nNode],*B_x=new int[nNode],*cut=new int[nNode+1],*map=new int[nNode];
	bool isCut,isRefine=false;

	float eval = Objectives( NULL,MLGP::MEASURE,-1,0x0 );

	for( i = 0; i < nNode; i++ )	{
		map[i]=-1;		
	}
//	part = wPart[1]>wPart[2] ? 1 : 2; 
	part = wPart[1]>wPart[2] ? 2 : 1;		//�ƺ�����
	nz = 0;
	for( i = 0; i < nNode; i++ )	{
		no = tag[i];
		if( no==part )
			continue;
		isCut = false;
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			if( tag[adj[j]]!=no )
			{	isCut=true;	break;	}
		}
		if( isCut==false )
			continue;
			
		map[i]=nS;		B_s[nS++]=i;
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			cur = adj[j];
			if( tag[cur]==part )	{
				if( map[cur]==-1 )	{
					B_x[nX++]=cur;
					map[cur]=1;
				}
				nz++;
			}
		/*	for( k = ptr[cur]; k < ptr[cur+1]; k++ )	{
				if( tag[adj[k]]==part && map[adj[k]]==-1)	{
					B_x[nX++]=adj[k];
					map[adj[k]]=1;
				}
			}
			nz += ptr[cur+1]-ptr[cur];*/
		}
	//	nz += ptr[i+1]-ptr[i];
	}
	for( i = 0; i < nX; i++ )	{
		map[B_x[i]]=i+nS;
		B_s[nS+i] = B_x[i];
	}
//	ASSERT( nS+nX==nV[0] );
//	goto END;
	if( nS+nX>500 )	{
//		goto END;
	}
	if( nz==nS && nz==nX )
		goto END;
	if( nS+nX==nNode )
		goto END;
	int F_nV=nS+nX+2;
	if( F_nV<5 )	
		goto END;
	nz = nz*2+2*(nS+nX);
	int *F_ptr=new int[F_nV+1],*F_adj=new int[nz],*F_cap=new int[nz],F_nz=0;
	F_ptr[0]=0;
	for( i = 0; i < nS; i++ )	{		//source node
		F_adj[F_nz]=i+1;		F_cap[F_nz]=(int)(wNode[B_s[i]]);
		F_nz++;					
		f_0 += (int)(wNode[B_s[i]]);
	}
	F_ptr[1]=F_nz;
	for( i = 0; i < nS+nX; i++ )	{
		cur =B_s[i];
		if( i<nS )	{			//separator
			F_adj[F_nz]=0;				F_cap[F_nz]=0;
		}else	{				//X part
			F_adj[F_nz]=nS+nX+1;		F_cap[F_nz]=(int)(wNode[cur]);
		}
		F_nz++;
		for( j = ptr[cur]; j < ptr[cur+1]; j++ )	{
			next = map[adj[j]];
			if( next==-1 )
				continue;
			if( (i<nS && next<nS) || (i>=nS && next>=nS) )
				continue;
			F_adj[F_nz]=next+1;		
			F_cap[F_nz]= i<nS ? CAPCITY_INFINITE : 0;
			F_nz++;
		}
		F_ptr[i+2]=F_nz;
	}
	for( i = nS; i < nS+nX; i++ )	{		//sink node
		F_adj[F_nz]=i+1;		F_cap[F_nz]=0;
		F_nz++;
	}
	F_ptr[nS+nX+2]=F_nz;
	ASSERT( F_nz <= nz );	
	flow.Init( F_nV,F_ptr,F_adj,F_cap,0x0 );
	f_1 = flow.Push_Relabel( 0x0,0x0 );
//	printf( "Flow_Refine_2: <%d,%d,%d>,%d=>%d\r\n",nS,nX,F_nz,f_0,f_1 );
	if( f_1<f_0 )	{
		isRefine=true;
		flow.GetCut( cut,NULL,0x0 );
		for( i = 0; i < nS+nX; i++ )	{
			cur = B_s[i];			
			if( cut[i+1]==1 )	{
				tag[cur]=0;
				f_2 += (int)(wNode[cur]);
			}else	{
/*				if( tag[cur]==0 )	{
					tag[cur]=part==2 ? 1 : 2;
				}else	{
				}*/
			}
		}
		ASSERT( f_2==f_1 );
	}
	CheckPartition( 0x0 );
	delete[] F_ptr,		delete[] F_adj,			delete[] F_cap;
	ret = GORD_OK;
END:
	if( !isRefine )	{
		for( i = 0; i < nS; i++ )	{
//		for( i = nS; i < nS+nX; i++ )	{
			cur =B_s[i];
			tag[cur] = 0;
		}

	}
	delete[] map;		
	delete[] B_s;		
	delete[] B_x;
	delete[] cut;
	return ret;
}

/*
	
	v0.1	cys
		5/14/2012
*/
int MLGP::UpdateGain_E( int no,int *x,float base,float *obj,int flag )	{
	int j,ret = GORD_OK,neibor,old=tag[no];
	ASSERT( x==NULL && nMaxPart==2 );
	float dCut=0.0;
//	float w[3]={wPart[0],wPart[1],wPart[2]};
	int from = tag[no],to = from==1 ? 2 : 1;
		
	for( j = ptr[no]; j < ptr[no+1]; j++ )	{
		neibor = adj[j];
		if( tag[neibor]==old )
			dCut+=wEdge[j];
		else
			dCut-=wEdge[j];
	/*	if( tag[neibor]==from )
		{	w[0]+=2*wEdge[j];		w[from]-=2*wEdge[j];			}
		else
		{	w[0]-=2*wEdge[j];		w[to]+=2*wEdge[j];				}*/
	}
/*	float size_1=MAX(w[1],w[2]),size_0=MIN(w[1],w[2]),w_1;
	if( size_0>0  && obj!=NULL  )	{
		w_1 = size_1/size_0;
		*obj = w[0]*(1.0+w_1*w_1 );
	}*/

	gain[no] = -dCut;
	gain[no] = gain[no]/gain_unit_edge+base;
	return ret;
}

/*
	
	v0.1	cys
		5/25/2012
*/
int MLGP::Tabu_SelectCand( SORT_BUCKET *moves,SORT_BUCKET **move_cand,int *modify,int *tenure,float *obj,float eval_0,int stp,int flag )	{
	int cand=-1,ring=1,part,head_bin,nM_0,nM,i,no,nAsp=0,f_M=0;
	SORT_BUCKET *hMove=NULL;
	float e_M,eval;

	*move_cand=NULL;
	if( stp%2==0 )	{
		part = wPart[1]>wPart[2] ? 2 : 1;
		ring = 1;
	}else	{
		int h_1=moves[0].Ceil_Bin( ),h_2=moves[1].Ceil_Bin( );
		part = h_1>h_2 ? 2 : 1;
		ring = 0;
	}
	hMove = moves+part-1;
	if( hMove->isEmpty( 0x0 ) )		{
		part = part==2 ? 1 : 2;
		hMove = moves+part-1;
		if( hMove->isEmpty( 0x0 ) )
			goto END;
	}

	if( CAND_TABU==2 )	{
		cand = hMove->First( );
		goto END;
	}
//selection
	head_bin = hMove->Ceil_Bin( );		ASSERT( head_bin>=0 );
	cand = -1;		e_M=FLT_MAX;
	do	{
		nM_0 = hMove->FirstRankList( &head_bin,modify,0x0 );		//ASSERT( nM_0>0 );
		nM = nM_0;
		nM = 0;			
		for( i = 0; i < nM_0; i++ )	{
			no = modify[i];
			if( tenure[no]>iter )	{	//check aspiration
				eval = Objectives( NULL,MLGP::MEASURE,no,0x0 );
				if( eval>=eval_0 )
				{	continue;		}	
				else
					nAsp++;
			}else	{
				eval = Objectives( NULL,MLGP::MEASURE,no,0x0 );
				if( eval==FLT_MAX )
					continue;
			}
		//	ASSERT( eval==obj[no] );
			if( eval<e_M )	{
				e_M =eval;		cand = no;
			}
			modify[nM++] = no;
		}
		head_bin--;
	}while( nM==0 && head_bin>=0 );
	if( nM==0 )
		goto END;
//	cand = nM>1 ? modify[G_RAND( )%nM] : modify[0];//�����������ƽ����ԣ�����֣�
//	ASSERT( gain[cand]>=G_thrsh );
	ASSERT( cand>=0 && cand<nNode );
END:
	*move_cand=hMove;
	return cand;
}
int MLGP::Tabu_Refine_0( float *tabu_obj,int flag )	{
	int ret = GORD_OK,i,j,no,cand=-1,T_move=0,from,to,nAsp=0,nCut=0;
	int *T_tag=new int[nNode],*modify=new int[nNode],stp=1,nM=0,nPertub=0,isMove;
	int *tenure=new int[nNode];//,*freq=new int[nNode];
	int gama = MAX( (int)(0.01*nNode),10 ),T_most=(int)MIN(1.1*nNode,1000),ring=0;
	int useGeneralSeed=!BIT_TEST( flag,TABU_ONLY_CUT );
	float g_1,g_0,eval_0=0,eval,G_thrsh=0.0,base=0;

	nTaRf++;
	nPertub = MAX((int)(0.02*nNode),5 );
	SORT_BUCKET *moves=new SORT_BUCKET[nMaxPart],*hMove=NULL;
	g_1=0,			g_0=FLT_MAX;
	for( i = 0; i < nNode; i++ )	{
		tenure[i]=-1;				
		gain[i] = -FLT_MAX;				//gainԽ��Խ��
		stmp[i] = stp;		
	}
	memcpy( T_tag,tag,sizeof(int)*nNode );
	base=float(gain_limit/gain_unit_edge);
	G_thrsh = base/2;
	for( i = 0; i < nMaxPart; i++ )
		moves[i].Init( (int)(base*2+1),nNode,0x0 );
	eval_0 = Objectives( NULL, MLGP::MEASURE,-1,0x0 );
	for( i = 0; i < nNode; i++ )	{
		from = tag[i];		to = from==1 ? 2 : 1;
		isMove = 0;
		if( useGeneralSeed )	{
			UpdateGain_E( i,NULL,base,NULL,0x0 );
			if( gain[i]>G_thrsh )
			{	isMove=1;		}	
		}	else	{
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
				if( tag[adj[j]]!=tag[i] )
				{	isMove=1;	break;}
			}
		}
		if( isMove==1 )	{
			UpdateGain_E( i,NULL,base,NULL,0x0 );
			moves[to-1].Insert( (int)(gain[i]),i,0x0 );		T_move++;
			g_1 = MAX( g_1,gain[i] );			g_0 = MIN( g_0,gain[i] );

		}else
			gain[i] = -100*gain_limit;
	}
	ASSERT( g_1<=gain_limit+base && g_0>=-gain_limit+base );

	T_most = nNode<100 ? nNode : (int)(TABU_rMove*nNode);
//	T_most = MIN( T_most,T_move*10 );
	iter = 0;				T_move=0;
	do	{
		iter++;
		stp++;				T_move++;
		cand = Tabu_SelectCand( moves,&hMove,modify,tenure,NULL,eval_0,stp,0x0 );
		if( cand==-1 )
			break;	
		from = tag[cand];				to = from==1 ? 2 : 1;
		eval = Objectives( NULL,MLGP::MEASURE,cand,UPDATE_WPART );	
		if( CAND_TABU==0 )
		{	ASSERT( tenure[cand]<=iter || eval<eval_0 );	}
//		tenure[cand]+=0.1*nNode+G_RAND()%2;		freq[cand]++;
		tenure[cand]+=T_most/2+G_RAND()%2;		//freq[cand]++;
		tag[cand] = to;
		hMove->Remove( (int)(gain[cand]),cand,0x0 );			gain[cand] = -100*gain_limit;
		if( eval<eval_0 )	{
			eval_0 = eval;			memcpy( T_tag,tag,sizeof(int)*nNode );
			gama = 0;
		}else	{
			gama++;
		}
		nM = 0;
		for( j = ptr[cand]; j < ptr[cand+1]; j++ )	{
			no = adj[j];
			if( stmp[no]==stp )
				continue;
			stmp[no]=stp;
			modify[nM++]=no;
		}
		if( gama>10 )	{	//Perturbation mechanism
		
		}		
		for( j=0; j < nM; j++ )	{
			no = modify[j];			//ASSERT( tag[no]==0 );
			if( no==151 )
				no = 151;
			from = tag[no];		to = from==1 ? 2 : 1;
			if( gain[no]>0 )
				moves[to-1].Remove( (int)(gain[no]),no,0x0 );
			UpdateGain_E( no,NULL,base,NULL,0x0 );			
			if( gain[no]>G_thrsh )
				moves[to-1].Insert( (int)(gain[no]),no,0x0 );
			else	{
				gain[no] = -100*gain_limit;
			}
		}
	}while( T_move<T_most );	
	MLGP::rMove += T_move*1.0/nNode;
	memcpy( tag,T_tag,sizeof(int)*nNode );	
#ifdef _DEBUG
	eval = Objectives( NULL,MLGP::MEASURE,-1,0x0 );	
	ASSERT( eval==eval_0 );
#endif
	*tabu_obj = eval_0;

	delete[] moves;
	delete[] T_tag;			delete[] modify;
	delete[] tenure;		

//	printf( "\tRefine:	time=%d,gain=%d,move=%d\r\n",clock( )-t_0,t_gain,nMove );

	return ret;
}


/*
	v0.1	cys
		10/5/2011
*/
int MLGP::NestSectionOrder( int *px,int x,int alg_,int level,int flag )	{
	int ret = GORD_OK,nCom=0,com=0,j,base=0,nSub=0;
	

	G_RANDm_init( );
	srand( 0x0 );		//ȷ��ÿ�εĽ��һ��  cys	6/6/2012
	MLGP **arrCom=NULL;
	int *map=new int[nNode],*sub=new int[nNode];	
	nLU = 0;
	double t_2 = G_TIC();
	nCom = SplitComponent( &arrCom,map,0x0 );	
	if( nCom>1 )
		nCom = nCom;
	for( com = 0; com < nCom; com++ )	{
		MLGP *child=NULL,*parent=NULL,*cur=NULL;	
		cur=arrCom[com];	
		bool isMMD = cur->nNode<MD_SWITCH;	
		if( !isMMD )	{
			MLGP *g_1=NULL,*g_2=NULL;
			int i,nS=0,nP_1=0,nP_2=0,code=0,no;
//			if( cur->BiSection( &nS,&nP_1,&nP_2,0x0,0x0 )==GORD_OK )	{
			if( cur->BiSection_2( &nS,&nP_1,&nP_2,0x0,level,0x0 )==GORD_OK )	{
//				int *s_0=new int[nNode],*s_1=new int[nNode],*s_2=new int[nNode];	
				MLGP::rEval += nS*1.0/nNode;
				ASSERT( nS+nP_1+nP_2==cur->nNode );
				ASSERT( /*nS>0 &&*/ nP_1>0 && nP_2>0 );				
				g_1 = cur->Sub( 1,cur->tag,&code,MLGP_CONNECTED );				ASSERT( code==GORD_OK );
				g_2 = cur->Sub( 2,cur->tag,&code,MLGP_CONNECTED );				ASSERT( code==GORD_OK );			
				if( g_1!=NULL && g_2!=NULL && MLGP::DUMP>1 )	{
					printf( "(%d-%d)=>[%d,%d,%d]\r\n",nNode,nEdge,nS,g_1->nNode,g_2->nNode );
				}
				g_1->NestSectionOrder( px,x,alg_,level+1,flag );
				g_2->NestSectionOrder( px,x,alg_,level+1,flag );
				int nzS=0,nzP_1=0,nzP_2=0;
				for( i = 0; i < cur->nNode; i++ )	{
					no = cur->tag[i];
					if( no==0 )	{
					//	cur->moveto[i]=nP_1+nP_2+g_0->moveto[nzS];
						cur->moveto[i]=(nP_1+nP_2+nzS);
						cur->kind[i] = 0;
					//	cur->moveto[i]=cur->nNode-1-(nzS);		//����
						nzS++;
					}else if( no==1 )	{
						cur->moveto[i] = g_1->moveto[nzP_1];
						cur->kind[i] = g_1->kind[nzP_1];
					//	cur->moveto[i] = nP_1-1-g_1->moveto[nzP_1];		//����
						nzP_1++;
					}else	{		//no=2
						cur->moveto[i] = g_2->moveto[nzP_2]+g_1->nNode;
						cur->kind[i] = g_2->kind[nzP_2];
					//	cur->moveto[i] = nP_1+nP_2-1-g_2->moveto[nzP_2];		//����
						nzP_2++;
					}
				}
				ASSERT( VERIFY_PERMUTE_VECTOR( cur->nNode,cur->moveto,sub ) == 0 );	
	//			delete[] s_0;		delete[] s_1;			delete[] s_2;
	//			nLU += g_1->nLU+g_2->nLU+nzS*(nzS+nzP_1+nzP_2);
				nLU += g_1->nLU+g_2->nLU+nzS*nzS;
			}else	
				isMMD = true;
			if( g_1!=NULL )		delete g_1;
			if( g_2!=NULL )		delete g_2;
		}
		if( isMMD  )	{
			nMDBlock++;
			for( int i = cur->nNode-1; i >=0; i-- )	{
				cur->moveto[i] = i;
				cur->kind[i] = nMDBlock;
			}
/*			cur->MdOrder( flag );
			if( nCom>1 )
				nLU += cur->nLU;
			else
				ASSERT( nLU>0 );*/
		}
		nSub=0;
		for( j = 0; j < nNode; j++ )	{
			if( map[j]==com+1 )
				sub[nSub++] = j;
		}
		ASSERT( nSub==cur->nNode );
		for( j = 0; j < nSub; j++ )	{
			moveto[sub[j]] = cur->moveto[j]+base;
			kind[sub[j]]=cur->kind[j];
		}
		base += cur->nNode;
		
		if( nCom>1 )
			delete cur;
	}

	if( arrCom!=NULL )	delete[] arrCom;
	delete[] map;			
	delete[] sub;		
	//printf( "\tNestSectionOrder:isMMD=%d,time=(%.3g,%.3g)\r\n",isMMD,G_TIC( )-t_0,t_r );
	
	return ret;
}

int MLGP::GetBlock( int B_dim,int **B_ptr,int **B_adj,int *map,int stp,int flag )	{
	ASSERT(map!=NULL);
	int ret=GORD_OK,nDad,i,j,no,no_j,nz,B_width,B_nz=0;
//	mask = temp;
	int *mask = tag;	//new int[nNode];	

	B_width = 0;		B_nz=0;
	for( i = 0; i < B_dim; i++ )	{
		no = map[i];
		B_nz +=ptr[no+1]-ptr[no];
		stmp[no]=stp;		mask[no] = B_width++;
	}
	int *C_ptr=new int[B_dim+1],*C_adj=new int[B_nz];
	B_nz = 0;
	C_ptr[0] = 0;
	for( i = 0; i < B_dim; i++ )	{
		no = map[i];
		for( j = ptr[no]; j < ptr[no+1]; j++ )	{
			no_j = adj[j];
			if( stmp[no_j]<stp )
			{	continue;	}
			C_adj[B_nz++] = mask[no_j];
			if(mask[no_j]<0 || mask[no_j]>=B_dim){
				ASSERT(0);
			}
		}
		C_ptr[i+1] = B_nz;
	}
	ASSERT( B_width==B_dim );
	*B_ptr = C_ptr;			*B_adj = C_adj;
// ASSERT(nz==B_ptr[nSet]);
//         for(i=0;i<nSet;i++){
//           for(j=B_ptr[i];j<B_ptr[i+1];j++){
//             if(B_adj[j]<0 || B_adj[j]>=nSet)
//               ASSERT(B_adj[j]>=0 && B_adj[j]<nSet);
//           }
//         }
	return B_nz;
}

/*
	v0.1	cys
		6/5/2012
	v0.2	cys
		6/20/2012
*/
int MLGP::GetXBlock_withcut( int B_dim,int *B_w,int **B_ptr,int **B_adj,int *map,int stp,int flag )	{
	int ret=GORD_OK,nDad,i,j,no,no_j,nz,B_width,B_nz=0;
	MLGP *hBasic=this,*hDad=NULL,*hCur=NULL;
	bool isMap=map==NULL;
	while( hBasic->parent!=NULL )	{
		hBasic = hBasic->parent;
	}
//	mask = temp;
	int *mask = tag;	//new int[nNode];
	if( isMap )	{
		map=new int[nNode];
		for( i = 0; i < B_dim; i++ )
			map[i]=i;

		hCur=this;
		while( hCur->parent!=NULL )	{
			hDad = hCur->parent;
			nDad = hDad->nNode;
			nz = 0;
			for( i = 0; i < nDad; i++ )		mask[i]=-1;
			for( i = 0; i < B_dim; i++ )	{
				if( i>0 )	{ASSERT( map[i]>map[i-1] );}
				mask[map[i]]=i;
			}
			for( i = 0; i < nDad; i++ )	{
				no = hCur->project[i];
				if( no==-1 )
					continue;
				if( mask[no]>=0 )	{
					map[nz++] = i;
				}
			}
			ASSERT( nz==B_dim );
			hCur = hDad;
		}
	}else	{
	}

//	for( i = 0; i < hBasic->nNode; i++ )	{
//		mask[i]=-1;
//	}
	B_width = 0;		B_nz=0;
	for( i = 0; i < B_dim; i++ )	{
		no = map[i];
		B_nz +=hBasic->ptr[no+1]-hBasic->ptr[no];
		stmp[no]=stp;		mask[no] = B_width++;
	}
	int *C_ptr=new int[B_dim+1],*C_adj=new int[B_nz];
	B_nz = 0;
	C_ptr[0] = 0;
	for( i = 0; i < B_dim; i++ )	{
		no = map[i];
		for( j = hBasic->ptr[no]; j < hBasic->ptr[no+1]; j++ )	{
			no_j = hBasic->adj[j];
			if( stmp[no_j]<stp )
			{	stmp[no_j] = stp;	mask[no_j] = B_width++;	}
			C_adj[B_nz++] = mask[no_j];
		}
		C_ptr[i+1] = B_nz;
	}
	ASSERT( B_width>=B_dim );
	*B_w = B_width;
	*B_ptr = C_ptr;			*B_adj = C_adj;

#ifdef _DEBUG
	/*	nz = 0;
		for( i = 0; i < B_dim; i++ )	{
			for( j = C_ptr[i]; j < C_ptr[i+1]; j++ )	{
				no = C_adj[j];
				if( no<B_dim )
					ASSERT( no==adj[nz++] );	
			}
		}*/
#endif
	if( isMap )
		delete[] map;
//	delete[] mask;

	return B_nz;
}

/*
	v0.1	cys
		6/15/2012
*/
int MLGP::Match_4( int *match,int *l_nIso,int loop_max,float z_thrsh,int flag )	{
	int nC=nNode,nRun,dim=nNode,i,j,cur,m_1,m_2,len,pos;
	int *c_m=tag;
	for( i = 0; i < nNode; i++ )	{
		match[i] = i;
	}
	do{
		nC=0;			nRun=dim;
		for( i = 0; i < dim; i++ )	{
			c_m[i]=-1;			
		}	
		while( nRun>0 )	{
			nRun--;
			cur = G_RAND( )%nNode;
			m_1 = match[cur];		
			if( c_m[m_1]!=-1 )
				continue;
			len = ptr[cur+1]-ptr[cur];
			for( j = ptr[cur]; j < ptr[cur+1]; j++ )	{
				pos = ptr[cur]+G_RAND( )%len;
				m_2 = match[adj[pos]];			
				if( c_m[m_2]==-1 )	{
					c_m[m_1] = c_m[m_2] = nC;
					nC++;
					break;
				}
			}
		}
		for( i = 0; i < dim; i++ )	{
			if( c_m[i]==-1 )
				c_m[i]=nC++;
		}
		for( i = 0; i < nNode; i++ )	{
			cur = match[i];
			ASSERT( cur>=0 && cur<dim );
			match[i] = c_m[cur];
			ASSERT( match[i]>=0 && match[i]<nC );
		}
		dim = nC;
	}while( nC>nNode*z_thrsh );
	return nC;
}

/*
	v0.1	cys
		6/15/2012
*/
bool MLGP::_HEM_match( int cur,int *m_1,int *m_2,int *match,int flag )	{
	*m_1=-1;		*m_2=-1;		
	float w_0=0.0,w_1=0.0,a,s;
	int cur_i = cur,j,next;
//	for( i = -3; i <=3; i++ )	{
//		cur_i = (cur+i+dim)%dim;
		if( match[cur_i]!=-1 )
			return false;		
		a = wNode[cur_i];
	/*	if( DOWN==3 )	{
			j = ptr[cur_i]+G_RAND( )%(ptr[cur_i+1]-ptr[cur_i]);
			next = adj[j];
			if( match[next]==-1 && a+wNode[next]<=wNode_max)	{
				*m_1=cur_i;				*m_2=next;	
				return true;
			}
		}*/
		for( j = ptr[cur_i]; j < ptr[cur_i+1]; j++ )	{
			next = adj[j];
			if( match[next]!=-1 )	{
				continue;
			}
			if( a+wNode[next]>wNode_max )
				continue;	
			if( DOWN==3 )	{		//RANDOM MATCH
				*m_1=cur_i;				*m_2=next;	
				break;
			}
			s = wEdge[j];	
			if( s>w_0 )	{
				*m_1=cur_i;				*m_2=next;	
				w_0 = s;				w_1=wNode[next];
			}else if( s==w_0 )	{
				if( wNode[next]<w_1 )	{
					*m_1=cur_i;				*m_2=next;	
					w_1 =  wNode[next];		
				}
			}
		}
//	}

	return *m_1!=-1;
}

int MLGP::PickConnectComponent(int seed,int C_no,int *C_id,int T_c,int flag){
	int nC = 0,i,j,cur,top=0,tail=0,no;
	int *stack=new int[nNode];
	stack[tail++] = seed;	
	while(top<tail){
		cur = stack[top++];
		C_id[cur] = C_no;
		if(top>T_c)
			break;
		for(j=ptr[cur];j<ptr[cur+1];j++){
			no = adj[j];
			if(C_id[no]!=-1)
				continue;

			stack[tail++] = no;	
		}
	}	
	delete[] stack;
	nC = top;
	return nC;
}

int MLGP::Match_N( int *l_match,int *l_nIso,int T_c,int flag )	{
	int dim=nNode,nRun=dim,nC=0,l_nC=0,cur,i,j,nIso=0,m_1,m_2,nz=0;
	int *left=new int[dim],*match=l_match,nLeft=0,loop,*l_project=NULL;
	double eval,eval_0=FLT_MAX;

	nLeft = 0;			nC=0;			nRun=dim;
	for( i = 0; i < dim; i++ )	{
		match[i]=-1;			
	}

	while( nRun>0 )	{//HEM matching	
		nRun--;
		cur = G_RAND( )%dim;
		if( match[cur]!=-1 )
			continue;		
	//	if( cur==1600 )
	//		cur = 1600;
		ASSERT( wNode[cur]<=wNode_max );	
		nz = PickConnectComponent(cur,nC,match,T_c,0x0);		
		nC++;		
	}
	for( i = 0; i < dim; i++ )	{
		if( match[i]!=-1 )
			continue;
		left[nLeft++]=i;
	}
	/*ASSERT( nLeft==dim-nC*2 );*/
	for( i = 0; i < nLeft; i++ )	{
		cur = G_RAND( )%nLeft;
		cur = left[cur];
		if( match[cur]!=-1 )
			continue;
		if( _HEM_match( cur,&m_1,&m_2,match,0x0 )  )	{
			ASSERT( match[m_1]==-1 && match[m_2]==-1 && cur==m_1 );
			match[m_1] = match[m_2] = nC;
			nC++;
		}else	{
			match[cur] = nC++;
		}
	}
	nIso = 0;
	for( i = 0; i < nLeft; i++ )	{
		cur = left[i];			ASSERT( cur>=0 && cur<dim );			
		if( match[cur]==-1 )
		{	match[cur] = nC++;		nIso++;	}
	}

	*l_nIso = nIso;

	
	l_nC = nC;

	delete[] left;			
	printf("\tMLGP::Match_N dim=%d nC=%d(%.5g) nIso=%d\n",dim,nC,dim*1.0/nC,nIso);
	return l_nC;
}

/*
	[REF]:	Comparison of coarsening schemes for multilevel	graph partitioning

	v0.1	cys
		6/13/2012
*/
MLGP* MLGP::ProjectDown_2( int *code,int loop_max,int flag )	{
	
	*code = GORD_OK;
	MLGP* child = NULL,*cur=NULL;
	int i,no,*match=NULL,*c_match=NULL,dim,nC=0,nIso=0;
	int nDense=0;
	float w_0=0,c_thrsh=0.95f,z_thrsh=0.2f;
//	ASSERT( isSymmetry( 0x0 )==1 );
	match = new int[nNode];		
	loop_max = MIN( loop_max,DWON_LOOP_MIN );
	for( i = 0; i < nNode; i++ )	{
		match[i] = i;
	}
	cur = this;
	do{
		
		dim = cur->nNode;
		c_match = new int[nNode];
	//	nC = project_core_( dim,cur->ptr,cur->adj,c_proj,cur->wNode,cur->wEdge,cur->wNode_max,0x0 );
		// if(nNode>300000){
		// 	//nC = cur->Match_N( c_match,&nIso,100,0x0 );
		// 	nC = cur->Match_( c_match,&nIso,loop_max,0x0 );
		// }else{
			if( DOWN==4 && nNode>THRSH_COARSEST*500 )	
				nC = cur->Match_4( c_match,&nIso,loop_max,z_thrsh,0x0);
			else	{
				nC = cur->Match_( c_match,&nIso,loop_max,0x0 );
			}
		//}
		
			
			
		if( nC==dim )
			break;
		if( STAR_0==-1.0 )	{
			STAR_0=nC*1.0/dim;
		}
		if( cur==this )
			STAR=MAX( STAR,nC*1.0/nNode );
		
		if( nC<=THRSH_COARSEST || nC<dim*c_thrsh )	{
			
			for( i = 0; i < nNode; i++ )	{
				no = match[i];
				ASSERT( no>=0 && no<dim );
				match[i] = c_match[no];
			}
			if( cur!=this )	{
				delete cur;
			}
			for( i = 0; i < nNode; i++ )	{
				no = match[i];
				ASSERT( no>=0 && no<nC );
				c_match[i] = no;
			}
			
			child = new MLGP( this,nC,c_match,type,0x0 );
			
			cur = child;
			
		}	else	{	
			break;
		}
	}while( child!=NULL && nC>nNode*z_thrsh && nC>THRSH_COARSEST );
	if( child==NULL )	{
		*code = GORD_DOWN_FAIL;		
	}
	
#ifdef _DEBUG
	//	printf( "====G=(%d,%d),ProjectDown cls_1=%g%%,cls_0=%g%%\r\n",nNode,nEdge,cls_1,cls_0 );
#endif//	if( nIso<nNode )	
	delete[] match;

	if(  child!=NULL )	{	
#ifdef _DEUB	//verify absorb down
			double w_sum = 0.0,w_0=0.0;
			for( i = 0; i < child->nNode; i++ )	w_sum += child->wNode[i];
			ASSERT( w_sum==wNode_p );
			w_sum = 0.0;		w_0=0.0;
			for( i = 0; i < nNode; i++ )		w_sum+=wAbsorb[i];
			for( i = 0; i < nEdge; i++ )		w_sum+=wEdge[i];
			for( i = 0; i < child->nNode; i++ )		w_0+=child->wAbsorb[i];
			for( i = 0; i < child->nEdge; i++ )		w_0+=child->wEdge[i];
			ASSERT( w_sum==w_0 /*&& w_sum==wEdge_p*/);
#endif
	}

	return child;
}

int MLGP::Match_( int *l_match,int *l_nIso,int loop_max,int flag )	{
	int dim=nNode,nRun=dim,nC=0,l_nC=0,cur,i,j,nIso=0,m_1,m_2;
	int *left=new int[dim],*match=NULL,nLeft=0,loop,*l_project=NULL;
	double eval,eval_0=FLT_MAX;

	if( loop_max>1 )	{
		match=new int[dim];
	}else	{
		match = l_match;
	}
	for( loop=0;loop <loop_max;loop++ )	{
		nLeft = 0;			nC=0;			nRun=dim;
		for( i = 0; i < dim; i++ )	{
			match[i]=-1;			
		}
	/*	for( i = dim-*l_nIso; i < dim; i++ )	{
			if( match[i]!=-1 )
				continue;		
			if( _HEM_match( i,&m_1,&m_2,match,0x0 )  )	{
				match[m_1] = match[m_2] = nC;
				nC++;
			}
		}*/
		while( nRun>0 )	{//HEM matching	
			nRun--;
			cur = G_RAND( )%dim;
			if( match[cur]!=-1 )
				continue;		
		//	if( cur==1600 )
		//		cur = 1600;
			ASSERT( wNode[cur]<=wNode_max );		
			if( _HEM_match( cur,&m_1,&m_2,match,0x0 )  )	{
				match[m_1] = match[m_2] = nC;
				nC++;
			}else	{
				m_1 = -1;
			}
		}
		for( i = 0; i < dim; i++ )	{
			if( match[i]!=-1 )
				continue;
			left[nLeft++]=i;
		}
		ASSERT( nLeft==dim-nC*2 );
		for( i = 0; i < nLeft; i++ )	{
			cur = G_RAND( )%nLeft;
			cur = left[cur];
			if( match[cur]!=-1 )
				continue;
			if( _HEM_match( cur,&m_1,&m_2,match,0x0 )  )	{
				ASSERT( match[m_1]==-1 && match[m_2]==-1 && cur==m_1 );
				match[m_1] = match[m_2] = nC;
				nC++;
			}else	{
				match[cur] = nC++;
			}
		}
		nIso = 0;
		for( i = 0; i < nLeft; i++ )	{
			cur = left[i];			ASSERT( cur>=0 && cur<dim );			
			if( match[cur]==-1 )
			{	match[cur] = nC++;		nIso++;	}
		}
		if( loop_max>1 )	{
			eval = 0.0;
			for( i = 0; i < nNode; i++ )		{
				ASSERT( match[i]>=0 &&  match[i]<nC  );
				for( j = ptr[i]; j < ptr[i+1]; j++ )	{
					if( match[i]!=match[adj[j]] )
						eval += wEdge[j];
				}
			}
			if( eval<eval_0 )	{
				eval_0 = eval;			
				l_nC = nC;				*l_nIso = nIso;
				memcpy( l_match,match,sizeof(int)*nNode );
			}
		}else	{
			*l_nIso = nIso;
		}
	}
	if( loop_max>1 )	{
		delete[] match;
	}else	{
		l_nC = nC;
	}
	delete[] left;			

	return l_nC;
}

/*
	v0.1	cys
		6/15/2012
*/
int MLGP::Pertub_2( int type,int *temp,int flag )	{
	int i,j,k,nP=0,nP_max=0,no;
	int *cut = temp,nCut=0;
	for( j = 0; j < nNode; j++ )	{
		if( tag[j]==0 )	{	
			//tag[j]=	G_RAND()%2+1;
			tag[j]=	wPart[2]>wPart[1] ? 1 : 2;
			nP_max++;
		}else	{
			for( k=ptr[j]; k < ptr[j+1]; k++ )	{
				if( tag[adj[k]]!=tag[j] )	{
					cut[nCut++] = j;
					nP_max++;		break;
				}
			}
		}
	}
	nP_max = MIN( nP_max,nNode/10 );
	nP_max = MAX( nP_max,50 );
	int part = wPart[1]>wPart[2] ? 1 : 2;		
	for( i = 0; i < nCut/10; i++ )	{
		no = G_RAND( )%nCut;
		no = cut[no];
		if( tag[no]==part )
			tag[no] = tag[no]==1 ? 2 : 1;
		nP++;
	}
	if( flag==0 )	{
		for( i = 0; i < nP_max; i++ )	{
			no = G_RAND( )%nNode;
		//	tag[no] = tag[no]==1 ? 2 : 1;
			nP++;
		}
	}

	return nP;
}

/*only valid for slim matrices
MLGP* MLGP::ProjectDown_slim( int *code,int slim_dim,int flag )	{	
	*code = GORD_OK;
	MLGP* child = NULL,*cur=NULL;
	int i,j,no,*match=NULL,*c_match=NULL,dim,nC=0,nIso=0;
	int nDense=0,nMatch=0;
	float w_0=0,c_thrsh=0.95f,z_thrsh=0.2f;
//	ASSERT( isSymmetry( 0x0 )==1 );
	match = new int[nNode];		
	for( i = 0; i < nNode; i++ )	{
		match[i] = i;
	}
	for( i = 0; i < nNode; i++ ) {
		for(j=ptr[i];j<ptr[i+1];j++){
			no = adj[j];
			if(match[no]==0){
				match[no]==1;		nMatch++;
			}
		}
		if( nMatch>nNode/2 )
			break;
	}
	for( i = 0; i < nNode; i++ )	{
		if(match[i] == 0 )
			match[i] = 2;
	}
	cur = this;
	do{
		
		dim = cur->nNode;
		c_match = new int[nNode];
	//	nC = project_core_( dim,cur->ptr,cur->adj,c_proj,cur->wNode,cur->wEdge,cur->wNode_max,0x0 );
		// if( DOWN==4 && nNode>THRSH_COARSEST*500 )	
		// 	nC = cur->Match_4( c_match,&nIso,loop_max,z_thrsh,0x0);
		// else
		// 	nC = cur->Match_( c_match,&nIso,loop_max,0x0 );
		if( nC==dim )
			break;
		if( STAR_0==-1.0 )	{
			STAR_0=nC*1.0/dim;
		}
		if( cur==this )
			STAR=MAX( STAR,nC*1.0/nNode );
		
		if( nC<=THRSH_COARSEST || nC<dim*c_thrsh )	{
			
			for( i = 0; i < nNode; i++ )	{
				no = match[i];
				ASSERT( no>=0 && no<dim );
				match[i] = c_match[no];
			}
			if( cur!=this )	{
				delete cur;
			}
			for( i = 0; i < nNode; i++ )	{
				no = match[i];
				ASSERT( no>=0 && no<nC );
				c_match[i] = no;
			}
			
			child = new MLGP( this,nC,c_match,type,0x0 );
			
			cur = child;
			
		}	else	{	
			break;
		}
	}while( child!=NULL && nC>nNode*z_thrsh && nC>THRSH_COARSEST );
	if( child==NULL )	{
		*code = GORD_DOWN_FAIL;		
	}
	
#ifdef _DEBUG
	//	printf( "====G=(%d,%d),ProjectDown cls_1=%g%%,cls_0=%g%%\r\n",nNode,nEdge,cls_1,cls_0 );
#endif//	if( nIso<nNode )	
	delete[] match;

	if(  child!=NULL )	{	
#ifdef _DEUB	//verify absorb down
			double w_sum = 0.0,w_0=0.0;
			for( i = 0; i < child->nNode; i++ )	w_sum += child->wNode[i];
			ASSERT( w_sum==wNode_p );
			w_sum = 0.0;		w_0=0.0;
			for( i = 0; i < nNode; i++ )		w_sum+=wAbsorb[i];
			for( i = 0; i < nEdge; i++ )		w_sum+=wEdge[i];
			for( i = 0; i < child->nNode; i++ )		w_0+=child->wAbsorb[i];
			for( i = 0; i < child->nEdge; i++ )		w_0+=child->wEdge[i];
#endif
	}

	return child;
}*/

/*
	v0.1	cys
		10/4/2011
	v0.2	cys
		6/15/2012
*/
int	MLGP::BiSection_2( int *nS_,int *nP_1_,int *nP_2_,int alg_,int level,int flag )	{
	double t_0=G_TIC( );
	int ret=GORD_OK,d_ret,i,j,loop_proj=MAX(level,1),o_best=0,i_best;
	int *tag_0 = new int[nNode],*tag_i=new int[nNode],T_coarsest=0;
	MLGP *child=NULL,*parent=NULL,*cur=this,*g_base;
	float s_0,eval,eval_0,obj_0,obj;
	bool isPartition=false;
	int inner_refine[10]={
		TABU_ONLY_CUT,
		TABU_ONLY_CUT,
		TABU_ONLY_CUT,
		TABU_ONLY_CUT,
		0x0,
		TABU_ONLY_CUT,0x0,TABU_ONLY_CUT,0x0,TABU_ONLY_CUT
	};
	ASSERT( TABU_INNER<=10 );

//	if( THRSH_COARSEST<0 )
//		T_coarsest=MD_SWITCH*2;
//	else
		T_coarsest=THRSH_COARSEST;
	nBisection++;
	s_0 = nEdge*1.0f/nNode;
	//printf("\tMLGP::BiSection_2:ProjectDown...\n");
	while( cur->nNode>T_coarsest )	{//down
		if( DOWN==0 )	{
			child = cur->ProjectDown( &d_ret,0x0 );		
		}	else
		{	child = cur->ProjectDown_2( &d_ret,loop_proj,0x0 );		loop_proj+=1;			}
		
		if( d_ret!=GORD_OK )
		{	ASSERT( child==NULL );		break;	}
		nDown++;
	//	s = child->nEdge*1.0/child->nNode;
	//	if( s > s_0*1.5 )
	//		break;
		cur = child;
	}
	//printf("\tMLGP::BiSection_2:Init Partition...\n");
//Init Partition
	g_base = cur;
	int P_dim=cur->nNode,nPass=TABU_OUTER,*P_tag=NULL;// new int[P_dim*TABU_OUTER];
	float *P_eval = NULL;	//new float[TABU_OUTER];
//	for( i = 0; i < TABU_OUTER; i++ )	P_eval[i]=FLT_MAX;
	if( MLGP::SEPRATOR==MLGP::SEPARATOR_E )	{
		isPartition=cur->Partition_2( &nPass,P_tag,P_eval,0x0 )==GORD_OK;
		ASSERT( isPartition );
	}else
		isPartition=cur->Partition_3( 0x0 )==GORD_OK;

	if( !isPartition && P_tag!=NULL )	{
		delete[] P_tag;		P_tag=NULL;
		delete[] P_eval;	P_eval=NULL;
		nPass=1;
	}
	eval_0 = FLT_MAX;	//cur->Objectives( NULL, MLGP::MEASURE,-1,0x0 );
	//printf("\tMLGP::BiSection_2:Partition_2...\n");
	for( i = 0; i < nPass; i++ )	{//up
	//	cur = g_base;
		if( i>=1 )	{
			cur->Pertub_2( 0x0,stmp,0x0 );
		}else 
			obj = FLT_MAX;
		if( 0 && P_eval!=NULL )	{
			if( P_eval[i]==FLT_MAX )
				continue;
			memcpy( cur->tag,P_tag+i*P_dim,sizeof(int)*P_dim );
			eval = cur->Objectives( NULL, MLGP::MEASURE,-1,0x0 );
			ASSERT( eval==P_eval[i] );
		}
		while( true )	{
			if( !isPartition )	{
				if( MLGP::SEPRATOR==MLGP::SEPARATOR_E )	{
					isPartition = cur->Partition_2( NULL,NULL,NULL,0x0 )==GORD_OK;
				}else	{
					isPartition = cur->Partition_3( 0x0 )==GORD_OK;
				}
			}else	{
				obj_0=FLT_MAX;		i_best=0;
				for( j = 0; j < TABU_INNER; j++ )	{
					if( MLGP::SEPRATOR==MLGP::SEPARATOR_E )	{
		
						cur->Tabu_Refine_0( &obj,inner_refine[j] );
		
				//		obj = cur->Objectives( NULL, MLGP::MEASURE,-1,0x0 );
						if( obj<obj_0 )	{
							obj_0 = obj;		i_best=j;
							memcpy( tag_i,tag,sizeof(int)*cur->nNode );							
						}			
	//					cur->Pertub_2( 0x0,stmp,1 );
						if( j>i_best )	{
							break;
						}
					}else	{
						cur->Refine_3( MLGP::VFM_REFINE,NULL,inner_refine[j] );		
					}	
				}
			}
			memcpy( tag,tag_i,sizeof(int)*cur->nNode );
		
			if( cur!=this )		{
				if( isPartition )
					cur->ProjectUp( 0x0 );	
				parent = cur->parent;		
				parent->CheckPartition( 0x0 );
				cur = parent;	
			}	else	{	
				eval = obj_0;
				if( eval<eval_0 )	{
					eval_0=eval;	
					memcpy( tag_0,tag,sizeof(int)*cur->nNode );
					o_best = i;
				}
//					Flow_Refine_2( 0x0 );
				if( i>o_best+1 )	{
					goto OUTER_NEXT;
				}
				break;
			}				
		}
	}	
OUTER_NEXT:
	
	memcpy( tag,tag_0,sizeof(int)*nNode );
	eval = Objectives( NULL, MLGP::MEASURE,-1,0x0 );//Eval( NULL, 0x0 );
	ASSERT( eval==eval_0 );
	double	t_1=G_TIC();	
	Flow_Refine_2( 0x0 );
	MLGP::tFlowRefine += G_TIC( )-t_1;
	
	Refine_3( MLGP::VFM_REFINE,stmp,0x0 );	
	
//	VFM_Refine( &eval,0x0 );
	cur = g_base;
	while( cur!=NULL && cur!=this )	{
		parent = cur->parent;	
		delete cur;					cur = parent;	
	}
	if( P_eval!=NULL )	delete[] P_eval;
	if( P_tag!=NULL )	delete[] P_tag;
	if( !isPartition )	{
		ret = GORD_BISEC_FAIL;
	}else	{
	}
	
/*	if( wPart[1]<wPart[2] )	{
		for( i = 0; i < nNode; i++ )	{
			if( tag[i]!=0 )
				tag[i] = tag[i]==2 ? 1 : 2;
		}
	}*/
	int no,nS=0,nP_1=0,nP_2=0;
	for( i = 0; i < nNode; i++ )	{
		no = cur->tag[i];
		if( no==0 )	{
			nS++;		//s_0[nS++] = i;
		}else if( no==1 )	{
			nP_1++;		//s_1[nP_1++] = i;
		}else	{
			ASSERT( no==2 );
			nP_2++;		//s_2[nP_2++] = i;
		}
	}	
	if( nS_!=NULL )			*nS_=nS;
	if( nP_1_!=NULL )		*nP_1_=nP_1;
	if( nP_2_!=NULL )		*nP_2_=nP_2;

	if( nP_1==0 || nP_2==0 )	{
		ret = GORD_BISEC_FAIL;
	}
	delete[] tag_0;		delete[] tag_i;
	MLGP::tBiSection += G_TIC()-t_0;
	if(level<=2){
		printf("\tMLGP::BiSection_2@L%d	T=%.3g\r\n",level, G_TIC()-t_0);
	}
	return ret;
}

