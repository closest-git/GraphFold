#include <memory.h> 
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "gord.h"
#include "Flow.h"

/*
	[1] Algorithms in C++,Part 5
	[2] Combinatorial Optimization Goldberg
*/
/*
	http://support.microsoft.com/?scid=kb;en-us;157623&amp;x=16&amp;y=12
	http://blog.csdn.net/sun_top/article/details/4213413
	typedef priority_queue<FLOW_NODE*,vector<FLOW_NODE*>,node_w_cmp> PQ_node;

class node_w_cmp{
public :
//strict weak ordering		http://support.microsoft.com/kb/949171
	bool operator()(FLOW_NODE* a,FLOW_NODE* b)
	{	return a->w>b->w;        }
};
*/
/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/24/2012	
*/
int FLOW_EDGE::R_Cap( FLOW_NODE * to )	{	
	return isBack(to) ? flow : cap-flow;	
}
/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/24/2012	
*/
void FLOW_EDGE::R_AddFlow( FLOW_NODE * to,int d )	{
	flow += isBack( to ) ? -d : d;
	ASSERT( flow<=cap );
};	
/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/24/2012	
*/
int FLOW_EDGE::R_Push( FLOW_NODE * from,int push_0,int flag )	{
	int push = 0;
	FLOW_NODE * to = from==pv ? pw : pv;
	bool isB = isBack(to);
	bool isDownHill = BIT_TEST(flag,MaxFlow_push::NO_CHECK_HEIGHT ) || from->height==to->height+1;
	int free = isB ? flow : cap-flow;		
	if( free>0 )	{
		ASSERT( from->height<=to->height+1 );		//valid distance labeling
		if( isDownHill )	{
			push = MIN( free,push_0 );
			flow += isB ? -push : push;
			from->w -= push;
			to->w += push;
		}
	}

	return push;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/24/2012	
*/
MaxFlow_push::MaxFlow_push( )	{
	FLOW_MAX=-1;

	nV=0,nE=0;
//	NETWORK_SOURCE=-1,		NETWORK_SINK=-1,		
//	ptr=NULL,				adj=NULL;	
//	height=NULL,			w=NULL;
	nodes=NULL;
	edges=NULL;
	actives.clear( );

	if( 0 )	{	//test the sample in 417
		int adj_[16]={1,2,0,3,4,0,3,4,1,2,5,1,2,5,3,4},
			cap_[]	={2,3,0,3,1,0,1,1,0,0,2,0,0,3,0,0};
		int ptr_[7]={0,2,5,8,11,14,16},flow;
		Init( 6,ptr_,adj_,cap_,0x0 );
		flow = Preflow_Push( 0x0,0x0 );
		flow = Push_Relabel( 0x0,0x0 );
	}
}
MaxFlow_push::~MaxFlow_push( )	{
//	if( ptr!=NULL )			delete[] ptr;
//	if( adj!=NULL )			delete[] adj;
//	if( height!=NULL )		delete[] height;
//	if( w!=NULL )			delete[] w;
	if( nodes!=NULL )		delete[] nodes;
	if( edges!=NULL )		delete[] edges;
}

/*
	v0.1	cys
		4/24/2012
*/
int MaxFlow_push::InitHeights( int flag )	{
	int ret = GORD_OK,base,i;
	queue<FLOW_NODE*> qn;
	FLOW_EDGE *hEdge=NULL;
	FLOW_NODE* hSink = SINK( ),*cur,*next;
	vector<FLOW_EDGE*>::iterator it;

	for( i = 0 ;i < nV; i++ )		{
		cur = nodes+i;
		cur->height=-1;
	}
	hSink->height=0;
	qn.push( hSink );
	while( !qn.empty( ) )	{
		cur = qn.front( );	
		qn.pop( );
		base = cur->height;
		for( it = cur->adj.begin( ); it != cur->adj.end( ); it++ )	{
			hEdge = *it;
			next = hEdge->Other( cur );
			if( next->height!=-1 )
				continue;
			next->height = base+1;
			qn.push( next );
		}
	}
	return ret;
}

/*
	[1]	P.418
	v0.1	cys
		4/24/2012
*/
int MaxFlow_push::Preflow_Push( int alg,int flag )	{
	int i,ret = GORD_OK,cap,push,flow=0,cut=0;
	FLOW_EDGE *hEdge=NULL;
	FLOW_NODE *hNode=NULL,*hSource=SOURCE(),*hSink=SINK( ),*cur=NULL,*next=NULL;
	vector<FLOW_EDGE*>::iterator it;
	InitHeights( 0x0 );
	for( i = 0; i < nV; i++ )	{
		hNode = nodes+i;
		hNode->w=0;
		hNode->flag=0x0;
		BIT_RESET( hNode->flag,FLOW_NODE::IN_QUEUE );
	}

	hSource->w = FLOW_MAX;	hSink->w = -hSource->w;
	actives.clear( );
	actives.push_back( hSource );		BIT_SET( hSource->flag,FLOW_NODE::IN_QUEUE );
	while( !actives.empty( ) )	{
		cur = actives.front( );	
		actives.pop_front( );			BIT_RESET( cur->flag,FLOW_NODE::IN_QUEUE );
		for( it = cur->adj.begin( ); it != cur->adj.end( ); it++ )	{
			hEdge = *it;
			next = hEdge->Other(cur);
			cap = hEdge->R_Cap( next );
			push = MIN( cap,cur->w );
			if( push>0 && ( hSource==cur || cur->height==next->height+1 ) )	{
				hEdge->R_AddFlow( next,push );
				cur->w-=push;		next->w+=push;
				if( hSource!=next && hSink!=next )	{
					if( !BIT_TEST( next->flag,FLOW_NODE::IN_QUEUE ) )
					{	actives.push_back( next );	BIT_SET( next->flag,FLOW_NODE::IN_QUEUE );	}
				}
			}
			if( cur->w==0 )
				break;
		}
		if( hSource!=cur && hSink!=cur && cur->w>0 )		{	
			cur->height++;			ASSERT( cur->height<2*nV );
			ASSERT( !BIT_TEST( cur->flag,FLOW_NODE::IN_QUEUE ) );
			actives.push_back( cur );	BIT_SET( cur->flag,FLOW_NODE::IN_QUEUE );
		}

	}
	ASSERT( hSink->w = -hSource->w );
	flow = FLOW_MAX-hSource->w;
	if( 1 )	{
		vector<FLOW_NODE*> path;
		GetAugmentPath( path,0x0,0x0 );
		cut = GetCut( NULL,NULL,0x0 );
		ASSERT( flow==cut );
	}
	return flow;
}

/*
	[2]	P.39
	v0.1	cys
		4/24/2012
*/
bool node_h_cmp(FLOW_NODE* a,FLOW_NODE* b){	return a->height<b->height;        }
int MaxFlow_push::Push_Relabel( int alg,int flag )	{
	int i,ret = GORD_OK,push,flow=0,cut=0,h_0;
	bool isRelabel=false;
	FLOW_EDGE *hEdge=NULL;
	FLOW_NODE *hSource=SOURCE(),*hSink=SINK( ),*cur=NULL,*next=NULL;
	deque<FLOW_NODE *> heap;
	vector<FLOW_EDGE*>::iterator it,itCur;

	for( i = 0; i < nV; i++ )	{
		nodes[i].Reset( );		
	}
	for( i = 0; i < nE; i++ )	{
		edges[i].Reset( );
	}
	flow = 0;			//Initialize		[2] P.29
	for( it = hSource->adj.begin( ); it != hSource->adj.end( ); it++ )	{
		hEdge = *it;
		flow += hEdge->R_Push( hSource,FLOW_MAX,NO_CHECK_HEIGHT );
		next = hEdge->Other(hSource);
		heap.push_back( next );		BIT_SET( next->flag,FLOW_NODE::P_ACTIVE );
	}
	hSource->height=nV;
	hSource->w = -flow;	hSink->w = 0;

	while( !heap.empty( ))	{
		make_heap ( heap.begin( ), heap.end( ),node_h_cmp );
		cur = heap[0];				BIT_RESET( cur->flag,FLOW_NODE::P_ACTIVE );		
		heap.pop_front( );
		isRelabel=false;
		h_0 = INT_MAX;
		itCur = cur->adj.begin( );
		for( it = itCur; it != cur->adj.end( ); it++ )	{
			hEdge = *it;
			next = hEdge->Other(cur);
			push = hEdge->R_Push( cur,cur->w,0x0 );
			if( push>0 && hSink!=next && next!=hSource )	{
				if( hSink==next )	{
				}else	if(next==hSource )	{
					ASSERT( hSource->w<0 );
				}else	if( !BIT_TEST( next->flag,FLOW_NODE::P_ACTIVE ) )
				{	heap.push_back( next );	BIT_SET( next->flag,FLOW_NODE::P_ACTIVE );	}
			}
			if( hEdge->R_Cap( next )>0 )				
				h_0 = MIN( h_0,next->height );
	/*		push = MIN( cap,cur->w );		ASSERT( push>0 );
			if( cur->height==next->height+1 )	{
				hEdge->R_AddFlow( next,push );
				cur->w-=push;		next->w+=push;
				if( hSink!=next && next!=hSource )	{
					if( !BIT_TEST( next->flag,FLOW_NODE::P_ACTIVE ) )
					{	heap.push_back( next );	BIT_SET( next->flag,FLOW_NODE::P_ACTIVE );	}
				}
			}*/
			if( cur->w==0 )
				break;
		}
		ASSERT( cur!=hSink );
		isRelabel = cur->w>0;
		if( isRelabel )	{
			ASSERT( h_0+1>cur->height && h_0+1<2*nV );
			cur->height = h_0+1;		
			heap.push_back( cur );		BIT_SET( cur->flag,FLOW_NODE::P_ACTIVE );
		}
	}
	flow = -hSource->w;			ASSERT( hSink->w==flow );
	if( 0 )	{
		vector<FLOW_NODE*> path;
		bool bAug = GetAugmentPath( path,0x0,0x0 );
		cut = GetCut( NULL,NULL,0x0 );
		ASSERT( flow==cut && bAug==false );
	}
	return flow;
}
/*
	http://msdn.microsoft.com/en-us/library/80tkwx9w.aspx
	defines sense in which one element is less than another
*/
bool node_w_cmp(FLOW_NODE* a,FLOW_NODE* b)
{	return a->w<b->w;        }
/*
	P.392
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/25/2012	
*/
bool MaxFlow_push::GetAugmentPath( vector<FLOW_NODE*>&path,int alg,int flag )	{
	deque<FLOW_NODE *> heap;
	vector<FLOW_EDGE *> st;
	vector<FLOW_EDGE*>::iterator it;
	FLOW_NODE *hSource=SOURCE(),*hSink=SINK(),*cur,*next;
	FLOW_EDGE *hEdge=NULL;
	int i,cap,P;
	for( i = 0; i < nV; i++ )	{
		cur = nodes+i;		cur->w=0;
		heap.push_back( cur );
	}
	hSource->w = FLOW_MAX;

	st.resize( nV );
	while( !heap.empty( ) )	{
		make_heap ( heap.begin( ), heap.end( ),node_w_cmp );
		cur = heap[0];			heap.pop_front( );
		if( cur==hSink  || (cur!=hSource && st[cur->no]==NULL ) )
			break;
		for( it=cur->adj.begin( );it!=cur->adj.end( );it++ )	{
			hEdge = (*it);
			next = hEdge->Other( cur );
			cap = hEdge->R_Cap( next );
			P = MIN( cap,cur->w );
			if( cap>0 && P>next->w )	{
				next->w=P;		st[next->no]=hEdge;
			}
		}
	}
	hEdge = st[hSink->no];
	bool bRet = hEdge!=NULL;
	cur = hSink;		path.push_back( cur );
	while( hEdge!=NULL )	{
		next = hEdge->Other( cur );		cur = next;
		hEdge = st[cur->no];		ASSERT( cur==hSource || hEdge!=NULL );
		path.push_back( cur );		
	}
	return bRet;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/25/2012	
*/
int MaxFlow_push::GetCut( int *c_nodes,int *c_edges,int flag )	{
	int i,n_1=0,n_c=0,cap=0;
	FLOW_NODE *hNode=NULL,*hSource=SOURCE(),*hSink=SINK( ),*cur=NULL,*next=NULL;
	vector<FLOW_EDGE*>::iterator it;
//	int *length=new int[nV],*part=new int[nV];
	bool isPush;
	for( i = 0; i < nV; i++ )	{	
		hNode = nodes+i;
		hNode->level=-1;		hNode->part=2;
//		part[i]=2;	length[i]=-1;		
	}

	stack<FLOW_NODE *> s_nodes;
	s_nodes.push( hSource );				hSource->level=0;
	FLOW_EDGE *hEdge=NULL;
	while( !s_nodes.empty( ) )	{
		cur = s_nodes.top( );				s_nodes.pop( );	
		cur->part=1;					n_1++;
		isPush = false;
		for( it = cur->adj.begin( ); it != cur->adj.end( ); it++ )	{
			hEdge = *it;
			next = hEdge->Other( cur );
			if( next->level==-1  && hEdge->R_Cap( next )>0 )				{	
				next->level=cur->level+1;
				s_nodes.push( next );		isPush=true;	
			}
		}
		if( isPush == false )	{
		}
	}
	ASSERT( hSink->part==2 );
	if( c_nodes!=NULL )	{
		for( i = 0; i < nV; i++ )
		{		c_nodes[i]=-1;	}
	}
	if( c_edges!=NULL )	{
		for( i = 0; i < nE; i++ )
		{		c_edges[i]=-1;	}
	}
	for( i = 0; i < nV; i++ )	{
		cur = nodes+i;
		for( it = cur->adj.begin( ); it != cur->adj.end( ); it++ )	{
			hEdge = *it;
			next = hEdge->Other( cur );
			if( cur->part!=next->part )			{	
				n_c++;		
				if( cur->part==1 && cur==hEdge->Source( ) )		{		//st edges
					ASSERT( hEdge->R_Cap( next )==0 );
					cap+=hEdge->Capcity( );		
				//	cap+=hEdge->Flow( );			
					if( c_edges!=NULL )
						c_edges[(int)(hEdge-edges)]=1;
				}else	{						//ts edges
				}
				if( c_nodes!=NULL )	
					c_nodes[cur->no]=1;						
			}
		}
	}

	return cap;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/25/2012	
*/
FLOW_EDGE *MaxFlow_push::FindEdge( FLOW_NODE * pv_,FLOW_NODE * pw_,int flag )	{
	FLOW_EDGE *hEdge=NULL;
	vector<FLOW_EDGE*>::iterator it;
	for( it = pv_->adj.begin( ); it != pv_->adj.end( ); it++ )	{
		if( (*it)->Other(pv_)==pw_ )	{
			hEdge = (*it);	break;
		}
	}
	ASSERT( hEdge!=NULL );
	return hEdge;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/25/2012	
*/
int MaxFlow_push::Check( int flag )	{
	int i,ret=GORD_OK;
	FLOW_NODE *v,*w;
	FLOW_EDGE *hEdge=NULL;
	for( i = 0; i < nE; i++ )	{
		hEdge = edges+i;
		v = hEdge->Source( );		w = hEdge->Destination( );
	/*	for( j = ptr[v]; j < ptr[v+1]; j++ )	{
			if( adj[j]==hEdge )	{
				v=-1;		break;
			}
		}
		for( j = ptr[w]; j < ptr[w+1]; j++ )	{
			if( adj[j]==hEdge )	{
				w=-1;		break;
			}
		}
		ASSERT( v==-1 && w==-1 );*/
	}

	return ret;
}


/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/25/2012	
*/
int MaxFlow_push::Init( int nG,int *G_p,int *G_adj,int *G_cap,int flag )	{
	int i,j,nzE=0,ret = GORD_OK;
	FLOW_NODE *hNode=NULL,*cur=NULL,*next=NULL;
	nV = nG;
	nodes = new FLOW_NODE[nV];
	for( i = 0; i < nV; i++ )	{
		hNode = nodes+i;
		hNode->no = i;
	}
/*	
	NETWORK_SOURCE=0;			NETWORK_SINK=nV-1;
	height = new int[nV];		w = new int[nV];
	for( i = 0; i < nV; i++ )	{
		height[i] = -1;
		w[i]=-1;
	}
	ptr = new int[nV+1];
	adj = new FLOW_EDGE*[G_p[nG]];
*/
	nE = G_p[nG]/2;
	edges = new FLOW_EDGE[nE];

	FLOW_MAX=0x0;
//	ptr[0] = 0;
	FLOW_EDGE *hEdge=NULL;
	for( i = 0; i < nG; i++ )	{	
		cur = nodes+i;
		for( j = G_p[i]; j < G_p[i+1]; j++ )	{
			ASSERT( G_adj[j]>=0 && G_adj[j]<nG );
			next = nodes+G_adj[j];
			if( G_cap[j]==0 )	{	//�ñ��ѳ���
				ASSERT( G_adj[j]<i );
				hEdge=FindEdge( next,cur,0x0 );						
			}else	{
				ASSERT( G_adj[j]>i);
				hEdge = edges+nzE;
				edges[nzE++].Init( G_cap[j],cur,next,0x0 );	
			}
			cur->adj.push_back(hEdge);				ASSERT( hEdge!=NULL );	
			if( G_cap[j]!=CAPCITY_INFINITE )
				FLOW_MAX+=G_cap[j];
			ASSERT( FLOW_MAX>0 );
		}
//		ptr[i+1]=G_p[i+1];
	}
	ASSERT( nzE==nE );
	FLOW_MAX*=100;

	Check( 0x0 );

	return ret;
}


