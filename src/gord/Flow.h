#pragma once
#include <queue>
#include <deque>
#include <stack>
#include <algorithm>
#include <functional>
//#include<bits/stdc++.h> 
using namespace std;

#define CAPCITY_INFINITE INT_MAX

/*
	�ο� 
		1 http://en.wikipedia.org/wiki/Flow_network
	ע��
		1 cap,flow������
*/
class FLOW_NODE;
class FLOW_EDGE{
	G_INT_64 flow;
//v==>w is Forward
	bool isBack( FLOW_NODE * to )	{	return to==pv;	}
public:
	enum{
		CAP_INFINITE=0x100,
	};
	FLOW_NODE *pv,*pw;
	G_INT_64 cap,flag;		//v->w

	FLOW_EDGE( )	{
		memset( this,0x0,sizeof(FLOW_EDGE) );
		cap=-1;	
	}

	void Init( G_INT_64 c,FLOW_NODE *v,FLOW_NODE *w,G_INT_64 f )	{
		cap=c,		pv=v,		pw=w,		flag=f;
	}

	FLOW_NODE * Source( )			{	return pv;	}
	FLOW_NODE * Destination( )		{	return pw;	}
	G_INT_64 Capcity( )			{	return cap;	}
	G_INT_64 Flow( )				{	return flow;	}
	void Reset( )			{	flow=0;			}
	G_INT_64 R_Push( FLOW_NODE * from,G_INT_64 push_0,G_INT_64 flag );

	FLOW_NODE * Other( FLOW_NODE * cur )	{
		ASSERT( cur==pv || cur==pw );
		return cur==pv ? pw : pv;
	}

	G_INT_64 R_Cap( FLOW_NODE * to );
	void R_AddFlow( FLOW_NODE * to,G_INT_64 d );	
};

/*
	ע��	w������ʹ��
	4/27/2012
*/
class FLOW_NODE{
public:
	enum{
		IN_QUEUE=0x10000,P_ACTIVE=0x2000,
	};
	G_INT_64 no;	
//union{
	G_INT_64 height;
	G_INT_64 level;
//};
//union{
	G_INT_64 w;
	G_INT_64 part;
//};
	G_INT_64 flag;
	vector<FLOW_EDGE*> adj;

	FLOW_NODE( )	{
		no=0,	flag=0;		height=0;		level=0;
		w=0;	part=0;
	}
	~FLOW_NODE( )	{
		adj.clear( );
	}

	void Reset( )	{	height=0;		w=0;		flag=0x0;	}
};

/*
	Preflow-Push Maxflow algorithms
*/
class MaxFlow_push
{
	G_INT_64 FLOW_MAX;

	G_INT_64 nV,nE;
//	G_INT_64 *ptr,*height,*w,NETWORK_SOURCE,NETWORK_SINK,;
	FLOW_NODE *nodes;
	FLOW_EDGE *edges;
	deque<FLOW_NODE *> actives;
protected:
//	bool isSource(G_INT_64 v)	{	return v==NETWORK_SOURCE;	}
//	bool isSink(G_INT_64 v)	{	return v==NETWORK_SINK;	}
	G_INT_64 InitHeights( G_INT_64 flag );
	FLOW_EDGE *FindEdge( FLOW_NODE * pv_,FLOW_NODE * pw_,G_INT_64 flag );

public:
	enum{
		NO_CHECK_HEIGHT=0x100
	};
	MaxFlow_push();
	~MaxFlow_push();
	
	FLOW_NODE *SOURCE( )	{	ASSERT( nV>=2 );	return nodes+0;	}
	FLOW_NODE *SINK( )	{	ASSERT( nV>=2 );	return nodes+nV-1;	}

	G_INT_64 Init( G_INT_64 nG,G_INT_64 *G_v,G_INT_64 *G_adj,G_INT_64 *G_cap,G_INT_64 flag );
	G_INT_64 Preflow_Push( G_INT_64 alg,G_INT_64 flag );
	G_INT_64 Push_Relabel( G_INT_64 alg,G_INT_64 flag );
	bool GetAugmentPath( vector<FLOW_NODE*>&path,G_INT_64 alg,G_INT_64 flag );
	G_INT_64 GetCut( G_INT_64 *c_nodes,G_INT_64 *c_edges,G_INT_64 flag );

	G_INT_64 Check( G_INT_64 flag );
//	G_INT_64 FindPath( G_INT_64 start,G_INT_64 end,G_INT_64 alg,G_INT_64 flag );
};

