#pragma once
#include <float.h>
#include "gord.h"
//#include "../util/grus_util.h"

unsigned long genrand_int32(void);
#define G_RAND( )	(	(G_INT_64)(genrand_int32()>>1)	)
//#define G_RAND( )	(	(G_INT_64)((rand()<<3)|(rand() ) ) )
/*
	�����б�����Bucket SortΪ��׼
	���ۺϿ��� list�̶ܶ�ȡֵ�ܴ�����
*/
class SORT_BUCKET	{		//����node���� bucket sort
	G_INT_64 rank,root,size,type;
	G_INT_64 *next,*seq;	
//	union {
		G_INT_64 *head;
		G_INT_64 *value;
//	};
	enum{
		VALUE_VECTOR,HEAD_BIN
	};
public:
	SORT_BUCKET( )	{
		memset( this,0x0,sizeof(SORT_BUCKET) );		root=-1;
	}
	~SORT_BUCKET( );
	void Clear( )	{
		G_INT_64 i;
		for( i = rank-1; i>=0; i-- )		{	head[i]=-1;					}
		for( i = size-1; i>=0; i-- )		{	next[i]=-1;			}
	}

	G_INT_64 Ceil_Bin( )					{	return root;		}
	G_INT_64 Floor_Bin( )					{	
		G_INT_64 i,bin=-1;
		for( i = 0; i <= root; i++ )	{
			if( head[i]==-1 )
				continue;
			bin = i;		break;
		}
		return bin;
	}
	bool isEmpty( G_INT_64 flag )	{	return root==-1;	}
	G_INT_64 First( )	{	return root==-1? -1:head[root];	}
	G_INT_64 FirstRankList( G_INT_64 *list,G_INT_64 flag )	{
		G_INT_64 nz = 0;
		if( root!=-1 )	{
			G_INT_64 cur = head[root];
			while( cur!=-1 )	{
				list[nz++] = cur;
				cur = next[cur];
			}
		}
		return nz;	
	}
	G_INT_64 FirstRankList( G_INT_64 *head_bin,G_INT_64 *list,G_INT_64 flag )	{
		G_INT_64 nz = 0,i,cur;
		ASSERT( *head_bin<=root );
		for( i = *head_bin; i >=0; i-- )	{
			cur =head[i];
			if( cur==-1 )
				continue;
			else
				break;
		}
		*head_bin = i;
		while( cur!=-1 )	{
			list[nz++] = cur;
			cur = next[cur];
		}
		
		return nz;	
	}
	G_INT_64 Bin_List( G_INT_64 bin,G_INT_64 *list,G_INT_64 flag )	{
		ASSERT( bin<=root );
		G_INT_64 nz = 0,cur = head[bin];
		while( cur!=-1 )	{
			list[nz++] = cur;
			cur = next[cur];
		}		
		return nz;	
	}

	void Init( G_INT_64 R,G_INT_64 len,G_INT_64 flag );
	G_INT_64 operator[](G_INT_64 pos)	{	
		ASSERT( pos>=0 && pos<size );
		ASSERT( seq!=NULL );
		return seq[pos];			
	}

	G_INT_64 Sequence( G_INT_64 flag )	{
		G_INT_64 last=-1,i,cur,nnz=0,ret=0x0;
		seq = new G_INT_64[size];
		root = -1;
		for( i = rank-1; i >=0; i-- )	{
			cur =head[i];
			if( cur==-1 )
				continue;
			if( root==-1 )	root=cur;
			while( cur!= -1 )	
			{	seq[nnz++]=cur;	cur = next[cur];	}
		}
		ASSERT( nnz==size );
		return ret;
	}

	void Insert( G_INT_64 h_pos,G_INT_64 i,G_INT_64 flag );
	void Remove( G_INT_64 h_pos,G_INT_64 no,G_INT_64 flag );	
};


/*
*/
class SORT_LIST	{		//����node���� bucket sort
	G_INT_64 size;
	G_INT_64 *next,*prev,head;		
public:
	float *weight;

	SORT_LIST( )	{
		memset( this,0x0,sizeof(SORT_LIST) );		
	}
	~SORT_LIST( )	{
		if( next!=NULL )		delete[] next;
		if( weight!=NULL )		delete[] weight;
		if( prev!=NULL )		delete[] prev;
	}

	void Init( G_INT_64 len,G_INT_64 flag )	{
		size=len;
		next = new G_INT_64[len];			prev = new G_INT_64[len];
		weight = new float[len];
		Clear( );
	}
	void Clear( )	{
		G_INT_64 i;
		for( i = 0; i<size; i++ )	{
			next[i]=-1;			prev[i]=-1;
			weight[i] = -FLT_MAX;
		}
		head = -1;
		ASSERT( IsValid(0x0) );
	}

	void Sort( G_INT_64 flag )	{
	}

	void Remove( G_INT_64 pos )	{
		G_INT_64 p,n;
		ASSERT( pos>=0 && pos<size );
		if( weight[pos] == -FLT_MAX )	{
			ASSERT( prev[pos]==-1 && next[pos]==-1 );
			return;
		}
		if( prev[pos]!=-1 )	{		//isolate
			p = prev[pos];		n = next[pos];
			next[p] = n;
			if( n!=-1 )
				prev[n] = p;
			prev[pos]=-1;	
		}
		if( pos==head )	{
			n = next[head];
			head = n;	
			if( n!=-1 )
				prev[n]=-1;
		}
		next[pos]=-1;	
		weight[pos] = -FLT_MAX;		
		ASSERT( IsValid( 0x0 ) );
	}

	void UpdateWeight( G_INT_64 pos,float w,G_INT_64 flag )	{
		ASSERT( pos>=0 && pos<size );
		if( weight[pos]==w )
			return;
		G_INT_64 cur = head,p;
//		if( head==pos )	{}
		Remove( pos );	
		weight[pos] = w;
		cur = head;				p=-1;
		while( cur!=-1 )	{
			ASSERT( prev[cur]==p );
			if( w>=weight[cur] )	{				
				break;
			}
			p = cur;		
			cur = next[cur];		
		}
		if( p==-1 )	{		//head
			head = pos;			
		}else	{
			next[p] = pos;
		}
		prev[pos] = p;
		next[pos] = cur;	
		if( cur!=-1 )
			prev[cur] = pos;
		ASSERT( IsValid( 0x0 ) );
	}

	G_INT_64 First( G_INT_64 *cands,G_INT_64 flag )	{
		G_INT_64 cur = head,nM=0,cand=-1;
		if( head==-1 )
			return nM;

		float w_0=weight[head];
		while( cur!=-1 )	{
			if( weight[cur]==w_0 )
				cands[nM++]=cur;
			else
				break;
			cur = next[cur];		
		}
		return nM;		
	}

	bool IsValid( G_INT_64 flag )	{
		G_INT_64 cur = head,p=-1;
		while( cur!=-1 )	{
			ASSERT( cur>=0 && cur<size );
			ASSERT( prev[cur]==p );	
			ASSERT( cur!=p );
			p = cur;
			cur = next[cur];		
		}
		return true;
	}
};

bool G_RANDs_init( G_INT_64 size,G_INT_64 flag );

bool G_RANDm_init( );

G_INT_64 CCS_Compress_( G_INT_64 nCol,G_INT_64* ptr,G_INT_64* ind,G_INT_64 *weight,G_INT_64 *comp,G_INT_64 alg,G_INT_64 flag );
G_INT_64 CCS_Compress_2( G_INT_64 dim,G_INT_64* ptr,G_INT_64* ind,G_INT_64 *weight,G_INT_64 *map,G_INT_64 alg,G_INT_64 flag );