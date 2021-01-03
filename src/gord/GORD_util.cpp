#include <memory.h> 
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "gord.h"
#include "GORD_util.h"
#include "MLGP.h"
extern clock_t g_tX;


/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/*
*/
bool G_RANDm_init( )	{
	unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
    init_by_array(init, length);

	return true;
}

void SORT_BUCKET::Remove( int h_pos,int no,int flag )	{
	if( head[h_pos]==no )			{	
		head[h_pos] = next[no];
		if( head[h_pos]==-1 )	{	
			if( h_pos==root )	{
				root = -1;
				for( int i = rank-1; i >= 0; i-- )	{
					if( head[i]!=-1 )	
					{	root = i;		break;	}
				}
			}
		}
	}	else	{
		int cur = head[h_pos];
	//	if( cur==-1 )	{	//5/15/2012		cys
	//	}else	{
			ASSERT( cur!=-1 );
			while( next[cur]!=-1 && next[cur]!=no )	{
				cur = next[cur];
			}
			ASSERT( next[cur]==no );
			next[cur] = next[no];
	//	}
	}
	next[no]=-1;
}

SORT_BUCKET::~SORT_BUCKET( )	{
	if( head!=NULL )	delete[] head;
	if( value!=NULL )	delete[] value;
	if( next!=NULL )	delete[] next;
	if( seq!=NULL )		delete[] seq;
}

/*
	v0.2	cys
		6/21/2012
*/
void SORT_BUCKET::Init( int R,int len,int flag )	{
	int i;
	rank = R;			root=-1;			size=len;
	next = new int[len];
	head = NULL;		value=NULL;
/*	if( 0 && R>len*20 )	{
		type = VALUE_VECTOR;
		value = new int[len];
		for( i = size-1; i>=0; i-- )	value[i]=-1;	
	}else*/	{
		type = HEAD_BIN;
		head = new int[rank];
		for( i = rank-1; i>=0; i-- )	head[i]=-1;
	}
	for( i = size-1; i>=0; i-- )	next[i]=-1;	
}

/*
	v0.2	cys
		6/21/2012
*/
void SORT_BUCKET::Insert( int h_pos,int i,int flag )	{
	ASSERT( i>=0 && i < size );
	ASSERT( h_pos>=0 && h_pos < rank );
	ASSERT( next[i]==-1 );
	if( head[h_pos]==-1 )	{	
		head[h_pos]=i;
		if( h_pos>root )
			root = h_pos;
	}else	{
		int cur = head[h_pos];
		next[i] = cur;
		head[h_pos] = i;
	}
}

void CCS2CRS_struc( int nRow,int nCol,int nnz,int* ccs_ptr, int *ccs_ind,int* rowptr,int *colind,int* temp	 )	{
	int i,col,row,isAlloc = temp==NULL;
	if( isAlloc )	{
		temp=new int[nRow];
	}
	memset( temp,0,sizeof(int)*nRow );
	for( i = 0; i < nnz; i++ )	
		temp[ccs_ind[i]]++;
	rowptr[0] = 0;
	for( i = 1; i < nRow+1; i++ )	{	
		rowptr[i] = rowptr[i-1]+temp[i-1];
		temp[i-1] = rowptr[i-1];
	}
	for( col = 0; col < nCol; col++ )	{
		for( i = ccs_ptr[col]; i < ccs_ptr[col+1]; i++ )	{
			row = ccs_ind[i];
			colind[temp[row]] = col;
			temp[row]++;
		}
	}
	if( isAlloc )	{
		delete[] temp;
	}
}

/* 
	get A+A' from A. CCS format
	ptr��֪!!!
	temp[nCol]

	v0.1	cys
		8/11/2004 CCS2AAT_struc( ) 
	v0.2	cys
		12/20/2004 can deal with no-diagonal
*/
int CCS2AAT_struc( int dim,int nnz,int* a_ptr,int* a_ind,int* ptr,int*ind,int*temp ) {
	int *cur,i,pL,pU,indL,indU,*colind,*rowptr;
	colind=(int*)malloc(nnz*sizeof(int)); 
	rowptr=(int*)malloc((dim+1)*sizeof(int));
	CCS2CRS_struc( dim,dim,nnz,a_ptr,a_ind,rowptr,colind,temp );

	cur=temp;
	for( i = 0;  i < dim ; i++)    cur[i] = ptr[i];

	for( i = 0; i < dim; i++ )	{
		pL = a_ptr[i]; 
		while( a_ind[pL] < i && pL < a_ptr[i+1] )	{ pL++;}
		pU=rowptr[i];
		while( colind[pU] < i && pU < rowptr[i+1] )	{ pU++;}
		if( colind[pU] > i || pU == rowptr[i+1])	{	//add diagonal
			ind[cur[i]]=i;	
			cur[i]++ ;	
		}

		while( pL < a_ptr[i+1] && pU < rowptr[i+1] )	{
			indL = a_ind[pL];	indU = colind[pU];	
			if( indL == indU )	{
				ind[cur[i]]=indL;	
				cur[i]++ ;	
				if( indL != i )		{	//non diagonal
					ind[cur[indU]]=i;					
					cur[indU]++ ;										
				}
				pL++;	pU++;
			}else if( indL < indU )	{	//Only in L
				ind[cur[indL]]=i;				ind[cur[i]]=indL;	
				cur[indL]++ ;					cur[i]++ ;						
				pL++;	
			}else	{					//Only in U
				ind[cur[indU]]=i;				ind[cur[i]]=indU;	
				cur[indU]++ ;					cur[i]++ ;						
				pU++;				
			}
		}
		while(  pL < a_ptr[i+1] )	{	//Only in L
				indL = a_ind[pL];	
				ind[cur[indL]]=i;				ind[cur[i]]=indL;	
				cur[indL]++ ;					cur[i]++ ;						
				pL++;	
		}
		while( pU < rowptr[i+1] )	{	//Only in U
				indU = colind[pU];	
				ind[cur[indU]]=i;				ind[cur[i]]=indU;	
				cur[indU]++ ;					cur[i]++ ;						
				pU++;				
		}

	}
	free(rowptr);	free(colind);
	
	for( i = 0; i < dim; i++ )
		ASSERT( cur[i] == ptr[i+1] );
    return 0;
}

/*
	v0.1	cys
		9/14/2005
	v0.2	cys
		6/26/2012
*/
int CCS_Compress_2( int dim,int* ptr,int* ind,int *weight,int *map,int alg,int flag )	{
	int i,j,k,key,cur,deg,nOff,stp=2,C_dim=0,thresh,top,nExtra=0;
	int *hash,*parent,*next,*mark,*stack,*C_ptr=NULL,*C_ind=NULL,C_nz=0,C_nz_max=0;
	double ksi = 0.05;

	hash=new int[dim*5];		
	parent=hash+dim;	next=parent+dim;		mark=next+dim;		stack=mark+dim;
	stp = 2;
	for( i = 0; i < dim; i++ )	{
		key=i;					//�Խ�Ԫ���Ǵ���
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			key += ind[j];
		}
		map[i]=-1;
		hash[i] = key;			mark[i] = stp;			
		parent[i] = -1;			next[i] = -1;
	}
	for(  i = 0; i < dim; i++ )	{
		if( parent[i]!=-1 )		
			continue;
		stp++;
		key = hash[i];			deg=ptr[i+1]-ptr[i];
		thresh = (int)(ksi*deg);
		mark[i] = stp;
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{ 
			cur = ind[j];		ASSERT( cur != i );
			mark[cur] = stp;
		}
		top=0;
		for( j = ptr[i]; j < ptr[i+1]; j++ )			{	
			if( parent[ind[j]]!=-1 )		continue;
			stack[top++] = ind[j];		parent[ind[j]]=i;	
		}
		while( top > 0 )	{
			cur = stack[--top];
			ASSERT( parent[cur]==i );
//			if( cur<=i || key!=hash[cur] || deg!=(ptr[cur+1]-ptr[cur]) )
			if( cur<=i || abs(deg-(ptr[cur+1]-ptr[cur]))>thresh )
			{		parent[cur]=-1;				continue;			}				
			nOff=0;
			for( k = ptr[cur]; k < ptr[cur+1]; k++ )	{
				if( mark[ind[k]] != stp )	{
					nOff++;
				}
			}
			if( nOff<=thresh )		{		//����CLUSTER NODE
				if( nOff > 0 )
					nExtra++;
				next[cur]=next[i];			next[i]=cur;			
				for( j = ptr[cur]; j < ptr[cur+1]; j++ )	{
					if( parent[ind[j]]!=-1 )		continue;
					stack[top++] = ind[j];			parent[ind[j]]=i;	
				}
				ASSERT( top <= dim );
			}else
				parent[cur]=-1;
		}
		if( parent[i] == -1 )	{
			C_dim++;
			C_nz_max += deg;
		}
	}
	if( C_dim==dim )	
	{	goto END;	}

	j = 0;
	for(  i = 0; i < dim; i++ )		{
		if( parent[i]==-1 )		
			map[i]=j++;
		else	{
			ASSERT( map[parent[i]] != -1 );
			map[i] = map[parent[i]];
		}
	}
	ASSERT( j==C_dim );
	
END:
	delete[] hash;	

	return 0;
}

/*
	v0.1	cys
		5/28/2012
*/
int CCS_Compress_( int nCol,int* ptr,int* ind,int *weight,int *comp,int alg,int flag )	{
	int i,j,k,no,same=0,nComp=0,nCls=0,nzMax=0,nz=0,nPass=0,slim[10]={0};
	int *stmp=new int[nCol],stp=0,*sum = new int[nCol];
	bool isSame;

	for( i = 0; i < nCol; i++ )	{	
		stmp[i]=stp;		comp[i]=-1;
		sum[i]=0;
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			no = ind[j];			
			sum[i]+=no;
		}
		if(ptr[i+1]-ptr[i]<10){		//some matrices is very slice!!!
			slim[ptr[i+1]-ptr[i]]++;
		}
	}
	float similar = 0.0,r_1,r_2;

	for( i = 0; i < nCol; i++ )	{
		if( comp[i]!=-1 )
			continue;
		nzMax += ptr[i+1]-ptr[i];
		comp[i]=nCls++;
		stp++;
		r_1=ptr[i+1]-ptr[i];
		if( r_1==0.0 )
			continue;	
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			no = ind[j];			
			stmp[no]=stp;
		}

		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			no = ind[j];
			if( no==i || comp[no]!=-1 )
				continue;
			r_2=ptr[no+1]-ptr[no];
			if( r_2==0.0 )
				continue;
	/*		if( sum[i]!=sum[no] )
				continue;
			if( r_1!=r_2 )
				continue;
			isSame=true;
			for( k = ptr[no]; k < ptr[no+1]; k++ )	{
				if( stmp[ind[k]]!=stp )
				{	isSame=false;		break;			}
			}*/
			same = 0;
			for( k = ptr[no]; k < ptr[no+1]; k++ )	{
				if( stmp[ind[k]]==stp )
				{	same++;			}
			}
			similar = ((same)/r_1+(same)/r_2)/2.0;		ASSERT( similar<=1.0);	
			isSame = similar > 0.95;
		//	isSame = similar ==1.0;
			if( isSame )	{
				comp[no]=comp[i];		nComp++;
			}else	{
				nPass++;
			}
			
		}
	};
	ASSERT( nComp+nCls==nCol );
	delete[] stmp;		delete[] sum;
	printf("CCS_Compress_(%d-%d) nSingle=(%d,%d，%d，%d）\r\n", nCol, ptr[nCol], slim[2],slim[3],slim[4],slim[5]);
	return nCls;
}