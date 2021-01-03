#pragma once

typedef long long G_INT_64;

#define BIT_SET( val,flag ) ((val) |= (flag))	
#define BIT_RESET( val,flag ) ((val) &= (~(flag)) ) 
#define BIT_TEST( val,flag ) (((val)&(flag))==(flag))
#define BIT_IS( val,flag ) (((val)&(flag))!=0)

#define MIN(a,b) ( (a)<(b)?(a):(b) )
#define MAX(a,b) ( (a)>(b)?(a):(b) )

#define SWAP(a,b,t) { (t)=(a);	(a)=(b);	(b)=(t); }

#define G_DOUBLE2INT(a)			(G_INT_64)( (a)>=0.0 ? (a)+0.5 : (a)-0.5 )

#define GORD_OK 0x0
#define	GORD_TOO_SMALL	-8
#define GORD_VSEP_FAIL	-10
#define GORD_BISEC_FAIL	-11
#define	GORD_DOWN_FAIL	-12
#define	GORD_PERTUB_FAIL	-13
#define	GORD_FLOW_REFINE_FAIL	-15
#define	GORD_UNKNOWN_ERROR	-100

#define GORD_MAT_A 0x100
#define GORD_MAT_L 0x200
#define GORD_MAT_SYMMETRY 0x400

#define _DEBUG
#ifdef _DEBUG
	#include <assert.h>
	#define ASSERT( a ) assert( (a) )
#else
	#define ASSERT( a )
#endif
/*
	//--------------------------------------------------------------------------------------
	// Thread safety 
	//--------------------------------------------------------------------------------------
	#include <Windows.h>
	CRITICAL_SECTION g_cs;  
	bool g_bThreadSafe = true; 
	class DXUTLock
	{
	public:
		inline DXUTLock()  { if( g_bThreadSafe ) EnterCriticalSection( &g_cs ); }
		inline ~DXUTLock() { if( g_bThreadSafe ) LeaveCriticalSection( &g_cs ); }
	};
*/
#ifdef __cplusplus
extern "C" {
#endif
#ifdef GORD_DLL_EXPORTS
	__declspec(dllimport) G_INT_64 ccs_MLGP_( G_INT_64 nCol,G_INT_64* ptr,G_INT_64* ind,G_INT_64 *weight,G_INT_64* pCol,G_INT_64 flag );
#else
	G_INT_64 ccs_MLGP_( G_INT_64 nCol,G_INT_64* ptr,G_INT_64* ind,G_INT_64 *weight,G_INT_64* pCol,G_INT_64 algorithm,G_INT_64 flag );
#endif


#ifdef __cplusplus
}
#endif