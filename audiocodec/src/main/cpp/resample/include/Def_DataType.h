//====================================================================
// File Name : Def_DataType.h
// Function  : Define the system Data type 
// Program   : Li, Shaoen (LSN)
// Date      : Sep 21, 2005
// Version   : 0.1
// History
//		Version		Author		Memo	
//		0.1			Lishaoen	Programming start
//====================================================================

#ifndef __DEF_DATATYPE_H_2005_09_21__
#define __DEF_DATATYPE_H_2005_09_21__


//====================================================================
//	Version definition
//====================================================================

#define		VER_MINUS	1

//====================================================================
//	Version definition end
//====================================================================

//====================================================================
//	fundamental datatype
//====================================================================
#define		CHAR		char
#define		INT8		char
#define		INT8S		char
#define		__int8		char
#define		UCHAR		unsigned char
#define		BYTE		unsigned char
#define		INT8U		unsigned char
#define		UINT8		unsigned char

#define		SHORT		short
#define		INT16		short
#define		INT16S		short
#define		__int16		short
#define		USHORT		unsigned short
#define		WORD		unsigned short
#define		INT16U		unsigned short
#define		UINT16		unsigned short

#define		DWORD		unsigned int
#define		INT			int
#define		INT32S		int
#define		__int32		int
#define		INT32		int
#define		UINT		unsigned int
#define		INT32U		unsigned int
#define		UINT32		unsigned int
#define		WPARAM		unsigned int
#define		LPARAM		unsigned int

#define		LONG		long
#define		ULONG		unsigned long

#ifdef	__MIPS__
#define		INT64		long long
#define		INT64S		long long
#define		INT64U		unsigned long long
#define		UINT64		unsigned long long
#define		LONGLONG	long long
#else
#define		INT64		__int64
#define		INT64S		__int64
#define		INT64U		unsigned __int64
#define		UINT64		unsigned __int64
#define		LONGLONG	__int64
#endif

#define		COLORREF	UINT

#define		FLOAT		float
#define		FP32		float
#define		DOUBLE		double
#define		FP64		double

#define		VOID		void
#define		VOIDPARA	void

//====================================================================
//	fundamental datatype end
//====================================================================


//====================================================================
//	pointer data type
//====================================================================

#define		LPADDR		char*
#define		LPINT		INT*
#define		LPUINT		UINT*
#define		PVOID		void*
#define		LPVOID		void*
#define		LPSTR		char*
#define		LPCSTR		const char*
#define		LPCTSTR		const char*
#ifndef		NULL
#define		NULL		0
#endif

//====================================================================
//	pointer data type end
//====================================================================

//====================================================================
//	functions result data type
//====================================================================

#define		LRESULT		int
#define		LR_OK 		0
#define		LR_ERR 		-1

//====================================================================
//	functions result data type end
//====================================================================

//====================================================================
//	logical data type
//====================================================================

#ifndef		BOOL
#define		BOOL		int
#endif
#define		TRUE 		1
#define		FALSE 		0

//====================================================================
//	logical data type end
//====================================================================

//====================================================================
//	macro of combine data
//====================================================================
#define		MAKEWORD(l, h)		((WORD)(((BYTE)(l)) | ((WORD)((BYTE)(h))) << 8))
#define		MAKELONG(l, h)		((LONG)(((WORD)(l)) | ((DWORD)((WORD)(h))) << 16))
#define		MAKEDWORD(l, h)		((DWORD)(((WORD)(l)) | ((DWORD)((WORD)(h))) << 16))
#define		LOWORD(l)			((WORD)(l))
#define		HIWORD(l)			((WORD)(((DWORD)(l) >> 16) & 0xFFFF))
#define		LOBYTE(w)			((BYTE)(w))
#define		HIBYTE(w)			((BYTE)(((WORD)(w) >> 8) & 0xFF))
#define		MAKEWPARAM(l, h)	(WPARAM)MAKELONG(l, h)
#define		MAKELPARAM(l, h)	(LPARAM)MAKELONG(l, h)
#define		MAKELRESULT(l, h)	(LRESULT)MAKELONG(l, h)
//====================================================================
//	macro of combine data end
//====================================================================

#endif /*__DEF_DATATYPE_H_2005_09_21__*/
