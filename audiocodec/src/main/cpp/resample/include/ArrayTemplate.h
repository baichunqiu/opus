// ArrayTemplate.h
//
//包括5个模板：
//		1.CArray1D		——	基本的1D矩阵，具有赋值、转置、矩阵Copy等基本功能，适用于所有数据类型
//								如果使用复杂数据结构，必须重载EmptyAry()以防内存泄漏
//
//*************************修改记录*********************************//
//
//		版本		修改者		修改日期			备注
//
//		1.0			李绍恩		2002-11-30			完成二维运算操作
//		2.0			李绍恩		2002-12-7			增加一维运算操作和二维的[][]功能实现
//
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TEMPLATE_ARRAY2D_INCLUDED_)
#define AFX_TEMPLATE_ARRAY2D_INCLUDED_

#if 0
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#endif

#include "Def_DataType.h"

#define	ASSERT(x)

extern "C" void* memset(void*, int, UINT);
extern "C" void* memcpy(void*, const void*, UINT);

//////////**************** CArray1D ***************///////////////////
//基本的1D矩阵，具有赋值、转置、矩阵Copy等基本功能，适用于所有数据类型
//如果使用复杂数据结构，必须重载EmptyAry()以防内存泄漏
//////////*****************************************///////////////////

template <typename TYPE> class CArray1D 
{
public:
//类型定义
	typedef	TYPE		value_type;
	typedef	TYPE*		iterator;
	typedef	const TYPE*	const_iterator;
	typedef	TYPE&		reference;
	typedef	const TYPE&	const_reference;
//	typedef ptrdiff_t	difference_type;

//构造、析构函数
	CArray1D();								//标准构造函数
	CArray1D(const CArray1D <TYPE>& oda);	//拷贝构造函数
	CArray1D(const UINT nColCount);			//初始化构造函数

	virtual ~CArray1D() {EmptyAry();}				//析构函数

//操作
	//返回矩阵数据的内存指针
	inline TYPE* GetDataBuf() const	{ return m_lpData; }
	
	//返回矩阵数据的尺寸
	inline int GetDataTypeSize() const { return sizeof(TYPE); }
	
	//从指定内存处复制数据
	inline void SetDataBuf(const TYPE* lpDataBuf)	{ ASSERT(lpDataBuf); ASSERT(m_lpData); memcpy(m_lpData, lpDataBuf, m_nColCount * sizeof(TYPE)); }
	
	//返回矩阵的空间大小
	inline const UINT GetLen() const {return m_nColCount;}

	//复制矩阵
	BOOL Copy(const CArray1D <TYPE>& oda);

	//安全复制矩阵（不受元素类型的限制）
	BOOL SafeCopy(const CArray1D <TYPE>& oda);

	//将所有数据重设为指定值
	void Reset(const TYPE& typeValue = TYPE());

	//交换两个矩阵的数据
	BOOL Swap(CArray1D <TYPE>& oda);

//运算符重载
	//重载运算符 []
	inline TYPE& operator [](const UINT nIndex);
	
	//重载运算符 =
	inline const CArray1D <TYPE>& operator = (const CArray1D <TYPE>& oda);

//清空数据
	//释放存储空间
	virtual void EmptyAry();

//初始化函数
	//初始化数据存储空间
	BOOL Init(const UINT nColCount);

protected:
	//数据存储空间指针
	TYPE*	m_lpData;
	
	//数据存储空间的大小
	UINT	m_nColCount;
};

////////////////////////////////////////////////////////////////////////
//函数名称:CArray1D()
//参数意义:
// 返回值 :
// 注  释 :标准构造函数
////////////////////////////////////////////////////////////////////////
template <typename TYPE>
CArray1D <TYPE>::CArray1D()
{
	m_lpData = NULL;
	m_nColCount = 0;
}

////////////////////////////////////////////////////////////////////////
//函数名称:CArray1D <TYPE>::CArray1D(const CArray1D <TYPE>& oda)
//参数意义:CArray1D <TYPE>& oda		- 常值一维矩阵
// 返回值 :
// 注  释 :拷贝构造函数,从另一个一维矩阵复制数据
////////////////////////////////////////////////////////////////////////
template <typename TYPE>
CArray1D <TYPE>::CArray1D(const CArray1D <TYPE>& oda)
{
	m_lpData = NULL;
	m_nColCount = 0;

	if(!Copy(oda))
	{
	#ifdef _AFXDLL
		AfxGetMainWnd()->MessageBox("一维数组Copy失败！", "CArray1D", MB_OK | MB_ICONSTOP);
	#endif
	}
}

////////////////////////////////////////////////////////////////////////
//函数名称:CArray1D <TYPE>::CArray1D(const UINT nColCount)
//参数意义:UINT nColCount			- 矩阵元素的个数
// 返回值 :
// 注  释 :初始化构造函数，初始化时指定矩阵存储空间的大小
////////////////////////////////////////////////////////////////////////
template <typename TYPE>
CArray1D <TYPE>::CArray1D(const UINT nColCount)
{
	m_lpData = NULL;
	m_nColCount = 0;

	if(!Init(nColCount))
	{
	#ifdef _AFXDLL
		AfxGetMainWnd()->MessageBox("一维数组初始化失败！", "CArray1D", MB_OK | MB_ICONSTOP);
	#endif
	}
}

////////////////////////////////////////////////////////////////////////
//函数名称:BOOL CArray1D <TYPE>::Init(const UINT nColCount)
//参数意义:UINT nColCount			- 矩阵元素的个数
// 返回值 :如果成功，返回TRUE；否则返回 FALSE
// 注  释 :初始化矩阵存储空间的大小，返回初始化是否成功
////////////////////////////////////////////////////////////////////////
template <typename TYPE>
BOOL CArray1D <TYPE>::Init(const UINT nColCount)
{
	if((m_lpData != NULL) && (m_nColCount == nColCount))
		return TRUE;

	EmptyAry();

#ifdef _AFXDLL
	//TRY CATCH STATEMENT
	try
	{
#endif
		m_lpData = new TYPE[nColCount];
#ifdef _AFXDLL
	}
	catch(CMemoryException* e)
	{
		char chErr[256];
		memset(chErr, 0, 256);
		e->GetErrorMessage(chErr, 255);
		e->Delete();
		AfxGetMainWnd()->MessageBox(chErr, "CArray1D", MB_OK | MB_ICONSTOP);
	}
	//CATCH END
#endif

	if(!m_lpData)
		return FALSE;

	m_nColCount = nColCount;

	return TRUE;
}

////////////////////////////////////////////////////////////////////////
//函数名称:void CArray1D <TYPE>::EmptyAry()
//参数意义:
// 返回值 :
// 注  释 :释放存储空间
////////////////////////////////////////////////////////////////////////
template <typename TYPE>
void CArray1D <TYPE>::EmptyAry()
{
	if(m_lpData)
	{
		delete [] m_lpData;
		m_lpData = NULL;
	}

	m_nColCount = 0;
}

////////////////////////////////////////////////////////////////////////
//函数名称:BOOL CArray1D <TYPE>::Copy(const CArray1D <TYPE>& oda)
//参数意义:CArray1D <TYPE>& oda		- 常值一维矩阵
// 返回值 :如果成功，返回TRUE；否则返回 FALSE
// 注  释 :从另一个一维矩阵复制数据，存储空间随之改变
////////////////////////////////////////////////////////////////////////
template <typename TYPE>
BOOL CArray1D <TYPE>::Copy(const CArray1D <TYPE>& oda)
{
	if(this == &oda)
		return TRUE;
	
	if(!Init(oda.GetLen()))
		return FALSE;

	SetDataBuf(oda.GetDataBuf());

	return TRUE;
}

////////////////////////////////////////////////////////////////////////
//函数名称:BOOL CArray1D <TYPE>::SafeCopy(const CArray1D <TYPE>& oda)
//参数意义:CArray1D <TYPE>& oda		- 常值一维矩阵
// 返回值 :如果成功，返回TRUE；否则返回 FALSE
// 注  释 :从另一个一维矩阵复制数据，存储空间随之改变，不受元素类型的限制
////////////////////////////////////////////////////////////////////////
template <typename TYPE>
BOOL CArray1D <TYPE>::SafeCopy(const CArray1D <TYPE>& oda)
{
	if(this == &oda)
		return TRUE;
	
	if(!Init(oda.GetLen()))
		return FALSE;

	for(UINT i=0; i<GetLen(); i++)
		(*this)[i]	= oda[i];

	return TRUE;
}

////////////////////////////////////////////////////////////////////////
//函数名称:void CArray1D <TYPE>::Reset(const TYPE& typeValue = TYPE())
//参数意义:TYPE& typeValue		- 新值
// 返回值 :
// 注  释 :将所有数据重设为指定值
////////////////////////////////////////////////////////////////////////
template <typename TYPE>
void CArray1D <TYPE>::Reset(const TYPE& typeValue)
{
	ASSERT(m_lpData);

	if((typeValue == 0) || (sizeof(TYPE) == 1))
		memset(m_lpData, (int)typeValue, sizeof(TYPE)*GetLen());
	else
	{
		for(UINT i=0; i<GetLen(); i++)
			(*this)[i] = typeValue;
	}
}

////////////////////////////////////////////////////////////////////////
//函数名称:BOOL CArray1D <TYPE>::Swap(CArray1D <TYPE>& oda)
//参数意义:CArray1D <TYPE>& oda		- 一维矩阵
// 返回值 :如果成功，返回TRUE；否则返回 FALSE
// 注  释 :交换两个矩阵的数据
////////////////////////////////////////////////////////////////////////
template <typename TYPE>
BOOL CArray1D <TYPE>::Swap(CArray1D <TYPE>& oda)
{
	if(this == &oda)
		return TRUE;

	TYPE* lpData	= oda.GetDataBuf();
	UINT  nLen		= oda.GetLen();

	oda.m_lpData	= m_lpData;
	oda.m_nColCount	= m_nColCount;

	m_lpData	= lpData;
	m_nColCount	= nLen;

	return TRUE;
}

////////////////////////////////////////////////////////////////////////
//函数名称:TYPE& CArray1D <TYPE>::operator [](const UINT nIndex)
//参数意义:UINT nIndex			-  0Base索引
// 返回值 :TYPE 引用
// 注  释 :重载操作符[]，返回指定内存的引用，可直接操作指定内存
////////////////////////////////////////////////////////////////////////
template <typename TYPE>
TYPE& CArray1D <TYPE>::operator [](const UINT nIndex)
{
	ASSERT(nIndex < m_nColCount);

	return m_lpData[nIndex];
}

////////////////////////////////////////////////////////////////////////
//函数名称:const CArray1D <TYPE>& CArray1D <TYPE>::operator = (const CArray1D <TYPE>& oda)
//参数意义:CArray1D <TYPE>& oda		- 常值一维矩阵
// 返回值 :返回赋值结果
// 注  释 :重载操作符=，作用同Copy
////////////////////////////////////////////////////////////////////////
template <typename TYPE>
const CArray1D <TYPE>& CArray1D <TYPE>::operator = (const CArray1D <TYPE>& oda)
{
	if(!Copy(oda))
	{
	#ifdef _AFXDLL
		AfxGetMainWnd()->MessageBox("一维数组赋值失败！", "CArray1D", MB_OK | MB_ICONSTOP);
	#endif
	}

	return *this;
}

#endif // !defined(AFX_TEMPLATE_ARRAY2D_INCLUDED_)
