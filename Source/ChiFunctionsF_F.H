#ifndef _CHIFUNCTIONSF_F_H_
#define _CHIFUNCTIONSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGETCHIPOINT 
#define GUARDGETCHIPOINT 
// Prototype for Fortran procedure getchipoint ...
//
void FORTRAN_NAME( GETCHIPOINT ,getchipoint )(
      CHFp_REAL(chi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_GETCHIPOINT FORTRAN_NAME( inlineGETCHIPOINT, inlineGETCHIPOINT)
#define FORTNT_GETCHIPOINT FORTRAN_NAME( GETCHIPOINT, getchipoint)

inline void FORTRAN_NAME(inlineGETCHIPOINT, inlineGETCHIPOINT)(
      CHFp_REAL(chi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_GETCHIPOINT");
 FORTRAN_NAME( GETCHIPOINT ,getchipoint )(
      CHFt_REAL(chi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDGETCHIPOINT 

#ifndef GUARDGETLOFCHI 
#define GUARDGETLOFCHI 
// Prototype for Fortran procedure getlofchi ...
//
void FORTRAN_NAME( GETLOFCHI ,getlofchi )(
      CHFp_FRA1(lofchi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_CONST_REAL(aCoef)
      ,CHFp_CONST_REAL(bCoef)
      ,CHFp_CONST_REAL(cCoef)
      ,CHFp_BOX(box) );

#define FORT_GETLOFCHI FORTRAN_NAME( inlineGETLOFCHI, inlineGETLOFCHI)
#define FORTNT_GETLOFCHI FORTRAN_NAME( GETLOFCHI, getlofchi)

inline void FORTRAN_NAME(inlineGETLOFCHI, inlineGETLOFCHI)(
      CHFp_FRA1(lofchi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_CONST_REAL(aCoef)
      ,CHFp_CONST_REAL(bCoef)
      ,CHFp_CONST_REAL(cCoef)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETLOFCHI");
 FORTRAN_NAME( GETLOFCHI ,getlofchi )(
      CHFt_FRA1(lofchi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REALVECT(probhi)
      ,CHFt_CONST_REAL(aCoef)
      ,CHFt_CONST_REAL(bCoef)
      ,CHFt_CONST_REAL(cCoef)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETLOFCHI 

#ifndef GUARDGETLOFCHIPOINT 
#define GUARDGETLOFCHIPOINT 
// Prototype for Fortran procedure getlofchipoint ...
//
void FORTRAN_NAME( GETLOFCHIPOINT ,getlofchipoint )(
      CHFp_REAL(lofchi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(aCoefmult)
      ,CHFp_CONST_REAL(bCoefmult)
      ,CHFp_CONST_REAL(cCoefmult) );

#define FORT_GETLOFCHIPOINT FORTRAN_NAME( inlineGETLOFCHIPOINT, inlineGETLOFCHIPOINT)
#define FORTNT_GETLOFCHIPOINT FORTRAN_NAME( GETLOFCHIPOINT, getlofchipoint)

inline void FORTRAN_NAME(inlineGETLOFCHIPOINT, inlineGETLOFCHIPOINT)(
      CHFp_REAL(lofchi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(aCoefmult)
      ,CHFp_CONST_REAL(bCoefmult)
      ,CHFp_CONST_REAL(cCoefmult) )
{
 CH_TIMELEAF("FORT_GETLOFCHIPOINT");
 FORTRAN_NAME( GETLOFCHIPOINT ,getlofchipoint )(
      CHFt_REAL(lofchi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(aCoefmult)
      ,CHFt_CONST_REAL(bCoefmult)
      ,CHFt_CONST_REAL(cCoefmult) );
}
#endif  // GUARDGETLOFCHIPOINT 

#ifndef GUARDGETGRADCHIPOINT 
#define GUARDGETGRADCHIPOINT 
// Prototype for Fortran procedure getgradchipoint ...
//
void FORTRAN_NAME( GETGRADCHIPOINT ,getgradchipoint )(
      CHFp_REALVECT(gradchi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_GETGRADCHIPOINT FORTRAN_NAME( inlineGETGRADCHIPOINT, inlineGETGRADCHIPOINT)
#define FORTNT_GETGRADCHIPOINT FORTRAN_NAME( GETGRADCHIPOINT, getgradchipoint)

inline void FORTRAN_NAME(inlineGETGRADCHIPOINT, inlineGETGRADCHIPOINT)(
      CHFp_REALVECT(gradchi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_GETGRADCHIPOINT");
 FORTRAN_NAME( GETGRADCHIPOINT ,getgradchipoint )(
      CHFt_REALVECT(gradchi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDGETGRADCHIPOINT 

}

#endif
