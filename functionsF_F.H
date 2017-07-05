#ifndef _FUNCTIONSF_F_H_
#define _FUNCTIONSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGETPHIPOINT 
#define GUARDGETPHIPOINT 
// Prototype for Fortran procedure getphipoint ...
//
void FORTRAN_NAME( GETPHIPOINT ,getphipoint )(
      CHFp_REAL(phi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_GETPHIPOINT FORTRAN_NAME( inlineGETPHIPOINT, inlineGETPHIPOINT)
#define FORTNT_GETPHIPOINT FORTRAN_NAME( GETPHIPOINT, getphipoint)

inline void FORTRAN_NAME(inlineGETPHIPOINT, inlineGETPHIPOINT)(
      CHFp_REAL(phi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_GETPHIPOINT");
 FORTRAN_NAME( GETPHIPOINT ,getphipoint )(
      CHFt_REAL(phi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDGETPHIPOINT 

#ifndef GUARDGETLOFPHI 
#define GUARDGETLOFPHI 
// Prototype for Fortran procedure getlofphi ...
//
void FORTRAN_NAME( GETLOFPHI ,getlofphi )(
      CHFp_FRA1(lofphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_CONST_REAL(aCoef)
      ,CHFp_CONST_REAL(bCoef)
      ,CHFp_BOX(box) );

#define FORT_GETLOFPHI FORTRAN_NAME( inlineGETLOFPHI, inlineGETLOFPHI)
#define FORTNT_GETLOFPHI FORTRAN_NAME( GETLOFPHI, getlofphi)

inline void FORTRAN_NAME(inlineGETLOFPHI, inlineGETLOFPHI)(
      CHFp_FRA1(lofphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_CONST_REAL(aCoef)
      ,CHFp_CONST_REAL(bCoef)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETLOFPHI");
 FORTRAN_NAME( GETLOFPHI ,getlofphi )(
      CHFt_FRA1(lofphi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REALVECT(probhi)
      ,CHFt_CONST_REAL(aCoef)
      ,CHFt_CONST_REAL(bCoef)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETLOFPHI 

#ifndef GUARDGETLOFPHIPOINT 
#define GUARDGETLOFPHIPOINT 
// Prototype for Fortran procedure getlofphipoint ...
//
void FORTRAN_NAME( GETLOFPHIPOINT ,getlofphipoint )(
      CHFp_REAL(lofphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(aCoefmult)
      ,CHFp_CONST_REAL(bCoefmult) );

#define FORT_GETLOFPHIPOINT FORTRAN_NAME( inlineGETLOFPHIPOINT, inlineGETLOFPHIPOINT)
#define FORTNT_GETLOFPHIPOINT FORTRAN_NAME( GETLOFPHIPOINT, getlofphipoint)

inline void FORTRAN_NAME(inlineGETLOFPHIPOINT, inlineGETLOFPHIPOINT)(
      CHFp_REAL(lofphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(aCoefmult)
      ,CHFp_CONST_REAL(bCoefmult) )
{
 CH_TIMELEAF("FORT_GETLOFPHIPOINT");
 FORTRAN_NAME( GETLOFPHIPOINT ,getlofphipoint )(
      CHFt_REAL(lofphi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(aCoefmult)
      ,CHFt_CONST_REAL(bCoefmult) );
}
#endif  // GUARDGETLOFPHIPOINT 

#ifndef GUARDGETGRADPHIPOINT 
#define GUARDGETGRADPHIPOINT 
// Prototype for Fortran procedure getgradphipoint ...
//
void FORTRAN_NAME( GETGRADPHIPOINT ,getgradphipoint )(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_GETGRADPHIPOINT FORTRAN_NAME( inlineGETGRADPHIPOINT, inlineGETGRADPHIPOINT)
#define FORTNT_GETGRADPHIPOINT FORTRAN_NAME( GETGRADPHIPOINT, getgradphipoint)

inline void FORTRAN_NAME(inlineGETGRADPHIPOINT, inlineGETGRADPHIPOINT)(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_GETGRADPHIPOINT");
 FORTRAN_NAME( GETGRADPHIPOINT ,getgradphipoint )(
      CHFt_REALVECT(gradphi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDGETGRADPHIPOINT 

}

#endif
