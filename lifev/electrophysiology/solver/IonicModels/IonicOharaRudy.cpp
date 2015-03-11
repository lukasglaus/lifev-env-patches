//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
 @file IonicOharaRudy
 @brief Ionic model O'hara Rudy
 
 @date 03-2015
 @author Thomas Kummer <kummerth@ethz.ch>
 
 @contributors
 @mantainer Thomas Kummer <kummerth@ethz.ch>
 @last update 03-2015
 */


#include <lifev/electrophysiology/solver/IonicModels/IonicOharaRudy.hpp>


namespace LifeV
{

// ===================================================
//! Constructors
// ===================================================
IonicOharaRudy::IonicOharaRudy()  :
    super       ( 41 , 31 ),
    M_nao       ( 140.0 ),
    M_cao       ( 1.8 ),
    M_ko        ( 5.4 ),
    M_BSRmax    ( 0.047 ),
    M_KmBSR     ( 0.00087 ),
    M_BSLmax    ( 1.124 ),
    M_KmBSL     ( 0.0087 ),
    M_cmdnmax   ( 0.05 ),
    M_kmcmdn    ( 0.00238 ),
    M_trpnmax   ( 0.07 ),
    M_kmtrpn    ( 0.0005 ),
    M_csqnmax   ( 10.0 ),
    M_kmcsqn    ( 0.8 ),
    M_aCaMK     ( 0.05 ),
    M_bCaMK     ( 0.00068 ),
    M_CaMKo     ( 0.05 ),
    M_KmCaM     ( 0.0015 ),
    M_KmCaMK    ( 0.15 ),
    M_R         ( 8314.0 ),
    M_T         ( 310.0 ),
    M_F         ( 96485.0 ),
    M_L         ( 0.01 ),
    M_rad       ( 0.0011 ),
    M_vcell     ( 3.7994e-5 ),
    M_Ageo      ( 7.66788e-5 ),
    M_Acap      ( 1.533576e-4 ),
    M_vmyo      ( 0.68*3.7994e-5 ),
    M_vmito     ( 0.26*3.7994e-5 ),
    M_vsr       ( 0.06*3.7994e-5 ),
    M_vnsr      ( 0.0552*3.7994e-5 ),
    M_vjsr      ( 0.0048*3.7994e-5 ),
    M_vss       ( 0.02*3.7994e-5 ),
    M_APD_flag  ( 0 ),
    M_vo        ( -87.5 ),
    M_dt        ( 0.005 ),
    M_t0        ( 0 ),
    M_t         ( 0 ),
    M_vdot      ( 0 ),
    M_p         ( 1 ),
    M_n         ( 0 ),
    M_count     ( 1 )
{
    M_restingConditions.at (0) = -87.5;                        // V
    M_restingConditions.at (1) = 7;                            // nai
    M_restingConditions.at (2) = M_restingConditions.at (1);   // nass
    M_restingConditions.at (3) = 145;                          // ki
    M_restingConditions.at (4) = M_restingConditions.at (3);   // kss
    M_restingConditions.at (5) = 1.0e-4;                       // cai
    M_restingConditions.at (6) = M_restingConditions.at (5);   // cass
    M_restingConditions.at (7) = 1.2;                          // cans
    M_restingConditions.at (8) = M_restingConditions.at (7);   // cajsr
    M_restingConditions.at (9) = 0;                            // m
    M_restingConditions.at (10) = 1;                           // hf
    M_restingConditions.at (11) = 1;                           // hs
    M_restingConditions.at (12) = 1;                           // j
    M_restingConditions.at (13) = 1;                           // hsp
    M_restingConditions.at (14) = 1;                           // jp
    M_restingConditions.at (15) = 0;                           // mL
    M_restingConditions.at (16) = 1;                           // hL
    M_restingConditions.at (17) = 1;                           // hLp
    M_restingConditions.at (18) = 0;                           // a
    M_restingConditions.at (19) = 1;                           // iF
    M_restingConditions.at (20) = 1;                           // iS
    M_restingConditions.at (21) = 0;                           // ap
    M_restingConditions.at (22) = 1;                           // iFp
    M_restingConditions.at (23) = 1;                           // iSp
    M_restingConditions.at (24) = 0;                           // d
    M_restingConditions.at (25) = 1;                           // ff
    M_restingConditions.at (26) = 1;                           // fs
    M_restingConditions.at (27) = 1;                           // fcaf
    M_restingConditions.at (28) = 1;                           // fcas
    M_restingConditions.at (29) = 1;                           // jca
    M_restingConditions.at (30) = 0;                           // nca
    M_restingConditions.at (31) = 1;                           // ffp
    M_restingConditions.at (32) = 1;                           // fcafp
    M_restingConditions.at (33) = 0;                           // xrf
    M_restingConditions.at (34) = 0;                           // xrs
    M_restingConditions.at (35) = 0;                           // xs1
    M_restingConditions.at (36) = 0;                           // xs2
    M_restingConditions.at (37) = 1;                           // xk1
    M_restingConditions.at (38) = 0;                           // Jrelnp
    M_restingConditions.at (39) = 0;                           // Jrelp
    M_restingConditions.at (40) = 0;                           // CaMKt
}

IonicOharaRudy::IonicOharaRudy ( Teuchos::ParameterList& parameterList     )   :
    super       ( 41, 31 )
{

}

IonicOharaRudy::IonicOharaRudy ( const IonicOharaRudy& model )
{

}

// ===================================================
//! Operator
// ===================================================
IonicOharaRudy& IonicOharaRudy::operator= ( const IonicOharaRudy& model )
{

}

// ===================================================
//! Methods
// ===================================================
void IonicOharaRudy::computeGatingRhs ( const   std::vector<Real>&  v,
                                       std::vector<Real>& rhs )
{
    Real V = v[0];
    Real nai = v[1];
    Real nass = v[2];
    Real ki = v[3];
    Real kss = v[4];
    Real cai = v[5];
    Real cass = v[6];
    Real cansr = v[7];
    Real cajsr = v[8];
    Real m = v[9];
    Real hf = v[10];
    Real hs = v[11];
    Real j = v[12];
    Real hsp = v[13];
    Real jp = v[14];
    Real mL = v[15];
    Real hL = v[16];
    Real hLp = v[17];
    Real a = v[18];
    Real iF = v[19];
    Real iS = v[20];
    Real ap = v[21];
    Real iFp = v[22];
    Real iSp = v[23];
    Real d = v[24];
    Real ff = v[25];
    Real fs = v[26];
    Real fcaf = v[27];
    Real fcas = v[28];
    Real jca = v[29];
    Real nca = v[30];
    Real ffp = v[31];
    Real fcafp = v[32];
    Real xrf = v[33];
    Real xrs = v[34];
    Real xs1 = v[35];
    Real xs2 = v[35];
    Real xk1 = v[37];
    Real Jrelnp = v[37];
    Real Jrelp = v[39];
    Real CaMKt = v[40];


    //m
    rhs[0] = dm (V, m);
    //h
    rhs[1] = dh (V, h);
    //j
    rhs[2] = dj (V, j);
    //d
    rhs[3] = dd (V, d);
    //f
    rhs[4] = df (V, f);
    //X
    rhs[5] = dX (V, X);
    //Ca
    rhs[6] = dCa (V, d, f, Ca);
}

void IonicOharaRudy::computeNonGatingRhs ( const   std::vector<Real>&  v,
                                          std::vector<Real>& rhs )
{
    Real V = v[0];
    Real d = v[4];
    Real f = v[5];
    Real Ca = v[7];

    //Ca
    rhs[0] = dCa (V, d, f, Ca);
}

void IonicOharaRudy::computeRhs ( const   std::vector<Real>&  v,
                                 std::vector<Real>& rhs )
{
    Real V = v[0];
    Real nai = v[1];
    Real nass = v[2];
    Real ki = v[3];
    Real kss = v[4];
    Real cai = v[5];
    Real cass = v[6];
    Real cansr = v[7];
    Real cajsr = v[8];
    Real m = v[9];
    Real hf = v[10];
    Real hs = v[11];
    Real j = v[12];
    Real hsp = v[13];
    Real jp = v[14];
    Real mL = v[15];
    Real hL = v[16];
    Real hLp = v[17];
    Real a = v[18];
    Real iF = v[19];
    Real iS = v[20];
    Real ap = v[21];
    Real iFp = v[22];
    Real iSp = v[23];
    Real d = v[24];
    Real ff = v[25];
    Real fs = v[26];
    Real fcaf = v[27];
    Real fcas = v[28];
    Real jca = v[29];
    Real nca = v[30];
    Real ffp = v[31];
    Real fcafp = v[32];
    Real xrf = v[33];
    Real xrs = v[34];
    Real xs1 = v[35];
    Real xs2 = v[35];
    Real xk1 = v[37];
    Real Jrelnp = v[37];
    Real Jrelp = v[39];
    Real CaMKt = v[40];

    revpots(nai, ki);

    //V
    rhs[0] = - Itot (V, m , h , j , d, f, X, Ca);
    //m
    rhs[1] = dm (V, m);
    //h
    rhs[2] = dh (V, h);
    //j
    rhs[3] = dj (V, j);
    //d
    rhs[4] = dd (V, d);
    //f
    rhs[5] = df (V, f);
    //X
    rhs[6] = dX (V, X);
    //Ca
    rhs[7] = dCa (V, d, f, Ca);
}

void IonicOharaRudy::computeGatingVariablesWithRushLarsen ( std::vector<Real>& v, const Real dt )
{
    /*
    Real V = v[0];
    Real m = v[1];
    Real h = v[2];
    Real j = v[3];
    Real d = v[4];
    Real f = v[5];
    Real X = v[6];

    v[1] = minf (V) - ( minf (V) - m ) * std::exp (- dt / tm (V) );
    v[2] = hinf (V) - ( hinf (V) - h ) * std::exp (- dt / th (V) );
    v[3] = jinf (V) - ( jinf (V) - j ) * std::exp (- dt / tj (V) );
    v[4] = dinf (V) - ( dinf (V) - d ) * std::exp (- dt / td (V) );
    v[5] = finf (V) - ( finf (V) - f ) * std::exp (- dt / tf (V) );
    v[6] = Xinf (V) - ( Xinf (V) - X ) * std::exp (- dt / tX (V) );
    */
}

Real IonicOharaRudy::computeLocalPotentialRhs ( const std::vector<Real>& v )
{
    Real dPotential (0.0);

    Real V = v[0];
    Real m = v[1];
    Real h = v[2];
    Real j = v[3];
    Real d = v[4];
    Real f = v[5];
    Real X = v[6];
    Real Ca = v[7];

    dPotential = - Itot (V, m , h , j , d, f, X, Ca);

    return dPotential;
}


void IonicOharaRudy::showMe()
{
    std::cout << "\n\n************************************";
    std::cout << "\n\tHi, I'm the Luo Rudy Phase I model";

    std::cout << "\nENa: " << M_ENa;
    std::cout << "\ngNa: " << M_gNa;
    std::cout << "\ngsi: " << M_gsi;
    std::cout << "\nK0: " << M_K0;
    std::cout << "\nEK: " << M_EK;
    std::cout << "\nEK1: " << M_EK1;
    std::cout << "\nEKp: " << M_EKp;
    std::cout << "\ngKp: " << M_gKp;
    std::cout << "\ngb: " << M_gb;
    std::cout << "\n************************************\n\n";
}

    
void IonicOharaRudy::revpots ( Real nai,
                               Real ki )
{
    M_ENa=(M_R*M_T/M_F)*log(M_nao/nai);
    M_EK=(M_R*M_T/M_F)*log(M_cko/ki);
    M_EKs=(M_R*M_T/M_F)*log((M_ko+0.01833*M_nao)/(ki+0.01833*nai));
}

    
Real IonicOharaRudy::dCaMKt ( Real CaMKt,
                              Real cass )
{
    M_CaMKb=M_CaMKo*(1.0-CaMKt)/(1.0+M_KmCaM/cass);
    M_CaMKa=M_CaMKb+CaMKt;
    return (aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt);
}
    
    


// updateConstants
    M_JdiffNa=(nass-nai)/2.0;
    M_JdiffK=(kss-ki)/2.0;
    M_Jdiff=(cass-cai)/0.2;

// dJrelnp
    double bt=4.75;
    double a_rel=0.5*bt;
    double Jrel_inf=a_rel*(-M_ICaL)/(1.0+pow(1.5/cajsr,8.0));


    if (celltype==2) // ?
    {
        Jrel_inf*=1.7;
    }

    double tau_rel=bt/(1.0+0.0123/cajsr);


    if (tau_rel<0.005)
    {
        tau_rel=0.005;
    }

    Jrelnp=Jrel_inf-(Jrel_inf-Jrelnp)*exp(-dt/tau_rel); // ?


// dJrelp
    double btp=1.25*bt;
    double a_relp=0.5*btp;
    double Jrel_infp=a_relp*(-M_ICaL)/(1.0+pow(1.5/cajsr,8.0));

    if (celltype==2) // ?
    {
        Jrel_infp*=1.7;
    }

    double tau_relp=btp/(1.0+0.0123/cajsr);
    if (tau_relp<0.005)
    {
        tau_relp=0.005;
    }

    Jrelp=Jrel_infp-(Jrel_infp-Jrelp)*exp(-dt/tau_relp); // ?


// updateConstants
    double fJrelp=(1.0/(1.0+KmCaMK/CaMKa));

    Jrel=(1.0-fJrelp)*Jrelnp+fJrelp*Jrelp; // ?
    
    double Jupnp=0.004375*cai/(cai+0.00092);
    double Jupp=2.75*0.004375*cai/(cai+0.00092-0.00017);

    if (celltype==1) // ?
    {
        Jupnp*=1.3;
        Jupp*=1.3;
    }

    double fJupp=(1.0/(1.0+KmCaMK/CaMKa));
    M_Jleak=0.0039375*cansr/15.0;
    Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-M_Jleak;
    Jtr=(cansr-cajsr)/100.0;


// dnai
    nai+=dt*(-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo);

// dnass
    nass+=dt*(-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa);

// dki
    ki+=dt*(-(Ito+IKr+IKs+IK1+IKb+Ist-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo);

// dkss
    kss+=dt*(-(ICaK)*Acap/(F*vss)-JdiffK);



// dcai
    double Bcai;
    if (celltype==1)
    {
        Bcai=1.0/(1.0+1.3*cmdnmax*kmcmdn/pow(kmcmdn+cai,2.0)+trpnmax*kmtrpn/pow(kmtrpn+cai,2.0));
    }
    else
    {
        Bcai=1.0/(1.0+cmdnmax*kmcmdn/pow(kmcmdn+cai,2.0)+trpnmax*kmtrpn/pow(kmtrpn+cai,2.0));
    }
    cai+=dt*(Bcai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo));



// dcass
Real dcass
    double Bcass=1.0/(1.0+M_BSRmax*M_KmBSR/pow(M_KmBSR+cass,2.0)+BSLmax*KmBSL/pow(KmBSL+cass,2.0));
    cass+=dt*(Bcass*(-(ICaL-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff));


Real dcansr IonicOharaRudy::dcansr ()
{
    return (M_Jup-M_Jtr*M_vjsr/M_vnsr);
}

    
Real dcajsr IonicOharaRudy::dcajsr ( Real cajsr )
{
    double Bcajsr=1.0/(1.0+M_csqnmax*M_kmcsqn/pow(M_kmcsqn+cajsr,2.0));
    return (Bcajsr*(M_Jtr-M_Jrel));
}

}