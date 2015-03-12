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
    Real Jrelnp = v[38];
    Real Jrelp = v[39];
    Real CaMKt = v[40];

    revpots ( nai, ki );
    updateConstants ( nass, nai, kss, ki, cass, cai, Jrelnp, Jrelp );
                     
    rhs[0] =
    rhs[1] = dnai ();
    rhs[2] = dnass ();
    rhs[3] = dki ();
    rhs[4] = dkss ();
    rhs[5] = dcai ( cai );
    rhs[6] = dcass ( cass );
    rhs[7] = dcansr ();
    rhs[8] = dcajsr ();
    rhs[9] =
    rhs[10] =
    rhs[11] =
    rhs[12] =
    rhs[13] =
    rhs[14] =
    rhs[15] =
    rhs[16] =
    rhs[17] =
    rhs[18] =
    rhs[19] =
    rhs[20] =
    rhs[21] =
    rhs[22] =
    rhs[23] =
    rhs[24] =
    rhs[25] =
    rhs[26] =
    rhs[27] =
    rhs[28] =
    rhs[29] =
    rhs[30] =
    rhs[31] =
    rhs[32] =
    rhs[33] =
    rhs[34] =
    rhs[35] =
    rhs[36] =
    rhs[37] =
    rhs[38] = dJrelnp ( Jrelnp, cajsr );
    rhs[39] = dJrelp ( Jrelp, cajsr );
    rhs[40] = dCaMKt ( CaMKt, cass );
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
    

void IonicOharaRudy::updateConstants( Real nass,
                                      Real nai,
                                      Real kss,
                                      Real ki,
                                      Real cass,
                                      Real cai,
                                      Real Jrelnp,
                                      Real Jrelp)
{
    M_JdiffNa=(nass-nai)/2.0;
    M_JdiffK=(kss-ki)/2.0;
    M_Jdiff=(cass-cai)/0.2;
    
    double fJrelp=(1.0/(1.0+M_KmCaMK/M_CaMKa));
    
    M_Jrel=(1.0-fJrelp)*Jrelnp+fJrelp*Jrelp;
    
    double Jupnp=0.004375*cai/(cai+0.00092);
    double Jupp=2.75*0.004375*cai/(cai+0.00092-0.00017);
    
    //if (celltype==1) // ?
    //{
     //   Jupnp*=1.3;
       // Jupp*=1.3;
    //}
    
    double fJupp=(1.0/(1.0+M_KmCaMK/M_CaMKa));
    M_Jleak=0.0039375*cansr/15.0;
    M_Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-M_Jleak;
    M_Jtr=(cansr-cajsr)/100.0;
}
    

Real IonicOharaRudy::dJrelnp ( Real Jrelnp,
                               Real cajsr )
{
    double bt=4.75;
    double a_rel=0.5*bt;
    double Jrel_inf=a_rel*(-M_ICaL)/(1.0+pow(1.5/cajsr,8.0));
    Jrel_inf*=1.7; // ? if celltype ==2 ?
    
    double tau_rel=bt/(1.0+0.0123/cajsr);
    if (tau_rel<0.005)
    {
        tau_rel=0.005;
    }

    return (Jrel_inf - Jrelnp) / tau_rel; // Jrelnp=Jrel_inf-(Jrel_inf-Jrelnp)*exp(-dt/tau_rel);
}

Real IonicOharaRudy::dJrelp ( Real Jrelp,
                              Real cajsr )
{
    double btp=1.25*bt;
    double a_relp=0.5*btp;
    double Jrel_infp=a_relp*(-M_ICaL)/(1.0+pow(1.5/cajsr,8.0));
    Jrel_infp*=1.7;  // ? if celltype ==2 ?

    double tau_relp=btp/(1.0+0.0123/cajsr);
    if (tau_relp<0.005)
    {
        tau_relp=0.005;
    }

    return (Jrel_infp - Jrelp) / tau_relp; //Jrelp=Jrel_infp-(Jrel_infp-Jrelp)*exp(-dt/tau_relp);
}

    
Real IonicOharaRudy::dnai ()
{
    return (-(M_INa+M_INaL+3.0*M_INaCa_i+3.0*M_INaK+M_INab)*M_Acap/(M_F*M_vmyo)+M_JdiffNa*M_vss/M_vmyo);
}
        
        
Real IonicOharaRudy::dnass ()
{
    return (-(M_ICaNa+3.0*M_INaCa_ss)*M_Acap/(M_F*M_vss)-M_JdiffNa);
}
            
            
Real IonicOharaRudy::dki ()
{
    return (-(M_Ito+M_IKr+M_IKs+M_IK1+M_IKb+M_Ist-2.0*M_INaK)*M_Acap/(M_F*M_vmyo)+M_JdiffK*M_vss/M_vmyo);
}
                
                
Real IonicOharaRudy::dkss ()
{
    return (-(M_ICaK)*M_Acap/(M_F*M_vss)-M_JdiffK);
}


Real IonicOharaRudy::dcai ( Real cai )
{
    double Bcai;
    //if (celltype==1)
    //{
    //    Bcai=1.0/(1.0+1.3*cmdnmax*kmcmdn/pow(kmcmdn+cai,2.0)+trpnmax*kmtrpn/pow(kmtrpn+cai,2.0));
    //}
    //else
    //{
        Bcai=1.0/(1.0+M_cmdnmax*M_kmcmdn/pow(M_kmcmdn+cai,2.0)+M_trpnmax*M_kmtrpn/pow(M_kmtrpn+cai,2.0));
    //}
    return (Bcai*(-(M_IpCa+M_ICab-2.0*M_INaCa_i)*M_Acap/(2.0*M_F*M_vmyo)-M_Jup*M_vnsr/M_vmyo+M_Jdiff*M_vss/M_vmyo));
}


Real IonicOharaRudy::dcass ( Real cass )
{
    double Bcass=1.0/(1.0+M_BSRmax*M_KmBSR/pow(M_KmBSR+cass,2.0)+M_BSLmax*M_KmBSL/pow(M_KmBSL+cass,2.0));
    return (Bcass*(-(M_ICaL-2.0*M_INaCa_ss)*M_Acap/(2.0*M_F*M_vss)+M_Jrel*M_vjsr/M_vss-M_Jdiff));
}
    

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






// update (CaMKt, cass)
    M_CaMKb=M_CaMKo*(1.0-CaMKt)/(1.0+M_KmCaM/cass);
    M_CaMKa=M_CaMKb+CaMKt;


// where to put them?
    double vffrt=M_v*M_F*M_F/(M_R*M_T);
    double vfrt=V*F/(M_R*M_T);

Real dm ( m )
{
    double mss=1.0/(1.0+exp((-(V+39.57))/9.871));
    double tm=1.0/(6.765*exp((V+11.64)/34.77)+8.552*exp(-(V+77.42)/5.955));
    return ( mss - m ) / tm; //m=mss-(mss-m)*exp(-dt/tm);
}
    double hss=1.0/(1+exp((v+82.90)/6.086));
    double thf=1.0/(1.432e-5*exp(-(v+1.196)/6.285)+6.149*exp((v+0.5096)/20.27));
    double ths=1.0/(0.009794*exp(-(v+17.95)/28.05)+0.3343*exp((v+5.730)/56.66));
    double Ahf=0.99;
    double Ahs=1.0-Ahf;
    hf=hss-(hss-hf)*exp(-dt/thf);
    hs=hss-(hss-hs)*exp(-dt/ths);
    double h=Ahf*hf+Ahs*hs;
    double jss=hss;
    double tj=2.038+1.0/(0.02136*exp(-(v+100.6)/8.281)+0.3052*exp((v+0.9941)/38.45));
    j=jss-(jss-j)*exp(-dt/tj);
    double hssp=1.0/(1+exp((v+89.1)/6.086));
    double thsp=3.0*ths;
    hsp=hssp-(hssp-hsp)*exp(-dt/thsp);
    double hp=Ahf*hf+Ahs*hsp;
    double tjp=1.46*tj;
    jp=jss-(jss-jp)*exp(-dt/tjp);
    double GNa=75;
    double fINap=(1.0/(1.0+KmCaMK/CaMKa));
    INa=GNa*(v-ENa)*m*m*m*((1.0-fINap)*h*j+fINap*hp*jp);
    
    double mLss=1.0/(1.0+exp((-(v+42.85))/5.264));
    double tmL=tm;
    mL=mLss-(mLss-mL)*exp(-dt/tmL);
    double hLss=1.0/(1.0+exp((v+87.61)/7.488));
    double thL=200.0;
    hL=hLss-(hLss-hL)*exp(-dt/thL);
    double hLssp=1.0/(1.0+exp((v+93.81)/7.488));
    double thLp=3.0*thL;
    hLp=hLssp-(hLssp-hLp)*exp(-dt/thLp);
    double GNaL=0.0075;
    if (celltype==1)
    {
        GNaL*=0.6;
    }
    double fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
    INaL=GNaL*(v-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);
    
    double ass=1.0/(1.0+exp((-(v-14.34))/14.82));
    double ta=1.0515/(1.0/(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+3.5/(1.0+exp((v+100.0)/29.3814)));
    a=ass-(ass-a)*exp(-dt/ta);
    double iss=1.0/(1.0+exp((v+43.94)/5.711));
    double delta_epi;
    if (celltype==1)
    {
        delta_epi=1.0-(0.95/(1.0+exp((v+70.0)/5.0)));
    }
    else
    {
        delta_epi=1.0;
    }
    double tiF=4.562+1/(0.3933*exp((-(v+100.0))/100.0)+0.08004*exp((v+50.0)/16.59));
    double tiS=23.62+1/(0.001416*exp((-(v+96.52))/59.05)+1.780e-8*exp((v+114.1)/8.079));
    tiF*=delta_epi;
    tiS*=delta_epi;
    double AiF=1.0/(1.0+exp((v-213.6)/151.2));
    double AiS=1.0-AiF;
    iF=iss-(iss-iF)*exp(-dt/tiF);
    iS=iss-(iss-iS)*exp(-dt/tiS);
    double i=AiF*iF+AiS*iS;
    double assp=1.0/(1.0+exp((-(v-24.34))/14.82));
    ap=assp-(assp-ap)*exp(-dt/ta);
    double dti_develop=1.354+1.0e-4/(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154));
    double dti_recover=1.0-0.5/(1.0+exp((v+70.0)/20.0));
    double tiFp=dti_develop*dti_recover*tiF;
    double tiSp=dti_develop*dti_recover*tiS;
    iFp=iss-(iss-iFp)*exp(-dt/tiFp);
    iSp=iss-(iss-iSp)*exp(-dt/tiSp);
    double ip=AiF*iFp+AiS*iSp;
    double Gto=0.02;
    if (celltype==1)
    {
        Gto*=4.0;
    }
    if (celltype==2)
    {
        Gto*=4.0;
    }
    double fItop=(1.0/(1.0+KmCaMK/CaMKa));
    Ito=Gto*(v-EK)*((1.0-fItop)*a*i+fItop*ap*ip);
    
    double dss=1.0/(1.0+exp((-(v+3.940))/4.230));
    double td=0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));
    d=dss-(dss-d)*exp(-dt/td);
    double fss=1.0/(1.0+exp((v+19.58)/3.696));
    double tff=7.0+1.0/(0.0045*exp(-(v+20.0)/10.0)+0.0045*exp((v+20.0)/10.0));
    double tfs=1000.0+1.0/(0.000035*exp(-(v+5.0)/4.0)+0.000035*exp((v+5.0)/6.0));
    double Aff=0.6;
    double Afs=1.0-Aff;
    ff=fss-(fss-ff)*exp(-dt/tff);
    fs=fss-(fss-fs)*exp(-dt/tfs);
    double f=Aff*ff+Afs*fs;
    double fcass=fss;
    double tfcaf=7.0+1.0/(0.04*exp(-(v-4.0)/7.0)+0.04*exp((v-4.0)/7.0));
    double tfcas=100.0+1.0/(0.00012*exp(-v/3.0)+0.00012*exp(v/7.0));
    double Afcaf=0.3+0.6/(1.0+exp((v-10.0)/10.0));
    double Afcas=1.0-Afcaf;
    fcaf=fcass-(fcass-fcaf)*exp(-dt/tfcaf);
    fcas=fcass-(fcass-fcas)*exp(-dt/tfcas);
    double fca=Afcaf*fcaf+Afcas*fcas;
    double tjca=75.0;
    jca=fcass-(fcass-jca)*exp(-dt/tjca);
    double tffp=2.5*tff;
    ffp=fss-(fss-ffp)*exp(-dt/tffp);
    double fp=Aff*ffp+Afs*fs;
    double tfcafp=2.5*tfcaf;
    fcafp=fcass-(fcass-fcafp)*exp(-dt/tfcafp);
    double fcap=Afcaf*fcafp+Afcas*fcas;
    double Kmn=0.002;
    double k2n=1000.0;
    double km2n=jca*1.0;
    double anca=1.0/(k2n/km2n+pow(1.0+Kmn/cass,4.0));
    nca=anca*k2n/km2n-(anca*k2n/km2n-nca)*exp(-km2n*dt);
    double PhiCaL=4.0*vffrt*(cass*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
    double PhiCaNa=1.0*vffrt*(0.75*nass*exp(1.0*vfrt)-0.75*nao)/(exp(1.0*vfrt)-1.0);
    double PhiCaK=1.0*vffrt*(0.75*kss*exp(1.0*vfrt)-0.75*ko)/(exp(1.0*vfrt)-1.0);
    double zca=2.0;
    double PCa=0.0001;
    if (celltype==1)
    {
        PCa*=1.2;
    }
    if (celltype==2)
    {
        PCa*=2.5;
    }
    double PCap=1.1*PCa;
    double PCaNa=0.00125*PCa;
    double PCaK=3.574e-4*PCa;
    double PCaNap=0.00125*PCap;
    double PCaKp=3.574e-4*PCap;
    double fICaLp=(1.0/(1.0+KmCaMK/CaMKa));
    ICaL=(1.0-fICaLp)*PCa*PhiCaL*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCap*PhiCaL*d*(fp*(1.0-nca)+jca*fcap*nca);
    ICaNa=(1.0-fICaLp)*PCaNa*PhiCaNa*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaNap*PhiCaNa*d*(fp*(1.0-nca)+jca*fcap*nca);
    ICaK=(1.0-fICaLp)*PCaK*PhiCaK*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaKp*PhiCaK*d*(fp*(1.0-nca)+jca*fcap*nca);
    
    double xrss=1.0/(1.0+exp((-(v+8.337))/6.789));
    double txrf=12.98+1.0/(0.3652*exp((v-31.66)/3.869)+4.123e-5*exp((-(v-47.78))/20.38));
    double txrs=1.865+1.0/(0.06629*exp((v-34.70)/7.355)+1.128e-5*exp((-(v-29.74))/25.94));
    double Axrf=1.0/(1.0+exp((v+54.81)/38.21));
    double Axrs=1.0-Axrf;
    xrf=xrss-(xrss-xrf)*exp(-dt/txrf);
    xrs=xrss-(xrss-xrs)*exp(-dt/txrs);
    double xr=Axrf*xrf+Axrs*xrs;
    double rkr=1.0/(1.0+exp((v+55.0)/75.0))*1.0/(1.0+exp((v-10.0)/30.0));
    double GKr=0.046;
    if (celltype==1)
    {
        GKr*=1.3;
    }
    if (celltype==2)
    {
        GKr*=0.8;
    }
    IKr=GKr*sqrt(ko/5.4)*xr*rkr*(v-EK);
    
    double xs1ss=1.0/(1.0+exp((-(v+11.60))/8.932));
    double txs1=817.3+1.0/(2.326e-4*exp((v+48.28)/17.80)+0.001292*exp((-(v+210.0))/230.0));
    xs1=xs1ss-(xs1ss-xs1)*exp(-dt/txs1);
    double xs2ss=xs1ss;
    double txs2=1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp((-(v+66.54))/31.0));
    xs2=xs2ss-(xs2ss-xs2)*exp(-dt/txs2);
    double KsCa=1.0+0.6/(1.0+pow(3.8e-5/cai,1.4));
    double GKs=0.0034;
    if (celltype==1)
    {
        GKs*=1.4;
    }
    IKs=GKs*KsCa*xs1*xs2*(v-EKs);
    
    double xk1ss=1.0/(1.0+exp(-(v+2.5538*ko+144.59)/(1.5692*ko+3.8115)));
    double txk1=122.2/(exp((-(v+127.2))/20.36)+exp((v+236.8)/69.33));
    xk1=xk1ss-(xk1ss-xk1)*exp(-dt/txk1);
    double rk1=1.0/(1.0+exp((v+105.8-2.6*ko)/9.493));
    double GK1=0.1908;
    if (celltype==1)
    {
        GK1*=1.2;
    }
    if (celltype==2)
    {
        GK1*=1.3;
    }
    IK1=GK1*sqrt(ko)*rk1*xk1*(v-EK);
    
    double kna1=15.0;
    double kna2=5.0;
    double kna3=88.12;
    double kasymm=12.5;
    double wna=6.0e4;
    double wca=6.0e4;
    double wnaca=5.0e3;
    double kcaon=1.5e6;
    double kcaoff=5.0e3;
    double qna=0.5224;
    double qca=0.1670;
    double hca=exp((qca*v*F)/(R*T));
    double hna=exp((qna*v*F)/(R*T));
    double h1=1+nai/kna3*(1+hna);
    double h2=(nai*hna)/(kna3*h1);
    double h3=1.0/h1;
    double h4=1.0+nai/kna1*(1+nai/kna2);
    double h5=nai*nai/(h4*kna1*kna2);
    double h6=1.0/h4;
    double h7=1.0+nao/kna3*(1.0+1.0/hna);
    double h8=nao/(kna3*hna*h7);
    double h9=1.0/h7;
    double h10=kasymm+1.0+nao/kna1*(1.0+nao/kna2);
    double h11=nao*nao/(h10*kna1*kna2);
    double h12=1.0/h10;
    double k1=h12*cao*kcaon;
    double k2=kcaoff;
    double k3p=h9*wca;
    double k3pp=h8*wnaca;
    double k3=k3p+k3pp;
    double k4p=h3*wca/hca;
    double k4pp=h2*wnaca;
    double k4=k4p+k4pp;
    double k5=kcaoff;
    double k6=h6*cai*kcaon;
    double k7=h5*h2*wna;
    double k8=h8*h11*wna;
    double x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
    double x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
    double x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
    double x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
    double E1=x1/(x1+x2+x3+x4);
    double E2=x2/(x1+x2+x3+x4);
    double E3=x3/(x1+x2+x3+x4);
    double E4=x4/(x1+x2+x3+x4);
    double KmCaAct=150.0e-6;
    double allo=1.0/(1.0+pow(KmCaAct/cai,2.0));
    double zna=1.0;
    double JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
    double JncxCa=E2*k2-E1*k1;
    double Gncx=0.0008;
    if (celltype==1)
    {
        Gncx*=1.1;
    }
    if (celltype==2)
    {
        Gncx*=1.4;
    }
    INaCa_i=0.8*Gncx*allo*(zna*JncxNa+zca*JncxCa);
    
    h1=1+nass/kna3*(1+hna);
    h2=(nass*hna)/(kna3*h1);
    h3=1.0/h1;
    h4=1.0+nass/kna1*(1+nass/kna2);
    h5=nass*nass/(h4*kna1*kna2);
    h6=1.0/h4;
    h7=1.0+nao/kna3*(1.0+1.0/hna);
    h8=nao/(kna3*hna*h7);
    h9=1.0/h7;
    h10=kasymm+1.0+nao/kna1*(1+nao/kna2);
    h11=nao*nao/(h10*kna1*kna2);
    h12=1.0/h10;
    k1=h12*cao*kcaon;
    k2=kcaoff;
    k3p=h9*wca;
    k3pp=h8*wnaca;
    k3=k3p+k3pp;
    k4p=h3*wca/hca;
    k4pp=h2*wnaca;
    k4=k4p+k4pp;
    k5=kcaoff;
    k6=h6*cass*kcaon;
    k7=h5*h2*wna;
    k8=h8*h11*wna;
    x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
    x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
    x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
    x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
    E1=x1/(x1+x2+x3+x4);
    E2=x2/(x1+x2+x3+x4);
    E3=x3/(x1+x2+x3+x4);
    E4=x4/(x1+x2+x3+x4);
    KmCaAct=150.0e-6;
    allo=1.0/(1.0+pow(KmCaAct/cass,2.0));
    JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
    JncxCa=E2*k2-E1*k1;
    INaCa_ss=0.2*Gncx*allo*(zna*JncxNa+zca*JncxCa);
    
    INaCa=INaCa_i+INaCa_ss;
    
    double k1p=949.5;
    double k1m=182.4;
    double k2p=687.2;
    double k2m=39.4;
    k3p=1899.0;
    double k3m=79300.0;
    k4p=639.0;
    double k4m=40.0;
    double Knai0=9.073;
    double Knao0=27.78;
    double delta=-0.1550;
    double Knai=Knai0*exp((delta*v*F)/(3.0*R*T));
    double Knao=Knao0*exp(((1.0-delta)*v*F)/(3.0*R*T));
    double Kki=0.5;
    double Kko=0.3582;
    double MgADP=0.05;
    double MgATP=9.8;
    double Kmgatp=1.698e-7;
    double H=1.0e-7;
    double eP=4.2;
    double Khp=1.698e-7;
    double Knap=224.0;
    double Kxkur=292.0;
    double P=eP/(1.0+H/Khp+nai/Knap+ki/Kxkur);
    double a1=(k1p*pow(nai/Knai,3.0))/(pow(1.0+nai/Knai,3.0)+pow(1.0+ki/Kki,2.0)-1.0);
    double b1=k1m*MgADP;
    double a2=k2p;
    double b2=(k2m*pow(nao/Knao,3.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
    double a3=(k3p*pow(ko/Kko,2.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
    double b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
    double a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
    double b4=(k4m*pow(ki/Kki,2.0))/(pow(1.0+nai/Knai,3.0)+pow(1.0+ki/Kki,2.0)-1.0);
    x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
    x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
    x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
    x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
    E1=x1/(x1+x2+x3+x4);
    E2=x2/(x1+x2+x3+x4);
    E3=x3/(x1+x2+x3+x4);
    E4=x4/(x1+x2+x3+x4);
    double zk=1.0;
    double JnakNa=3.0*(E1*a3-E2*b3);
    double JnakK=2.0*(E4*b1-E3*a1);
    double Pnak=30;
    if (celltype==1)
    {
        Pnak*=0.9;
    }
    if (celltype==2)
    {
        Pnak*=0.7;
    }
    INaK=Pnak*(zna*JnakNa+zk*JnakK);
    
    double xkb=1.0/(1.0+exp(-(v-14.48)/18.34));
    double GKb=0.003;
    if (celltype==1)
    {
        GKb*=0.6;
    }
    IKb=GKb*xkb*(v-EK);
    
    double PNab=3.75e-10;
    INab=PNab*vffrt*(nai*exp(vfrt)-nao)/(exp(vfrt)-1.0);
    
    double PCab=2.5e-8;
    ICab=PCab*4.0*vffrt*(cai*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
    
    double GpCa=0.0005;
    IpCa=GpCa*cai/(0.0005+cai);



