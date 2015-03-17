#ifndef _XBridge4SM_H_
#define _XBridge4SM_H_


namespace LifeV
{

class XBridge4SM
{

public:

    XBridge4SM();

    XBridge4SM ( const XBridge4SM& model );
    
    virtual ~XBridge4SM() {}

    XBridge4SM& operator= ( const XBridge4SM& model );



private:

    // Constants
    Real M_K2;
    Real M_K3;
    Real M_K4;
    Real M_Kd;
    Real M_Kdd;
    Real M_a1;
    Real b1;
    Real aa;
    Real ba;

}; // XBridge4SM

}

#endif
