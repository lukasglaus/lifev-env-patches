//
//  Circulation.h
//  
//
//  Created by Thomas Kummer on 16.03.15.
//
//

#ifndef ____Circulation__
#define ____Circulation__

#include <stdio.h>


namespace LifeV
{
    class Circulation
    {
        public:
        
        Circulation()
        {
            std::cout << "Circulation.h" << std::endl;
        }
        
        static Real computePressure(Real Volume)
        {
            std::cout << "Circulation - Compute Pressure B.C. ... ";
            
            std::cout << "Done." << std::endl;
            
            return 1/Volume;
        }
        
        
        private:


    };

}

#endif /* defined(____Circulation__) */