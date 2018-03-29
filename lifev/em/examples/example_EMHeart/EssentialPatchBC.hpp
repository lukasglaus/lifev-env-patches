




class EssentialPatchBCCircular
{
public:
    EssentialPatchBCCircular(){}
    ~EssentialPatchBCCircular(){}
    
    CreatePatch();
    setupPatchBC();
    modifyPatchBC();


protected:

    const std::string m_name;
    const int m_prevFaceFlag;
    const int m_currentPatchFlag;
    
    // BCFunctionBase m_bcFunctionBase;
    // BCFunctionDirectional m_bcFunctionDirectional;
    
    // PatchBCFunctionBaseCreator m_patchBCFunctionBaseCreator;
    
    Vector3D m_center { 0. , 0. , 0. };
    Real m_radius { 0. };
}