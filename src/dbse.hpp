#pragma once 
#include <SOP/SOP_Node.h>

namespace facedeform {


class DBSE 
{
    typedef std::vector<const GU_Detail*> BlendShapesV;
public: DBSE();
        bool init(const SOP_Node * node, const BlendShapesV & shapes);
        bool isInitialized() { return initialized; }
        bool build();
        bool isBuilt() { return built; }
private:
    const SOP_Node * mySop;
    uint myNpoints = 0;
    BlendShapesV myBlendShapesV;
    Eigen::MatrixXd myBlendShapesM;
    Eigen::VectorXd myDeltaV;
    Eigen::VectorXd myWeights;
    bool initialized = false;
    bool built = false;

};





} // end of facedeform namespace