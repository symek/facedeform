#pragma once 
#include <SOP/SOP_Node.h>

namespace facedeform {


class DBSE 
{
public:
    typedef std::vector<const GU_Detail*> BlendShapesV;
private:
    typedef Eigen::HouseholderQR<Eigen::MatrixXd> QRMatrix;
public: DBSE() {};
        ~DBSE();
        bool init(const SOP_Node * node, const BlendShapesV & shapes);
        bool isInitialized() { return myInitialized; }
        bool compute(const GA_Attribute * rest_attrib);
        bool isComputed() { return myComputed; }
        void displaceVector(const GA_Index & ptidx, UT_Vector3 & disp);
private:
    uint myNpoints = 0;
    bool myInitialized = false;
    bool myComputed = false;
    const SOP_Node * mySop;

    Eigen::MatrixXd myBlendShapesM;
    QRMatrix myQrMatrix;
    Eigen::VectorXd myDeltaV;
    Eigen::VectorXd myWeights;
};





} // end of facedeform namespace