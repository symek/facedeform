#pragma once 
#include <SOP/SOP_Node.h>

namespace facedeform {


class DBSE 
{
public:
    typedef std::vector<const GU_Detail*> ShapesVector;
private:
    typedef Eigen::HouseholderQR<Eigen::MatrixXd> QRMatrix;
    typedef std::unique_ptr<QRMatrix> QRMatrixPtr;
    typedef std::unique_ptr<Eigen::VectorXd> WeightsVector;
public: 
        bool init(const GU_Detail * gdp, const ShapesVector & shapes);
        bool isInitialized() { return myInitialized && myShapesMatrix.rows(); }
        bool computeWeights(const GA_Attribute * rest_attrib);
        bool isComputed() { return myComputed && myQrMatrix && myWeights; }
        void displaceVector(const GA_Index & ptidx, UT_Vector3 & disp);
        bool getWeights(UT_FprealArray & weights_array);
private:
    //uint myNpoints = 0;
    bool myInitialized = false;
    bool myComputed    = false;
    const GU_Detail * myGdp;

    Eigen::MatrixXd myShapesMatrix;
    Eigen::VectorXd myDeltaV;
    QRMatrixPtr myQrMatrix = NULL;
    WeightsVector myWeights = NULL;
};





} // end of facedeform namespace