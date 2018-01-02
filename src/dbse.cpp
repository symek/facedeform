#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <memory>
#include "dbse.hpp"

namespace facedeform
{

bool DirectBSEdit::init(const GU_Detail * gdp, const ShapesVector & shapes)  {
    myGdp = gdp;
    const GA_Size npoints = gdp->getNumPoints();
    myShapesMatrix.conservativeResize(npoints*3, shapes.size());
    myDeltaV.conservativeResize(npoints*3);

    unsigned col = 0;
    GA_Offset ptoff;
    ShapesVector::const_iterator it;
    for(it=shapes.begin(); it != shapes.end(); it++, ++col) {
        const GU_Detail * shape = *it;
        GA_FOR_ALL_PTOFF(shape, ptoff) {
            const GA_Index   rest_itx  = myGdp->pointIndex(ptoff);
            const UT_Vector3 rest_pos  = myGdp->getPos3(ptoff);
            const GA_Offset  shape_off = shape->pointOffset(rest_itx);
            const UT_Vector3 shape_pos = shape->getPos3(shape_off);
            const UT_Vector3 shape_delta(shape_pos - rest_pos);
            myShapesMatrix(3*rest_itx + 0, col) = shape_delta.x();
            myShapesMatrix(3*rest_itx + 1, col) = shape_delta.y(); 
            myShapesMatrix(3*rest_itx + 2, col) = shape_delta.z();
        }
    }
    myQrMatrix = std::move(QRMatrixPtr(new QRMatrix(myShapesMatrix)));
    myInitialized = true;
    myComputed    = false; // we need to recompute weights.
    return myInitialized;
}

bool DirectBSEdit::computeWeights(const GA_Attribute * rest_attrib) {
    GA_ROHandleV3 rest_h(rest_attrib);
    if (rest_h.isInvalid()) {
        return false;
    }
    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(myGdp, ptoff) {
        const GA_Index ptidx  = myGdp->pointIndex(ptoff);
        const UT_Vector3 pos  = myGdp->getPos3(ptoff);
        const UT_Vector3 rest = rest_h.get(ptoff);
        myDeltaV(3*ptidx + 0) = pos.x() - rest.x();
        myDeltaV(3*ptidx + 1) = pos.y() - rest.y();
        myDeltaV(3*ptidx + 2) = pos.z() - rest.z();
    }
    // scalar product of delta and Q's columns:
    // Get weights out of this: 
    Eigen::MatrixXd weights_mat = myDeltaV.asDiagonal() * myQrMatrix->matrixQR();
    Eigen::VectorXd * weights = new Eigen::VectorXd(weights_mat.colwise().sum());
    myWeights = std::move(WeightsVector(weights));
    myComputed = true;
    return myComputed;
}

void DirectBSEdit::displaceVector(const GA_Index & ptidx, UT_Vector3 & disp) 
{
    // TODO: add clamping
    const Eigen::VectorXd & weights = *(myWeights.get());
    UT_ASSERT(ptidx < myShapesMatrix.rows());
    for(int col=0; col<myShapesMatrix.cols(); ++col) {
        const float xd = myShapesMatrix(3*ptidx + 0, col);
        const float yd = myShapesMatrix(3*ptidx + 1, col);
        const float zd = myShapesMatrix(3*ptidx + 2, col);
        const float w  = weights(col); //!!!???
        disp += UT_Vector3(xd, yd, zd) * w * 3;
    }
    
}

bool DirectBSEdit::getWeights(UT_FprealArray & weights_array) 
{
    if (!myComputed) {
        return false;
    }

    const Eigen::VectorXd & weights = *(myWeights.get());
    weights_array.setSize(weights.rows());
    for(int i=0; i< weights.size(); ++i)
        weights_array(i) = weights(i); 
    return true;
}



} // end of namespace facedeform