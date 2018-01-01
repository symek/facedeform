#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include "dbse.hpp"

namespace facedeform
{

bool DBSE::init(const SOP_Node * node, const BlendShapesV & shapes)  {
    mySop = node;
    const GU_Detail * gdp = mySop->curGdp(0);
    myNpoints = gdp->getNumPoints();
            
    myBlendShapesM.resize(myNpoints*3, shapes.size());
    myDeltaV.resize(myNpoints*3);

    unsigned col = 0;
    GA_Offset ptoff;
    BlendShapesV::const_iterator it;
    for(it=shapes.begin(); it != shapes.end(); it++, ++col) {
        const GU_Detail * shape = *it;
        GA_FOR_ALL_PTOFF(shape, ptoff) {
            const GA_Index   rest_itx  = gdp->pointIndex(ptoff);
            const UT_Vector3 rest_pos  = gdp->getPos3(ptoff);
            const GA_Offset  shape_off = shape->pointOffset(rest_itx);
            const UT_Vector3 shape_pos = shape->getPos3(shape_off);
            const UT_Vector3 shape_delta(shape_pos - rest_pos);
            myBlendShapesM(3*rest_itx + 0, col) = shape_delta.x();
            myBlendShapesM(3*rest_itx + 1, col) = shape_delta.y(); 
            myBlendShapesM(3*rest_itx + 2, col) = shape_delta.z();
        }
    }

    myQrMatrix.compute(myBlendShapesM);
    myInitialized = true;
    return true;
}

bool DBSE::compute(const GA_Attribute * rest_attrib) {
    const GU_Detail * gdp = mySop->curGdp(0);
    GA_ROHandleV3 rest_h(rest_attrib);
    if (rest_h.isInvalid()) {
        return false;
    }

    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(gdp, ptoff) {
        const GA_Index ptidx  = gdp->pointIndex(ptoff);
        const UT_Vector3 pos  = gdp->getPos3(ptoff);
        const UT_Vector3 rest = rest_h.get(ptoff);
        myDeltaV(3*ptidx + 0) = pos.x() - rest.x();
        myDeltaV(3*ptidx + 1) = pos.y() - rest.y();
        myDeltaV(3*ptidx + 2) = pos.z() - rest.z();
    }
    // scalar product of delta and Q's columns:
    Eigen::MatrixXd weights_mat = myDeltaV.asDiagonal() * myQrMatrix.matrixQR();
    // Get weights out of this: 
    myWeights = weights_mat.colwise().sum();
    myComputed = true;
    return true;
}

void DBSE::displaceVector(const GA_Index & ptidx, UT_Vector3 & disp) 
{
    // TODO: add clamping
    for(int col=0; col<myBlendShapesM.size(); ++col) {
        const float xd = myBlendShapesM(3*ptidx + 0, col);
        const float yd = myBlendShapesM(3*ptidx + 1, col);
        const float zd = myBlendShapesM(3*ptidx + 2, col);
        const float w  = myWeights(col)*3; //!!!???
        disp += UT_Vector3(xd, yd, zd) * w;
    }
    
}



} // end of namespace facedeform