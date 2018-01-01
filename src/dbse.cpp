#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include "dbse.hpp"


namespace facedeform
{

bool DBSE::init(const SOP_Node * node, const BlendShapesV & shapes) 
{
	mySop = node;
	const GU_Detail * gdp = mySop->curGdp(0);
    myNpoints = gdp->getNumPoints();
    GA_ROHandleV3 rest_h(gdp, GA_ATTRIB_POINT, "rest");
    	
    if (rest_h.isInvalid()) {
    	return false;
    }

    myBlendShapesV = shapes;
    myBlendShapesM.resize(myNpoints*3, shapes.size());
    myDeltaV.resize(myNpoints*3);

    unsigned col = 0;
    GA_Offset ptoff;
    BlendShapesV::const_iterator it;
    for(it=shapes.begin(); it != shapes.end(); it++, ++col) {
        const GU_Detail * shape = *it;
        GA_FOR_ALL_PTOFF(shape, ptoff) {
            const UT_Vector3 rest_pos  = rest_h.get(ptoff);
            const GA_Index   rest_itx  = gdp->pointIndex(ptoff);
            const GA_Offset  shape_off = shape->pointOffset(rest_itx);
            const UT_Vector3 shape_pos = shape->getPos3(shape_off);
            const UT_Vector3 shape_delta(shape_pos - rest_pos);
            myBlendShapesM(3*rest_itx + 0, col) = shape_delta.x();
            myBlendShapesM(3*rest_itx + 1, col) = shape_delta.y(); 
            myBlendShapesM(3*rest_itx + 2, col) = shape_delta.z();
        }
    }

    GA_FOR_ALL_PTOFF(gdp, ptoff) {
        const GA_Index ptidx  = gdp->pointIndex(ptoff);
        const UT_Vector3 pos  = gdp->getPos3(ptoff);
        const UT_Vector3 rest = rest_h.get(ptoff);
        myDeltaV(3*ptidx + 0) = pos.x() - rest.x();
        myDeltaV(3*ptidx + 1) = pos.y() - rest.y();
        myDeltaV(3*ptidx + 2) = pos.z() - rest.z();
    }

    // for (unsigned i=3; i < mySop->nConnectedInputs(); ++i) {
    //     const GU_Detail* shape = mySop->inputGeo(i);
    //     if (shape->getNumPoints() != myNpoints) {
    //         continue;
    //     }
    //     myBlendShapesV.push_back(shape);
    // }
    initialized = true;
    return true;
}

bool DBSE::build()
{
	Eigen::HouseholderQR<Eigen::MatrixXd> orthonormal_mat(myBlendShapesM);
    // scalar product of delta and Q's columns:
    Eigen::MatrixXd weights_mat = myDeltaV.asDiagonal() * orthonormal_mat.matrixQR();
    // Get weights out of this: 
    myWeights = weights_mat.colwise().sum();
    built = true;
    return true;
}



}