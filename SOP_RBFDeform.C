#include "SOP_RBFDeform.h"

#include <GU/GU_Detail.h>
#include <OP/OP_Operator.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Matrix3.h>
#include <UT/UT_Matrix4.h>
#include <SYS/SYS_Math.h>
#include <stddef.h>


#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interpolation.h"

using namespace RBFDeform;

void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator(
        "RBFDeform",
        "RBF Deform",
        SOP_RBFDeform::myConstructor,
        SOP_RBFDeform::myTemplateList,
        3,
        3,
        0));
}

static PRM_Name  modelChoices[] =
{
    PRM_Name("0", "QNN"),
    PRM_Name("1", "Multilayer"),
    PRM_Name(0)
};

static PRM_Name  termChoices[] =
{
    PRM_Name("0", "Linear"),
    PRM_Name("1", "Constant"),
    PRM_Name("2", "Zero"),
    PRM_Name(0)
};


static PRM_ChoiceList  modelMenu(PRM_CHOICELIST_SINGLE, modelChoices);
static PRM_ChoiceList  termMenu(PRM_CHOICELIST_SINGLE,  termChoices);


static PRM_Name names[] = {
    PRM_Name("model",   "Model"),
    PRM_Name("term",    "RBF Term"),
    PRM_Name("radius",  "Radius"),
    PRM_Name("layers",  "Layers"),
    PRM_Name("lambda",  "Lambda"),
};

PRM_Template
SOP_RBFDeform::myTemplateList[] = {
    PRM_Template(PRM_STRING,    1, &PRMgroupName, 0, &SOP_Node::pointGroupMenu, 0, 0, \
        SOP_Node::getGroupSelectButton(GA_GROUP_POINT)),
    PRM_Template(PRM_ORD,   1, &names[0], 0, &modelMenu, 0, 0),
    PRM_Template(PRM_ORD,   1, &names[1], 0, &termMenu, 0, 0),
    PRM_Template(PRM_FLT_J,	1, &names[2], PRMoneDefaults, 0, &PRMscaleRange),
    PRM_Template(PRM_INT_J,	1, &names[3], PRMfourDefaults),
    PRM_Template(PRM_FLT_J,	1, &names[4], PRMpointOneDefaults),
    PRM_Template(),
};


OP_Node *
SOP_RBFDeform::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_RBFDeform(net, name, op);
}

SOP_RBFDeform::SOP_RBFDeform(OP_Network *net, const char *name, OP_Operator *op)
    : SOP_Node(net, name, op), myGroup(NULL)
{
   
    mySopFlags.setManagesDataIDs(true);
}

SOP_RBFDeform::~SOP_RBFDeform() {}

OP_ERROR
SOP_RBFDeform::cookInputGroups(OP_Context &context, int alone)
{
    
    return cookInputPointGroups(
        context, // This is needed for cooking the group parameter, and cooking the input if alone.
        myGroup, // The group (or NULL) is written to myGroup if not alone.
        alone,   // This is true iff called outside of cookMySop to update handles.
                 // true means the group will be for the input geometry.
                 // false means the group will be for gdp (the working/output geometry).
        true,    // (default) true means to set the selection to the group if not alone and the highlight flag is on.
        0,       // (default) Parameter index of the group field
        -1,      // (default) Parameter index of the group type field (-1 since there isn't one)
        true,    // (default) true means that a pointer to an existing group is okay; false means group is always new.
        false,   // (default) false means new groups should be unordered; true means new groups should be ordered.
        true,    // (default) true means that all new groups should be detached, so not owned by the detail;
                 //           false means that new point and primitive groups on gdp will be owned by gdp.
        0        // (default) Index of the input whose geometry the group will be made for if alone.
    );
}

OP_ERROR
SOP_RBFDeform::cookMySop(OP_Context &context)
{
    
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT)
        return error();

    duplicatePointSource(0, context);

    // Get rest and deform geometry:
    const GU_Detail *rest_gdp   = inputGeo(1);
    const GU_Detail *deform_gdp = inputGeo(2);

    // Point count should match:
    if (rest_gdp->getNumPoints() != deform_gdp->getNumPoints())
        return error();
    
    alglib::real_2d_array rbf_data_model;
    int numpoints = rest_gdp->getNumPoints();
    rbf_data_model.setlength(numpoints, 6);

    // Construct model data:
    GA_Offset ptoff;
    {   
        GA_FOR_ALL_PTOFF(rest_gdp, ptoff)
        {
            const UT_Vector3 restP   = rest_gdp->getPos3(ptoff);
            const UT_Vector3 deformP = deform_gdp->getPos3(ptoff);
            const UT_Vector3 delta   = UT_Vector3(deformP - restP);
            double data[6] = {restP.x(), restP.y(), restP.z(), delta.x(), delta.y(), delta.z()};
            for (int i=0; i<6; ++i)
                rbf_data_model[static_cast<int>(ptoff)][i] = data[i]; 
        }
    }

    fpreal t = context.getTime();

    // Params:
    UT_String modelName, termName;
    MODEL(modelName);
    TERM(termName);
    float radius = RADIUS(t);
    int   layers = LAYERS(t);
    float lambda = LAMBDA(t);
    int modelIndex = SYSatoi(modelName.buffer());
    int termIndex  = SYSatoi(termName.buffer());

    if (error() >= UT_ERROR_ABORT)
        return error();

    // Build model:
    alglib::rbfmodel model;
    alglib::rbfreport report;
    alglib::rbfcreate(3, 3, model);
    alglib::rbfsetpoints(model, rbf_data_model);

    // Select RBF model:
    if (modelIndex == 0)
        alglib::rbfsetalgoqnn(model);
    if (modelIndex == 1)
        alglib::rbfsetalgomultilayer(model, radius, layers, lambda);

    // Select RBF term:
    if (termIndex == 0)
        alglib::rbfsetlinterm(model);
    if (termIndex == 1)
        alglib::rbfsetconstterm(model);
    if (termIndex == 2)
        alglib::rbfsetzeroterm(model);


    // Finally build model:
    alglib::rbfbuildmodel(model, report);

    // Debug:
    const char *info;
    char info_buffer[200];
    sprintf(info_buffer, "Termination type: %d, Iterations: %d", static_cast<int>(report.terminationtype),
    static_cast<int>(report.iterationscount));
    info = &info_buffer[0];
    addMessage(SOP_MESSAGE, info);

    // Early quit if model wasn't built properly (singular matrix etc):
    if (static_cast<int>(report.terminationtype) == 0)
        return error();

    // Execute storage:
    alglib::real_1d_array coord;
    coord.setlength(3);
    alglib::real_1d_array result;
    result.setlength(3);


    // Here we determine which groups we have to work on.  We only
    // handle point groups.
    if (cookInputGroups(context) >= UT_ERROR_ABORT)
        return error();


    // Execute model
    GA_FOR_ALL_GROUP_PTOFF(gdp, myGroup, ptoff)
    {
        UT_Vector3 pos = gdp->getPos3(ptoff);
        const double dp[3] = {pos.x(), pos.y(), pos.z()};
        coord.setcontent(3, dp);
        alglib::rbfcalc(model, coord, result);
        const UT_Vector3 delta = UT_Vector3(result[0], result[1],result[2]);
        gdp->setPos3(ptoff, pos+delta);
    }

    // If we've modified P, and we're managing our own data IDs,
    // we must bump the data ID for P.
    if (!myGroup || !myGroup->isEmpty())
        gdp->getP()->bumpDataId();

    return error();
}
