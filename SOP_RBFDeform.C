#include <UT/UT_DSOVersion.h>
#include "SOP_RBFDeform.h"

#include <GU/GU_Detail.h>
#include <OP/OP_Operator.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <UT/UT_Matrix3.h>
#include <UT/UT_Matrix4.h>
#include <SYS/SYS_Math.h>
#include <stddef.h>


#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


using namespace RBFDeform;

void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator(
        "rbfDeform",
        "Radial Basis Deformer",
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
static PRM_Range       radiusRange(PRM_RANGE_PRM, 0, PRM_RANGE_PRM, 10);

static PRM_Name names[] = {
    PRM_Name("model",   "Model"),
    PRM_Name("term",    "RBF Term"),
    PRM_Name("qcoef",   "Q (Smoothness)"),
    PRM_Name("zcoef",   "Z (Deviation)"),
    PRM_Name("radius",  "Radius"),
    PRM_Name("layers",  "Layers"),
    PRM_Name("lambda",  "Lambda"),
    PRM_Name("tangent", "Tangent space"),
};

PRM_Template
SOP_RBFDeform::myTemplateList[] = {
    PRM_Template(PRM_STRING,    1, &PRMgroupName, 0, &SOP_Node::pointGroupMenu, 0, 0, \
        SOP_Node::getGroupSelectButton(GA_GROUP_POINT)),
    PRM_Template(PRM_ORD,   1, &names[0], 0, &modelMenu, 0, 0),
    PRM_Template(PRM_ORD,   1, &names[1], 0, &termMenu, 0, 0),
    PRM_Template(PRM_FLT_J, 1, &names[2], PRMoneDefaults),
    PRM_Template(PRM_FLT_J, 1, &names[3], PRMfiveDefaults),
    PRM_Template(PRM_FLT_J,	1, &names[4], PRMoneDefaults, 0, &radiusRange),
    PRM_Template(PRM_INT_J,	1, &names[5], PRMfourDefaults),
    PRM_Template(PRM_FLT_J,	1, &names[6], PRMpointOneDefaults),
    PRM_Template(PRM_TOGGLE,1, &names[7], PRMzeroDefaults),
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

    fpreal t = context.getTime();
    duplicatePointSource(0, context);

    // Get rest and deform geometry:
    const GU_Detail *rest_gdp   = inputGeo(1);
    const GU_Detail *deform_gdp = inputGeo(2);

    // Point count should match:
    if (rest_gdp->getNumPoints() != deform_gdp->getNumPoints())
    {
        addError(SOP_ERR_MISMATCH_POINT, "Rest and deform geometry should match.");
        return error();
    }

    #ifdef DEBUG
        Timer timer;
        timer.start();
    #endif

    alglib::real_2d_array rbf_data_model;
    int numpoints = rest_gdp->getNumPoints();
    rbf_data_model.setlength(numpoints, 6);

    #ifdef DEBUG
        std::cout << "Storage allocated: " << timer.current() << std::endl;
    #endif

    // Construct model data:
    GA_Offset ptoff;
    {   
        GA_FOR_ALL_PTOFF(rest_gdp, ptoff)
        {
            const UT_Vector3 restP   = rest_gdp->getPos3(ptoff);
            const UT_Vector3 deformP = deform_gdp->getPos3(ptoff);
            const UT_Vector3 delta   = UT_Vector3(deformP - restP);
            double data[6] = {restP.x(), restP.y(), restP.z(), delta.x(), delta.y(), delta.z()};
            // std::cout << "POINT " << static_cast<uint>(ptoff) << ": ";
            // TODO: Why we have more ptoffs then numpoints in Houdini 15?
            // Does it have something to do with changes to ptoffs?
            if (static_cast<uint>(ptoff) < numpoints)
            {
                for (int i=0; i<6; ++i)
                    rbf_data_model[static_cast<uint>(ptoff)][i] = data[i]; 
            }
        }
    }

    #if 1

    #ifdef DEBUG
    std::cout << "Data built: " << timer.current() << std::endl;
    #endif

    // Parms:
    UT_String modelName, termName;
    MODEL(modelName);
    TERM(termName);
    float qcoef  = SYSmax(0.1,  SYSabs(QCOEF(t)));
    float zcoef  = SYSmax(0.1,  SYSabs(ZCOEF(t)));
    float radius = SYSmax(0.01, SYSabs(RADIUS(t)));
    int   layers = SYSmax(1,    SYSabs(LAYERS(t)));
    float lambda = SYSmax(0.01, SYSabs(LAMBDA(t)));
    int   tangent= TANGENT(t);
    int modelIndex = atoi(modelName.buffer());
    int termIndex  = atoi(termName.buffer());

    if (error() >= UT_ERROR_ABORT)
        return error();

    // Do we have tangents?
    GA_ROHandleV3       tangentu_h(gdp, GA_ATTRIB_POINT, "tangentu");
    GA_ROHandleV3       tangentv_h(gdp, GA_ATTRIB_POINT, "tangentv");
    GA_RWHandleV3       normals_h(gdp, GA_ATTRIB_POINT,  "N");

    if (tangent == 1 && (!tangentu_h.isValid() || !tangentv_h.isValid()))
         addWarning(SOP_MESSAGE, "Can't deform in tangent space without tangent[u/v] attribs.");

    if (!normals_h.isValid()) 
        gdp->normal(normals_h);

    // Create model objects and ralated items:
    std::string str_model;
    alglib::rbfmodel model;
    alglib::rbfreport report;
    alglib::rbfcreate(3, 3, model);
    alglib::rbfsetpoints(model, rbf_data_model);

    // Select RBF model:
    switch(modelIndex)
    {
        case ALGLIB_MODEL_QNN:
            alglib::rbfsetalgoqnn(model, qcoef, zcoef);
            break;
        case ALGLIB_MODEL_ML:
            alglib::rbfsetalgomultilayer(model, radius, layers, lambda);
            break;
    }

    // Select RBF term:
    switch(termIndex)
    {
        case ALGLIB_TERM_LINEAR:
            alglib::rbfsetlinterm(model);
            break;
        case ALGLIB_TERM_CONST:
            alglib::rbfsetconstterm(model);
            break;
        case ALGLIB_TERM_ZERO:
            alglib::rbfsetzeroterm(model);
            break;
    }


    // Finally build model:
    alglib::rbfbuildmodel(model, report);

    #ifdef DEBUG
        std::cout << "Model built: " << timer.current() << std::endl;
    #endif

    // Debug:
    char info_buffer[200];
    sprintf(info_buffer, "Termination type: %d, Iterations: %d", \
    static_cast<int>(report.terminationtype), static_cast<int>(report.iterationscount));
    const char *info = &info_buffer[0];
    addMessage(SOP_MESSAGE, info);

    // Early quit if model wasn't built properly (singular matrix etc):
    if (static_cast<int>(report.terminationtype) != 1)
    {
        addError(SOP_ERR_NO_DEFORM_EFFECT, "Bad matrix.");
        return error();
    }

    // We won't use this model directly (unless sigle threaded path was chosen in compile time)
    // Instead we serialize it as send std::string to threads to be recreated there for 
    // further calculation.
    alglib::rbfserialize(model, str_model);

    #ifdef DEBUG
        std::cout << "Model serialized: " << timer.current() << std::endl;
    #endif

    // Here we determine which groups we have to work on.  We only
    // handle point groups.
    if (cookInputGroups(context) >= UT_ERROR_ABORT)
        return error();

    // Execute mode directly:
    #ifdef NO_RBF_THREADS
    #ifdef DEBUG
    std::cout << "Single thread mode." << std::endl;
    #endif
    // Execute storage:
    alglib::real_1d_array coord("[0,0,0]");
    alglib::real_1d_array result("[0,0,0]");

    // Execute model
    GA_FOR_ALL_GROUP_PTOFF(gdp, myGroup, ptoff)
    {
        UT_Vector3 pos = gdp->getPos3(ptoff);
        const double dp[3] = {pos.x(), pos.y(), pos.z()};
        coord.setcontent(3, dp);
        alglib::rbfcalc(model, coord, result);

        UT_Vector3 displace = UT_Vector3(result[0], result[1], result[2]);

        const float distance = displace.length();
        if (tangent && tangentu_h.isValid() && \
            tangentv_h.isValid() && normals_h.isValid())
        {

            UT_Vector3 tangentu = tangentu_h.get(ptoff); 
            UT_Vector3 tangentv = tangentv_h.get(ptoff);
            UT_Vector3 normal   = normals_h.get(ptoff);
            tangentu.normalize(); tangentv.normalize(); normal.normalize();
            UT_Matrix3  tangent_space(tangentu.x(), tangentu.y(), tangentu.z(),
                                      tangentv.x(), tangentv.y(), tangentv.z(),
                                      normal.x(), normal.y(), normal.z());

            tangent_space = tangent_space.transposedCopy() * tangent_space;

            UT_Vector3 a1 = tangentu * tangent_space;
            UT_Vector3 a2 = tangentv * tangent_space;
            a1 = a1.normalize();
            a2 = a2.normalize();

            float da1 = displace.dot(a1);
            float da2 = displace.dot(a2);
            displace  = (a1 * da1 + a2 * da2);
        }

        gdp->setPos3(ptoff, pos + displace);
    }

    #else

    // or try it in parallel:
    const GA_Range range(gdp->getPointRange());
    rbfDeformThreaded(range, str_model, gdp);

    #endif

    #ifdef DEBUG
    std::cout << "Model executed: " << timer.current() << std::endl;
    #endif

    // If we've modified P, and we're managing our own data IDs,
    // we must bump the data ID for P.
    #endif

    if (!myGroup || !myGroup->isEmpty())
        gdp->getP()->bumpDataId();

    return error();
}
