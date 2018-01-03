#include "interpolation.h"
#include <UT/UT_DSOVersion.h>

#include <GU/GU_Detail.h>
#include <OP/OP_Operator.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <UT/UT_Matrix3.h>
#include <UT/UT_Matrix4.h>
#include <SYS/SYS_Math.h>
#include <stddef.h>
#include <GQ/GQ_Detail.h>
#include <GEO/GEO_PointTree.h>
#include <UT/UT_Color.h>
#include <GA/GA_AIFCopyData.h>

#include <unordered_map>
#include <memory>

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include "dbse.hpp"
#include "capture.hpp"
#include "SOP_FaceDeform.hpp"

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


using namespace facedeform;

void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator(
        "facedeform",
        "Face Deform",
        SOP_FaceDeform::myConstructor,
        SOP_FaceDeform::myTemplateList,
        3,
        1000,
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

static PRM_Default ZeroOneDefaults[] =
{
    PRM_Default(0),
    PRM_Default(1),
};

const char * model_help = "QNN and Multilayer are different algorithms to perform RBF \
 interpolation in ALGLIB. Multilayer is more robust and thus more expensive.";

const char * term_help = "By appending small linear or constant term to RBF system one \
 can stabalize it and help to solve smooth solution.";

const char * radius_help = "Radius controls not only RBF solution (how far to reach for a scattered data), \
 but also radius of deformation applied to geometry.";

const char * maxedges_help = "Number of edges deformation affects geometry. This is applied before radius.";

const char * tangent_help = "Project deformation into tangential space of a rest geometry. This helps to remove \
 extreme deformations. ";

const char * morphspace_help = "Projects deformation into subspace defined by blendshapes of a base mesh. \
 They should be connected after second and third input and match rest mesh topology (unlike 2d and 3rd input which are\
    typically sparser than first input (rest pose geo)).";

const char * weightrange_help = "Clamps total blendshape weights, so that deformation \
 will be constrained strictly to blends' poses.";

 const char * falloff_help = "This is exponent of distance ratio (distance / radius) \
  with which displacement falls off.";

static PRM_ChoiceList  modelMenu(PRM_CHOICELIST_SINGLE, modelChoices);
static PRM_ChoiceList  termMenu(PRM_CHOICELIST_SINGLE,  termChoices);
static PRM_Range       radiusRange(PRM_RANGE_RESTRICTED, 0.0, PRM_RANGE_UI, 10.0);
static PRM_Range       falloffRange(PRM_RANGE_RESTRICTED, 0.0, PRM_RANGE_UI, 2.0);


static PRM_Name names[] = {
    PRM_Name("model",   "Model"),
    PRM_Name("term",    "RBF Term"),
    PRM_Name("qcoef",   "Q (Smoothness)"),
    PRM_Name("zcoef",   "Z (Deviation)"),
    PRM_Name("radius",  "Radius"),
    PRM_Name("layers",  "Layers"),
    PRM_Name("lambda",  "Lambda"),
    PRM_Name("tangent", "Tangent space"),
    PRM_Name("morphspace","Blendshapes subspace"),
    PRM_Name("maxedges",      "Max edges"),
    PRM_Name("doclampweight", "Clamp weights"),
    PRM_Name("weightrange", "Range"),
    PRM_Name("dofalloff",   "Falloff"),
    PRM_Name("falloffradius", "Falloff radius"),
    PRM_Name("falloffrate", "Falloff rate (exponent)"),
};

PRM_Template
SOP_FaceDeform::myTemplateList[] = {
    PRM_Template(PRM_STRING,    1, &PRMgroupName, 0, &SOP_Node::pointGroupMenu, 0, 0, \
        SOP_Node::getGroupSelectButton(GA_GROUP_POINT)),
    PRM_Template(PRM_ORD,   1, &names[0], 0, &modelMenu, 0, 0, 0, 0, model_help),
    PRM_Template(PRM_ORD,   1, &names[1], 0, &termMenu, 0, 0, 0, 0, term_help),
    PRM_Template(PRM_FLT_J, 1, &names[2], PRMoneDefaults),
    PRM_Template(PRM_FLT_J, 1, &names[3], PRMfiveDefaults),
    PRM_Template(PRM_FLT_LOG,1,&names[4], PRMoneDefaults, 0, &radiusRange, 0, 0, 0, radius_help), // radius
    PRM_Template(PRM_INT_J, 1, &names[9], PRMfourDefaults, 0, 0, 0, 0, 0, maxedges_help), // maxedges
    PRM_Template(PRM_INT_J, 1, &names[5], PRMfourDefaults), // layers
    PRM_Template(PRM_FLT_J, 1, &names[6], PRMpointOneDefaults), // lambda
    PRM_Template(PRM_TOGGLE,1, &names[7], PRMzeroDefaults, 0, 0, 0, 0, 0, tangent_help), // tangent
    PRM_Template(PRM_TOGGLE,1, &names[8], PRMzeroDefaults, 0, 0, 0, 0, 0, morphspace_help), // morphspace
    PRM_Template(PRM_TOGGLE,1, &names[10], PRMzeroDefaults, 0, 0, 0, 0, 0, weightrange_help),
    PRM_Template(PRM_FLT_J, 2, &names[11], ZeroOneDefaults, 0, 0, 0, 0, 0, weightrange_help),
    PRM_Template(PRM_TOGGLE,1, &names[12], PRMzeroDefaults, 0, 0, 0, 0, 0),
    PRM_Template(PRM_FLT_LOG,1,&names[13], PRMoneDefaults, 0, &radiusRange, 0, 0, 0, radius_help), // falloffradius
    PRM_Template(PRM_FLT_J, 1, &names[14], PRMoneDefaults, 0, &falloffRange, 0, 0, 0, falloff_help), // falloff rate
    PRM_Template(),
};


OP_Node *
SOP_FaceDeform::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_FaceDeform(net, name, op);
}

SOP_FaceDeform::SOP_FaceDeform(OP_Network *net, const char *name, OP_Operator *op)
    : SOP_Node(net, name, op), myGroup(NULL)
{
   
    mySopFlags.setManagesDataIDs(true);
}

SOP_FaceDeform::~SOP_FaceDeform() {}

OP_ERROR
SOP_FaceDeform::cookInputGroups(OP_Context &context, int alone)
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

void SOP_FaceDeform::setupBlends(OP_Context &context)
{
     // Copy rest only if rest pose changed.
    int rest_pose_changed = 0;
    GA_Attribute * rest = gdp->findFloatTuple(GA_ATTRIB_POINT, "rest", 3);
    checkChangedSourceFlags(0, context, &rest_pose_changed);
    if (rest_pose_changed || !rest) {   
        rest = gdp->addRestAttribute(GA_ATTRIB_POINT);
        const GA_Attribute * pos = gdp->getP();
        rest->replace(*pos);
    }
    // Check if there are any changes to blendshapes
    int blends_changed = 0;
    for (unsigned input=3; input < nConnectedInputs(); ++input) {
        int changed = 0;
        checkChangedSourceFlags(input, context, &changed);
        if (changed) {
            blends_changed = 1;
            break;
        }
    }
    // Gather blendshapes
    if (blends_changed || !m_direct_blends.isInitialized()) {
        DirectBSEdit::ShapesVector shapes;
        for (unsigned i=3; i < nConnectedInputs(); ++i) {
            const GU_Detail* shape = inputGeo(i);
            if (shape->getNumPoints() != gdp->getNumPoints()) {
                addWarning(SOP_ERR_MISMATCH_POINT, \
                    "Some blendshapes don't match rest pose point count. Ignoring them.");
                continue;
            }
            shapes.push_back(shape);
        }
        // Initialize DirectBSEdit matrix 
        if (!m_direct_blends.init(gdp, shapes)) {
             addWarning(SOP_MESSAGE, "Can't proceed with morph space deformation. Ingoring it.");
        }
    }
}

OP_ERROR
SOP_FaceDeform::cookMySop(OP_Context &context)
{
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT)
        return error();
    // In case we'd like to trigger recapture mesh once rest pose changed.
    int rest_pose_changed = 0;
    checkChangedSourceFlags(0, context, &rest_pose_changed);
    DEBUG_PRINT("cookMySop mesh changed? %i\n", rest_pose_changed );
    fpreal t = context.getTime();
    duplicatePointSource(0, context);
    // Get rest and deform geometry:
    const GU_Detail *rest_control_rig   = inputGeo(1);
    const GU_Detail *deform_control_rig = inputGeo(2);
    // Points count in control rig should match:
    if (rest_control_rig->getNumPoints() != deform_control_rig->getNumPoints()) {
        addError(SOP_ERR_MISMATCH_POINT, "Rest and deform geometry should match.");
        return error();
    }
    // Setup tracker to keep track of changes in rig.
    int rest_rig_changed = 0;
    if (m_input_tracker.size() != 2) {
        m_input_tracker.push_back(InputGeoID(rest_control_rig));
        m_input_tracker.push_back(InputGeoID(deform_control_rig));
        rest_rig_changed = 1;
    }
    /// UI
    // ALGLIB parms:
    UT_String modelName, termName;
    MODEL(modelName); 
    TERM(termName);
    const int model_index  = atoi(modelName.buffer());
    const int term_index   = atoi(termName.buffer());
    const float qcoef  = SYSmax(0.1,  QCOEF(t));
    const float zcoef  = SYSmax(0.1,  ZCOEF(t));
    const float radius = SYSmax(0.01, RADIUS(t));
    const int   layers = SYSmax(1,    LAYERS(t));
    const float lambda = SYSmax(0.01, LAMBDA(t));
    // Application parms:
    const int tangent_disp = TANGENT(t);
    const int morph_space  = MORPHSPACE(t);
    const int max_edges    = SYSmax(1, MAXEDGES(t));
    UT_Vector2 weightrange;
    const int doclampweight = DOCLAMPWEIGHT(t);
    WEIGHTRANGE(t, weightrange);
    const int   dofalloff     = DOFALLOFF(t);
    const float falloffrate   = FALLOFFRATE(t);
    const float falloffradius = FALLOFFRADIUS(t);

    if (error() >= UT_ERROR_ABORT)
        return error();

    const int rest_npoints = rest_control_rig->getNumPoints();
    alglib::real_2d_array rbf_data_model;
    rbf_data_model.setlength(rest_npoints, 6);
    // Construct model data:
    GA_Offset ptoff;
    {   
        GA_FOR_ALL_PTOFF(rest_control_rig, ptoff)
        {
            const UT_Vector3 restP   = rest_control_rig->getPos3(ptoff);
            const UT_Vector3 deformP = deform_control_rig->getPos3(ptoff);
            const UT_Vector3 delta   = UT_Vector3(deformP - restP);
            double data[6] = {restP.x(), restP.y(), restP.z(), delta.x(), 
                delta.y(), delta.z()};
            const GA_Index ptidx = rest_control_rig->pointIndex(ptoff); 
            if (ptidx < rest_npoints) {
                for (int i=0; i<6; ++i)
                    rbf_data_model[ptidx][i] = data[i]; 
            }
        }
    }
    // Do we have tangents?
    GA_ROHandleV3  tangentu_h(gdp, GA_ATTRIB_POINT, "tangentu");
    GA_ROHandleV3  tangentv_h(gdp, GA_ATTRIB_POINT, "tangentv");
    GA_ROHandleV3  normals_h(gdp, GA_ATTRIB_POINT,  "N");
    // Do we do tangent projection?
    const bool do_tangent_disp = tangent_disp  && tangentu_h.isValid() && \
        tangentv_h.isValid() && normals_h.isValid();
    if (tangent_disp  && !do_tangent_disp) {
         addWarning(SOP_MESSAGE, "Append PolyFrameSOP and enable tangent[u/v] \
            and N attribute to allow tangent displacement.");
    }
    // Proximity capture.
    // We keep track of dataid in m_input_tracker, at 0 is our control rig...
    if (!(m_input_tracker.at(0) == rest_control_rig)) {
        m_input_tracker.at(0).update();
        m_input_tracker.at(1).update();
        rest_rig_changed = 1;
    }
    DEBUG_PRINT("rest rig changed? %i\n", rest_rig_changed );
    DEBUG_PRINT("m_mesh_capture init: %i, captured: %i\n", \
        m_mesh_capture.isInitialized() , m_mesh_capture.isCaptured());

    if (rest_pose_changed || rest_rig_changed \
        || !m_mesh_capture.isInitialized() || !m_mesh_capture.isCaptured()) {

        if(!m_mesh_capture.init(gdp, rest_control_rig)) {
            addError(SOP_MESSAGE, "Can't initialize geometry to capture with a rig!");
            return error();
        }
        if(!m_mesh_capture.capture(max_edges, radius, dofalloff, falloffrate)) {
            addError(SOP_MESSAGE, "Can't capture geometry with a rig!");
            return error();
        }
    }
    // Morph space deformation.
    // we need to store rest P.
    if (morph_space && (nConnectedInputs() > 3)) {
        setupBlends(context);
    } else if (morph_space) {
        addWarning(SOP_MESSAGE, "No blendshapes found. Ignoring morphspace deformation.");
    }
    // Create model objects and ralated items:
    std::string str_model;
    alglib::rbfmodel model;
    alglib::rbfreport report;
    try {
        alglib::rbfcreate(3, 3, model);
        alglib::rbfsetpoints(model, rbf_data_model);
    } catch (alglib::ap_error err) {
        addError(SOP_ERR_NO_DEFORM_EFFECT, "Can't build RBF model.");
        return error();
    }
    // Select RBF model:
    switch(model_index) {
        case ALGLIB_MODEL_QNN:
            alglib::rbfsetalgoqnn(model, qcoef, zcoef);
            break;
        case ALGLIB_MODEL_ML:
            alglib::rbfsetalgomultilayer(model, radius, layers, lambda);
            break;
    }
    // Select RBF term:
    switch(term_index) {
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
    // Early quit if model wasn't built properly (singular matrix etc):
    if (static_cast<int>(report.terminationtype) != 1) {
        addError(SOP_ERR_NO_DEFORM_EFFECT, "Can't solve the problem.");
        return error();
    }
    // Debug:
    char info_buffer[200];
    sprintf(info_buffer, "Termination type: %d, Iterations: %d", \
        static_cast<int>(report.terminationtype),static_cast<int>(report.iterationscount));
    addMessage(SOP_MESSAGE, &info_buffer[0]);
    // We won't use this model directly (unless sigle threaded path was chosen in compile time)
    // Instead we serialize it as send std::string to threads to be recreated there for 
    // further calculation.
    // alglib::rbfserialize(model, str_model);
    // Here we determine which groups we have to work on.  We only
    // handle point groups.
    if (cookInputGroups(context) >= UT_ERROR_ABORT)
        return error();
    #if 1 
    // Execute storage:
    alglib::real_1d_array coord("[0,0,0]");
    alglib::real_1d_array result("[0,0,0]");
    GA_Attribute * cd_a = gdp->addFloatTuple(GA_ATTRIB_POINT, "Cd", 3, \
        GA_Defaults(GA_STORE_REAL32, 3, 1.f, 1.f, 1.f));
    GA_RWHandleV3 cd_h(cd_a);
    UT_Color affected_clr(UT_RGB, 1.f,1.f,1.f);
    // FIXME This doesn't work.
    // const GA_Attribute * cd_capture = m_mesh_capture.getColorAttribute();
    // if (getPicked() && cd_a && cd_capture) {
    //     const GA_AIFCopyData * copy_attrib = cd_a->getAIFCopyData();
    //     copy_attrib->copy(*cd_a, gdp->getPointRange(), *cd_capture, gdp->getPointRange());
    // }
    GA_Attribute * distance_a = m_mesh_capture.getDistanceAttribute();
    if (!distance_a) {
        addWarning(SOP_MESSAGE, "Can't find distance capture attribute. Won't apply radius nor falloff.");
    }
    GA_ROHandleF distance_h(distance_a);
    GA_RWHandleF fd_falloff_h(gdp->addFloatTuple(GA_ATTRIB_POINT, "fd_falloff", 1));
    const float radius_sqrt = radius*radius;
    GA_FOR_ALL_PTOFF(gdp, ptoff) {
        float distance_sqrt = 0.f;
        if (distance_h.isValid())
            distance_sqrt = distance_h.get(ptoff);
        if (distance_sqrt > radius_sqrt) {
            continue;
        }
        const UT_Vector3 pos = gdp->getPos3(ptoff);
        const double dp[3] = {pos.x(), pos.y(), pos.z()};
        coord.setcontent(3, dp);
        alglib::rbfcalc(model, coord, result);
        UT_Vector3 displace(result[0], result[1], result[2]);
        if (do_tangent_disp) {
            UT_Vector3 u = tangentu_h.get(ptoff); 
            UT_Vector3 v = tangentv_h.get(ptoff);
            UT_Vector3 n = normals_h.get(ptoff);
            u.normalize(); v.normalize(); n.normalize();
            project_to_tangents(u, v, n, displace);
        }
        float falloff = SYSmin(distance_sqrt/radius_sqrt, 1.f);
        falloff = SYSpow(1.f - falloff, falloffrate);
        fd_falloff_h.set(ptoff, falloff);
        // Set pseudo color only when SOP_Node is selected. FIXME
        if (getPicked()) {
            float hue = SYSfit(falloff, 0.f, 1.f, 200.f, 360.f);
            affected_clr.setHSV(hue, 1.f, 1.f);
        } else {
            affected_clr.setRGB(1.f,1.f,1.f);
        }

        if (cd_h.isValid()) {
            cd_h.set(ptoff, affected_clr.rgb());
        }
        displace *= falloff;
        gdp->setPos3(ptoff, pos + displace);
    }
    #else
    ////
    #endif
    // Any inputs above 2 is considered as morph targets...
    if (morph_space && (nConnectedInputs() > 3) && m_direct_blends.isInitialized() ) {
        const GA_Attribute * rest = gdp->findFloatTuple(GA_ATTRIB_POINT, "rest", 3);
        GA_ROHandleV3 rest_h(rest);
        bool weights_done = false;
        if(!m_direct_blends.isComputed()) {
            weights_done = m_direct_blends.computeWeights(rest);
        }
        if (!weights_done) {
            addWarning(SOP_MESSAGE, "Can't compute weights for morphspace deformation. Ingoring it.");
        } else {
            GA_FOR_ALL_PTOFF(gdp, ptoff) {
                const GA_Index ptidx = gdp->pointIndex(ptoff);
                UT_Vector3 displace(0,0,0);
                m_direct_blends.displaceVector(ptidx, displace);
                const UT_Vector3 rest = rest_h.get(ptoff);
                const UT_Vector3 pos  = gdp->getPos3(ptoff);
                if (dofalloff) {
                    UT_Vector3 delta(pos - rest);
                    displace += delta*falloffradius;
                }
                gdp->setPos3(ptoff, pos + displace);
            }
            // copy blendshape's weights into detail attribute
            GA_Attribute * w_attrib = gdp->addFloatArray(GA_ATTRIB_DETAIL, "weights", 1);
            const GA_AIFNumericArray * w_aif = w_attrib->getAIFNumericArray();
            UT_FprealArray weights_array;
            if(m_direct_blends.getWeights(weights_array)) {
                w_aif->set(w_attrib, 0, weights_array);
                w_attrib->bumpDataId();
            }
        }
    }
    // If we've modified P, and we're managing our own data IDs,
    // we must bump the data ID for P.
    gdp->destroyAttribute(GA_ATTRIB_POINT, "distance");
    if (!myGroup || !myGroup->isEmpty())
        gdp->getP()->bumpDataId();

    return error();
}
