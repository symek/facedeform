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

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include "SOP_FaceDeform.hpp"

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unordered_map>
#include <memory>


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
    PRM_Template(PRM_FLT_LOG,1, &names[4], PRMoneDefaults, 0, &radiusRange, 0, 0, 0, radius_help), // radius
    PRM_Template(PRM_INT_J, 1, &names[9], PRMfourDefaults, 0, 0, 0, 0, 0, maxedges_help), // maxedges
    PRM_Template(PRM_INT_J, 1, &names[5], PRMfourDefaults), // layers
    PRM_Template(PRM_FLT_J, 1, &names[6], PRMpointOneDefaults), // lambda
    PRM_Template(PRM_TOGGLE,1, &names[7], PRMzeroDefaults, 0, 0, 0, 0, 0, tangent_help), // tangent
    PRM_Template(PRM_TOGGLE,1, &names[8], PRMzeroDefaults, 0, 0, 0, 0, 0, morphspace_help), // morphspace
    PRM_Template(PRM_TOGGLE,1, &names[10], PRMzeroDefaults, 0, 0, 0, 0, 0, weightrange_help),
    PRM_Template(PRM_FLT_J, 2, &names[11], ZeroOneDefaults, 0, 0, 0, 0, 0, weightrange_help),
    PRM_Template(PRM_TOGGLE,1, &names[12], PRMzeroDefaults, 0, 0, 0, 0, 0),
    PRM_Template(PRM_FLT_J, 1, &names[13], PRMoneDefaults, 0, &falloffRange, 0, 0, 0, falloff_help), // falloff rate
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

OP_ERROR
SOP_FaceDeform::cookMySop(OP_Context &context)
{
    
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT)
        return error();

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

    const int   dofalloff   = DOFALLOFF(t);
    const float falloffrate = FALLOFFRATE(t);

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

    // In case we do blendspace projection we need to store rest P.
    if (morph_space && (nConnectedInputs() > 3)) {
         GA_Attribute * rest = gdp->addRestAttribute(GA_ATTRIB_POINT);
         const GA_Attribute * pos = gdp->getP();
         rest->replace(*pos);
    } else if (morph_space) {
        addWarning(SOP_MESSAGE, "No Blendshapes found. Ignoring morphspace deformation.");
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

    // Add temp color to visualize affected regions
    // FIXME: this doesn't work correctly.
    UT_Color affected_clr;
    GA_RWHandleV3 cd_h;
    if (getPicked()) {
        cd_h = GA_RWHandleV3(gdp->findDiffuseAttribute(GA_ATTRIB_POINT));
        if (!cd_h.isValid()) {
            cd_h = GA_RWHandleV3(gdp->addDiffuseAttribute(GA_ATTRIB_POINT));
        }
    }

    // Helper distance attribute per target point ot closest rig point.
    GA_RWHandleF dist_h(gdp->addFloatTuple(GA_ATTRIB_POINT, "distance", 3));

    // Poor's man geodesic distance
    GQ_Detail     gq_detail(gdp);
    GEO_PointTree gdp_tree, rest_tree;
    gdp_tree.build(gdp);
    rest_tree.build(rest_control_rig);
  
    { 
        typedef std::unique_ptr<GA_PointGroup> GA_PointGroupPtr;
        typedef std::unordered_map<int, GA_PointGroupPtr> HandlerGroupMap;
       
        HandlerGroupMap handlers_map;
        GA_PointGroup partial_group(*gdp);

        GA_ROHandleI class_h(rest_control_rig->findIntTuple(GA_ATTRIB_POINT, "class", 1));
        if(class_h.isInvalid()) {
            GA_PointGroupPtr handle_group(new GA_PointGroup(*gdp));
            handlers_map.insert(std::make_pair<int, \
                GA_PointGroupPtr>(0, std::move(handle_group)));
        }
 
        GA_FOR_ALL_PTOFF(rest_control_rig, ptoff) { 
            const UT_Vector3 anchor_pos  = rest_control_rig->getPos3(ptoff);
            const GA_Index  target_idx   = gdp_tree.findNearestIdx(anchor_pos);
            const GA_Offset target_ptoff = gdp->pointOffset(target_idx);

            int handle_id = 0;
            if (class_h.isValid()) {
                handle_id = class_h.get(ptoff); 
            }
            if(handlers_map.find(handle_id) == handlers_map.end()) {
                GA_PointGroupPtr handle_group(new GA_PointGroup(*gdp));
                handlers_map.insert(std::pair<int, \
                    GA_PointGroupPtr>(handle_id, std::move(handle_group)));  
            }
            gq_detail.groupEdgePoints(target_ptoff, max_edges, partial_group);
            GA_PointGroup * handle_group = handlers_map[handle_id].get();
            handle_group->combine(&partial_group);
            partial_group.clear();
        }
        

        // Execute storage:
        alglib::real_1d_array coord("[0,0,0]");
        alglib::real_1d_array result("[0,0,0]");
        HandlerGroupMap::const_iterator it;
        for (it = handlers_map.begin(); it != handlers_map.end(); it++) {
            const GA_PointGroupPtr & affected_group = it->second;
            for (GA_Size i=0; i<affected_group->entries(); ++i)  {
                const GA_Offset  target_ptoff = affected_group->findOffsetAtGroupIndex(i);
                const UT_Vector3 target_pos   = gdp->getPos3(target_ptoff);
                const GA_Index   anchor_idx   = rest_tree.findNearestIdx(target_pos);
                const GA_Offset  anchor_ptoff = rest_control_rig->pointOffset(anchor_idx);
                const UT_Vector3 anchor_pos   = rest_control_rig->getPos3(anchor_ptoff); 
                const float      distance     = distance3d(target_pos, anchor_pos);
                dist_h.set(target_ptoff, distance);

                if (distance > radius)
                    continue;

                const double dp[3] = {target_pos.x(), target_pos.y(), target_pos.z()};

                coord.setcontent(3, dp);
                alglib::rbfcalc(model, coord, result);
                UT_Vector3 displace(result[0], result[1], result[2]);

                float falloff = 1.f;
                if (dofalloff) {
                    falloff = SYSpow(1.f - (distance/radius), falloffrate);
                }
                displace *= falloff;

                if (getPicked() && cd_h.isValid()) {
                    float hue = SYSfit(1.f - falloff, 0.f, 1.f, 360.f, 200.f);
                    uint seed = static_cast<uint>(it->first);
                    hue += SYSfastRandom(seed)*10;
                    affected_clr.setHSV(hue, 1.f, 1.f);
                    cd_h.set(target_ptoff, affected_clr.rgb());
                }

                if (do_tangent_disp) {
                    UT_Vector3 u = tangentu_h.get(ptoff); 
                    UT_Vector3 v = tangentv_h.get(ptoff);
                    UT_Vector3 n = normals_h.get(ptoff);
                    u.normalize(); v.normalize(); n.normalize();
                    project_to_tangents(u, v, n, displace);
                }
                gdp->setPos3(target_ptoff, target_pos + displace);
            }
        }
    }


    // 
    // Any input above 2 is considered as morph targets...
    if (morph_space && (nConnectedInputs() > 3) )
    {
        const GA_Size npoints = gdp->getNumPoints();
        GA_ROHandleV3 rest_h(gdp, GA_ATTRIB_POINT, "rest");
        std::vector<const GU_Detail*> shapes;
        for (unsigned i=3; i < nConnectedInputs(); ++i) {
            const GU_Detail* shape = inputGeo(i);
            if (shape->getNumPoints() != npoints) {
                addWarning(SOP_ERR_MISMATCH_POINT, \
                    "Some blendshapes doesn't match rest pose point count!");
                continue;
            }
            shapes.push_back(shape);
        }

        Eigen::MatrixXd blends_mat(npoints*3, shapes.size());
        Eigen::VectorXd delta(npoints*3);
        create_blendshape_matrix(gdp, shapes, blends_mat, delta);
        // Orthonormalize blends
        Eigen::HouseholderQR<Eigen::MatrixXd> orthonormal_mat(blends_mat);
        // scalar product of delta and Q's columns:
        Eigen::MatrixXd weights_mat = delta.asDiagonal() * orthonormal_mat.matrixQR();
        // Get weights out of this: 
        Eigen::VectorXd weights = weights_mat.colwise().sum();
        // copy blendshape's weights into detail attribute
        GA_Attribute * w_attrib = gdp->addFloatArray(GA_ATTRIB_DETAIL, "weights", 1);
        const GA_AIFNumericArray * w_aif = w_attrib->getAIFNumericArray();
        UT_FprealArray weights_array(shapes.size());
        for(int i=0;i<shapes.size(); ++i) {
            weights_array.append(weights(i));
        }
        w_aif->set(w_attrib, 0, weights_array);
        w_attrib->bumpDataId();

        {
            GA_FOR_ALL_PTOFF(gdp, ptoff) {
                const GA_Size ptidx = gdp->pointIndex(ptoff);
                UT_Vector3 disp(0,0,0);
                for(int col=0; col<shapes.size(); ++col) {
                    const float xd = blends_mat(3*ptidx + 0, col);
                    const float yd = blends_mat(3*ptidx + 1, col);
                    const float zd = blends_mat(3*ptidx + 2, col);
                    const float w  = weights(col)*3; //!!!???
                    const float cw = doclampweight == 1 ? \
                        SYSclamp(w, weightrange.x(), weightrange.y()) : w; 
                    disp += UT_Vector3(xd, yd, zd) * cw;

                }
                UT_Vector3 pos = rest_h.get(ptoff);// gdp->getPos3(ptoff);
                gdp->setPos3(ptoff, pos + disp);
            }
        }

    }

    // Execute model
    // GA_FOR_ALL_GROUP_PTOFF(gdp, myGroup, ptoff)
    // {
    //     const UT_Vector3 pos = gdp->getPos3(ptoff);
    //     const double dp[3]   = {pos.x(), pos.y(), pos.z()};
    //     coord.setcontent(3, dp);
    //     alglib::rbfcalc(model, coord, result);
    //     UT_Vector3 displace = UT_Vector3(result[0], result[1], result[2]);

    //     if (do_tangent_disp)
    //     {
    //         UT_Vector3 u = tangentu_h.get(ptoff); 
    //         UT_Vector3 v = tangentv_h.get(ptoff);
    //         UT_Vector3 n = normals_h.get(ptoff);
    //         u.normalize(); v.normalize(); n.normalize();
    //         UT_Matrix3  b(u.x(), u.y(), u.z(),
    //                       v.x(), v.y(), v.z(),
    //                       n.x(), n.y(), n.z());

    //         b = b.transposedCopy() * b;
    //         UT_Vector3 a1(u * b); a1.normalize();
    //         UT_Vector3 a2(v * b); a2.normalize();
    //         const float da1 = displace.dot(a1);
    //         const float da2 = displace.dot(a2);
    //         displace        = UT_Vector3(a1 * da1 + a2 * da2);
    //     }

    //     gdp->setPos3(ptoff, pos + displace);
    // }


    // #else

    // or try it in parallel (no groups support yet)
    // const GA_Range range(gdp->getPointRange());
    // rbfDeformThreaded(range, str_model, do_tangent_disp, gdp);

    // #endif // end of NO_RBF_THREADS

    // If we've modified P, and we're managing our own data IDs,
    // we must bump the data ID for P.

    gdp->destroyAttribute(GA_ATTRIB_POINT, "distance");

    if (!myGroup || !myGroup->isEmpty())
        gdp->getP()->bumpDataId();

    return error();
}
