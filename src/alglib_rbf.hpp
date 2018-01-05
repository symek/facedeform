#pragma once
#include "interpolation.h"

namespace facedeform {

template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
    return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}

enum ALGLIB_Model {
    ALGLIB_MODEL_QNN,
    ALGLIB_MODEL_ML,
};

enum ALGLIB_Term {
    ALGLIB_TERM_LINEAR,
    ALGLIB_TERM_CONST,
    ALGLIB_TERM_ZERO,
};

struct ALGLIB_Parms {
        ALGLIB_Model model = ALGLIB_MODEL_QNN;
        ALGLIB_Term   term = ALGLIB_TERM_ZERO;
        float        qcoef = 1;
        float        zcoef = 5; 
        float       radius = 1;
        float       layers = 4;
        float       lambda =.1;
};

template<class Geometry_Type, class Parm_Type>
class InterpolatorBase {
public:
    virtual bool init(const Geometry_Type *, 
        const Geometry_Type *)             = 0;
    virtual bool build(const Parm_Type &)  = 0;
    virtual bool is_valid() const          = 0;
    virtual bool interpolate(const double *, UT_Vector3 &) = 0;
};

template<class Geometry_Type, class Parm_Type>
class ALGLIB_RadialBasisFunc_Impl : 
    public InterpolatorBase<Geometry_Type, Parm_Type>
{
    typedef alglib::real_1d_array                   RealArray;
    typedef alglib::real_2d_array                   Real2DArray;
public:
    using InputParametersT = Parm_Type;
    bool init(const Geometry_Type * rest, const Geometry_Type * deform) {
        assert(rest != nullptr);
        assert(deform != nullptr);
        m_rest_geo   = rest;
        m_deform_geo = deform;
        m_data_model.setlength(m_rest_geo->getNumPoints(), 6);
        GA_Offset ptoff; 
        GA_FOR_ALL_PTOFF(rest, ptoff) {
            const UT_Vector3 restP   = rest->getPos3(ptoff);
            const UT_Vector3 deformP = deform->getPos3(ptoff);
            const GA_Index ptidx     = rest->pointIndex(ptoff); 
            const UT_Vector3 delta(deformP - restP);
            double data[6] = {restP.x(), restP.y(), restP.z(),
                delta.x(), delta.y(), delta.z()};
            if (ptidx < m_rest_geo->getNumPoints()) {
                for (int i=0; i<6; ++i)
                    m_data_model[ptidx][i] = data[i]; 
            }
        }
        alglib::rbfcreate(3, 3, m_model);
        alglib::rbfsetpoints(m_model, m_data_model);
        m_init  = true;
        m_built = false;
        return m_init;
    }

    bool build(const Parm_Type & parms) {
        if (!m_init)
            return false;
        // Select RBF model:
        switch(parms.model) {
            case ALGLIB_MODEL_QNN: alglib::rbfsetalgoqnn(m_model,  parms.qcoef, parms.zcoef);
                break;
            case ALGLIB_MODEL_ML: alglib::rbfsetalgomultilayer(m_model, parms.radius,
                parms.layers, parms.lambda);
                break;
        }
        // Select RBF term:
        switch(parms.term) {
            case ALGLIB_TERM_LINEAR: alglib::rbfsetlinterm(m_model);
                break;
            case ALGLIB_TERM_CONST: alglib::rbfsetconstterm(m_model);
                break;
            case ALGLIB_TERM_ZERO: alglib::rbfsetzeroterm(m_model);
                break;
        }
        // Finally build model:
        alglib::rbfbuildmodel(m_model,  m_report);
        if (static_cast<int>(m_report.terminationtype) == 1) {
            m_built = true;
        }
        return m_built;
    }
    bool is_init()  const { return m_init; }
    bool is_built() const { return m_built; }
    bool is_valid() const { return m_init && m_built; }
    bool interpolate( const double * pos, UT_Vector3 & output) {
    	coord.setcontent(3, pos);
        alglib::rbfcalc(m_model, coord, result);
        output.x() = result[0];
        output.y() = result[1];
        output.z() = result[2];
        return true;
    }
private:
    const Geometry_Type * m_rest_geo   = nullptr;
    const Geometry_Type * m_deform_geo = nullptr;
    bool m_init  = false;
    bool m_built = false;
    char               info_buffer[200];
    alglib::rbfmodel   m_model;
    alglib::rbfreport  m_report;
    Real2DArray        m_data_model; 
    RealArray coord  = RealArray("[0,0,0]");
    RealArray result = RealArray("[0,0,0]");
};


typedef ALGLIB_RadialBasisFunc_Impl<GU_Detail, ALGLIB_Parms> ALGLIB_RadialBasisFuncT;
} // end of facedeform namespace
