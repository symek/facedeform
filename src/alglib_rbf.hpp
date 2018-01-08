#pragma once
#include "interpolation.h"
#include <GA/GA_Handle.h>

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
        
        int         in_dim = 3;
        int        out_dim = 3;
        int     model_size = 0; 
};

template<class Parm_Type>
class InterpolatorBase {
public:
    virtual bool init(const Parm_Type & )     = 0;
    virtual bool build(const GA_Attribute *a, 
                       const GA_Attribute *b) = 0;
    virtual bool is_valid() const             = 0;
    virtual void interpolate(const double *, double *) = 0;
};



template<class DataIn_Type, 
         class DataOut_Type, 
         class ModelData_Type 
         >
class DistanceGEO_Model
{
public:
    DistanceGEO_Model(const size_t & in_dim, const size_t & out_dim, const size_t & model_size) {
        // m_data = std::make_unique(in_dim+out_dim);
        m_in_dim = in_dim;
        m_out_dim = out_dim;
        m_model_size = model_size;
    }
    void operator()(const GA_Attribute * a, const GA_Attribute * b, ModelData_Type & output) {
        const GA_ROHandleT<DataIn_Type> h_a(a);
        const GA_ROHandleT<DataIn_Type> h_b(b);
        const GA_IndexMap & m_a = a->getIndexMap();
        const GA_IndexMap & m_b = b->getIndexMap();
        const GA_Index      end = m_a.indexSize();
        assert(end <= m_model_size);
        for (GA_Index idx(0); idx != end; ++idx) {
            const GA_Offset ptoffa = m_a.offsetFromIndex(idx);
            const GA_Offset ptoffb = m_b.offsetFromIndex(idx);
            const DataIn_Type va   = h_a.get(ptoffa);
            const DataIn_Type vb   = h_b.get(ptoffb);
            const DataIn_Type        delta(vb - va);
            for(size_t i=0; i<m_in_dim+m_out_dim; ++i) {
                output[idx][i]          = static_cast<double>(va(i));
                output[idx][m_in_dim+i] = static_cast<double>(vb(i));
            }    
        }
    }
private:
    size_t m_in_dim     = 0;
    size_t m_out_dim    = 0;
    size_t m_model_size = 0;
};


template<class ModelBuilder_Type, 
         class Geometry_Type, 
         class Parm_Type
         >
class ALGLIB_RadialBasisFunc_Impl : 
    public InterpolatorBase<Parm_Type>
{
    typedef alglib::real_1d_array RealArray;
    typedef alglib::real_2d_array Real2DArray;
public:
    using InputParametersT = Parm_Type;
    bool init(const Parm_Type & parms) {
        m_parms = parms;
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
        // Set dimentionality of a problem:
        coord.setlength(parms.in_dim);
        result.setlength(parms.out_dim);
        m_data_model.setlength(parms.model_size, parms.in_dim+parms.out_dim);
        m_init  = true;
        m_built = false;
        return m_init;
    }

    bool build(const GA_Attribute *a, const GA_Attribute *b) {
        if (!m_init)
            return false;
        assert(a != nullptr); assert(b != nullptr);
        alglib::rbfcreate(3, 3, m_model);
        ModelBuilder_Type model_builer(m_parms.in_dim, m_parms.out_dim, 
            m_parms.model_size);
        model_builer(a, b, m_data_model);
        alglib::rbfsetpoints(m_model, m_data_model);
        alglib::rbfbuildmodel(m_model,  m_report);
        if (static_cast<int>(m_report.terminationtype) == 1) {
            m_built = true;
        }
        return m_built;
    }
    bool is_init()  const { return m_init; }
    bool is_built() const { return m_built; }
    bool is_valid() const { return m_init && m_built; }
    void interpolate( const double * pos, double * output) {
    	coord.setcontent(3, pos);
        alglib::rbfcalc(m_model, coord, result);
        output[0] = result[0];
        output[1] = result[1];
        output[2] = result[2];
    }
private:
    Parm_Type m_parms;
    const Geometry_Type * m_rest_geo   = nullptr;
    const Geometry_Type * m_deform_geo = nullptr;
    bool m_init  = false;
    bool m_built = false;
    char               info_buffer[200];
    alglib::rbfmodel   m_model;
    alglib::rbfreport  m_report;
    Real2DArray        m_data_model; 
    RealArray          coord;
    RealArray          result;
};

typedef DistanceGEO_Model<UT_Vector3, UT_Vector3, alglib::real_2d_array> DistanceGEO_ModelT;
typedef ALGLIB_RadialBasisFunc_Impl<DistanceGEO_ModelT, GA_Attribute, ALGLIB_Parms> ALGLIB_RadialBasisFuncT;
} // end of facedeform namespace
