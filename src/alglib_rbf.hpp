#pragma once
#include "interpolation.h"

enum ALGLIB_Model {
    ALGLIB_MODEL_QNN,
    ALGLIB_MODEL_ML,
};

enum ALGLIB_Term {
    ALGLIB_TERM_LINEAR,
    ALGLIB_TERM_CONST,
    ALGLIB_TERM_ZERO,
};

class InterpolatorBase {
public:
    virtual bool init()     = 0;
    virtual bool build()    = 0;
    virtual bool is_valid() = 0;
};


class ALGLIB_Interpolator : public InterpolatorBase
{
	typedef std::unique_ptr<alglib::real_2d_array>  RbfModelDataPtr;
	typedef std::unique_ptr<alglib::rbfmodel model> RbfModelPtr;
public:
	bool init()     {}
	bool built()    {}
	bool is_valid() {}
private:
	bool m_init  = false;
	bool m_built = false;
	RbfModelDataPtr model_data = nullptr; 
	RbfModelPtr     model      = nullptr;
};


