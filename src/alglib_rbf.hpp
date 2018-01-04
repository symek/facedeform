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
	typedef std::unique_ptr<alglib::rbfreport>      RbfReport;
public:
	bool init()     {}
	bool built()    {}
	bool is_valid() {}
private:
	bool m_init  = false;
	bool m_built = false;
	RbfModelDataPtr model_data = nullptr; 
	RbfModelPtr     model      = nullptr;
	RbfReport       report     = nullptr;
	char            info_buffer[200];
	alglib::real_1d_array coord("[0,0,0]");
    alglib::real_1d_array result("[0,0,0]");
};


}
