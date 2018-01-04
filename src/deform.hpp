#pragma once 
#include <SOP/SOP_Node.h>

namespace facedeform {

#ifndef NDEBUG
#define DEBUG_PRINT(fmt, ...) fprintf(stderr, fmt, __VA_ARGS__)
#else
#define DEBUG_PRINT(fmt, ...) do {} while (0)
#endif

class DeformerBase
{
public:
    virtual bool init(const GU_Detail *, const GU_Detail *, 
        const GU_Detail *)                     ;// = delete;
    virtual bool is_init()  const              ;// = delete;
    virtual bool build(const float & radius)   ;// = delete;
    virtual bool is_built() const              ;// = delete;
    virtual void * interpolant() const         ;// = delete;
};

struct DummyInterpolator {
    int data = 1234;
};

template <typename T>
class CageDeformer : DeformerBase
{
public: 
        typedef std::unique_ptr<T> InterpolatorPtr;
        CageDeformer<T>() { m_interpolator = std::move(InterpolatorPtr(new T())); }
        virtual bool init(const GU_Detail *, const GU_Detail *, const GU_Detail *);
        virtual bool is_init() const { return m_init; }  
        virtual bool build(const float & radius);
        virtual bool is_built() const { return m_built; }
        virtual T *  interpolant() const { return m_interpolator.get(); }
private:
    bool            m_init  = false;
    bool            m_built = false;
    InterpolatorPtr m_interpolator = nullptr;
};

typedef CageDeformer<DummyInterpolator> DummyDeformer;



} // end of facedeform namespace