#include <memory>
#include <unordered_map>
#include <GU/GU_Detail.h>
#include <GQ/GQ_Detail.h>
#include <GEO/GEO_PointTree.h>
#include "deform.hpp"

namespace facedeform
{
template<typename T>
bool CageDeformer<T>::init(const GU_Detail * mesh, const GU_Detail * rest_rig, const GU_Detail * deform_rig)  {
    // m_interpolator->reset(nullptr);
    // m_interpolator->init(mesh, rest_rig, deform_rig);
    m_init = true;
    m_built =false;
    return m_init;
}
template <typename T>
bool CageDeformer<T>::build(const float & radius) {
    if (!m_init) {
        return false;
    }
    // m_interpolator->build(radius);
    m_built = true;
    return m_built;  
}



} // end of namespace facedeform