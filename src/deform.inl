
namespace facedeform
{

template<class Interpol_Type, class Geometry_Type>
bool CageDeformer<Interpol_Type, Geometry_Type>::init(const Geometry_Type *mesh,
    const Geometry_Type *rest_rig, const Geometry_Type *deform_rig)  {
    if (m_init) {
        m_interpolator.reset(nullptr);
        m_interpolator = std::move(InterpolatorPtr(new Interpol_Type()));
        m_init = false;
    }
    m_interpolator->init(rest_rig, deform_rig);
    m_init = true;
    m_built =false;
    return m_init;
}
template <class Interpol_Type, class Geometry_Type>
bool CageDeformer<Interpol_Type, Geometry_Type>::build(const float & radius) {
    if (!m_init) {
        return false;
    }
    typename Interpol_Type::InputParametersT parms;
    if (m_interpolator->build(parms)) {
        m_built = true;
    }
    return m_built;  
}

} // end of namespace facedeform