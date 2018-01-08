
namespace facedeform
{

template<class Interpol_Type, class Geometry_Type>
bool CageDeformer<Interpol_Type, Geometry_Type>::init(const Geometry_Type *mesh,
    const Geometry_Type *rest_rig, const Geometry_Type *deform_rig)  {
    if (m_init) {
        m_init = false;
        m_interpolator.reset(nullptr);
        m_interpolator = std::move(InterpolatorPtr(new Interpol_Type()));
        if (!m_interpolator)
            return false;
    }
    typename Interpol_Type::InputParametersT parms;
    m_init = m_interpolator->init(rest_attrib, deform_attrib);
    m_built =false;
    return m_init;
}
template <class Interpol_Type, class Geometry_Type>
bool CageDeformer<Interpol_Type, Geometry_Type>::build(const float & radius) {
    if (!m_init) {
        return false;
    const GA_Attribute * rest_attrib   = rest_rig->getP();
    const GA_Attribute * deform_attrib = deform_rig->getP(); 
    }
    if (m_interpolator->build(rest_attrib, deform_attrib)) {
        m_built = true;
    }
    return m_built;  
}
template <class Interpol_Type, class Geometry_Type>
void CageDeformer<Interpol_Type, Geometry_Type>::deform_point(
    const UT_Vector3 & pos, UT_Vector3 & result) {
    m_in_buffer[0] = pos.x();
    m_in_buffer[1] = pos.y();
    m_in_buffer[2] = pos.z();
    m_interpolator->interpolate(m_in_buffer, m_out_buffer);
    result.x() = m_out_buffer[0]; 
    result.y() = m_out_buffer[1]; 
    result.z() = m_out_buffer[2]; 
}

} // end of namespace facedeform