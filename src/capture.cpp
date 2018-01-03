#include <memory>
#include <unordered_map>
#include <GQ/GQ_Detail.h>
#include <GEO/GEO_PointTree.h>
#include "capture.hpp"

namespace facedeform
{

bool ProximityCapture::init(GU_Detail * mesh, const GU_Detail * rest_rig)  {

    DEBUG_PRINT("ProximityCapture::init: %i\n", init_counter);
    m_gdp = mesh; m_rig = rest_rig;
    // mesh point tree
    m_gdp_tree.clear();  m_gdp_tree.build(m_gdp);
    // rig point tree
    m_rest_tree.clear(); m_rest_tree.build(m_rig);
    // rig ray intersect
    m_handler_ray_cache.clear(); m_handler_ray_cache.init(m_rig);
    // GQ Edge structure 
    if (m_gq_detail) {
        m_gq_detail->clearAndDestroy();
    }
    m_gq_detail.reset(nullptr); m_gq_detail = std::move(GQ_DetailPtr(new GQ_Detail(m_gdp)));
    // Helper detached attrbis. 
    cd_a.reset(nullptr); dist_a.reset(nullptr);
    // std::map of affected groups
    m_handlers_map.clear();
    // own attribs
    cd_a.reset(m_gdp->createDetachedTupleAttribute(GA_ATTRIB_POINT, GA_STORE_REAL32, 3));
    dist_a.reset(m_gdp->createDetachedTupleAttribute(GA_ATTRIB_POINT, GA_STORE_REAL32, 1));
    // quality check
    if (!cd_a || !dist_a /*more conditions related to trees and ray cache*/) {
        m_init = false; 
        m_capture = false;  
        init_counter = 0;
        capture_counter = 0;  
    }  else {
        init_counter++;
        m_init    = true;
        m_capture = false; 
    }
    return m_init;
}

bool ProximityCapture::capture(const int & max_edges, const float & radius,
    const int & dofalloff, const float & falloffrate ) {

    DEBUG_PRINT("ProximityCapture::capture: %i\n", capture_counter);
    if (!m_init) {
        return false;
    }
    // At least one islads should be found...
    if(findIslands(max_edges) == 0) {
        return false;
    }
    UT_Color affected_clr(UT_RGB, 1.f,1.f,1.f);
    const UT_Vector3 white(1.f, 1.f, 1.f);
    GA_RWHandleV3 cd_h(cd_a.get());
    GA_RWHandleF  fd_dist_h(dist_a.get());

    const float radius_sqrt = radius*radius;
    GU_MinInfo  m_closest_pt_info(radius_sqrt);
    // Iterate other islands per handle.
    HandlerGroupMap::const_iterator it;
    for (it = m_handlers_map.begin(); it != m_handlers_map.end(); it++) {
        const GA_PointGroupPtr & affected_group = it->second;
        for (GA_Size i=0; i<affected_group->entries(); ++i)  {
            const GA_Offset ptoff = affected_group->findOffsetAtGroupIndex(i);
            // Get out from here, if not doing any falloffs
            if (!dofalloff) {
                fd_dist_h.set(ptoff, 0.f);
                cd_h.set(ptoff, white);
                continue;
            }
            float distance_sqrt = -1.f;
            m_closest_pt_info.init(radius_sqrt);
            const UT_Vector3 pos = m_gdp->getPos3(ptoff);
            // We should actually have special case for all points bellow handle
            // to have distance = 0;
            if (m_handler_ray_cache.minimumPoint(pos, m_closest_pt_info)) {
                UT_Vector4 anchor_pos;
                m_closest_pt_info.prim->evaluateInteriorPoint(anchor_pos, \
                    m_closest_pt_info.u1, m_closest_pt_info.v1);
                distance_sqrt = m_closest_pt_info.d;
            }
            // save distance (squered)
            fd_dist_h.set(ptoff, distance_sqrt);
            if (distance_sqrt > radius_sqrt || distance_sqrt < 0) {
                cd_h.set(ptoff, white);
                continue;
            }
            // continue only for distsances smaller than radius
            float falloff = SYSmin(distance_sqrt/radius_sqrt, 1.f);
            falloff = SYSpow(1.f - falloff, falloffrate);
            float hue = SYSfit(falloff, 0.f, 1.f, 200.f, 250.f);
            affected_clr.setHSV(hue, 1.f, 1.f);
            cd_h.set(ptoff, affected_clr.rgb());
        }
    }

    capture_counter++;
    m_capture = true;
    return m_capture;
}

int ProximityCapture::findIslands(const int & max_edges) {

    if (!m_init) {
        return 0;
    }
    GA_PointGroup partial_group(*m_gdp);
    GA_ROHandleI class_h(m_rig->findIntTuple(GA_ATTRIB_POINT, "class", 1));
    if(class_h.isInvalid()) {
        GA_PointGroupPtr handle_group(new GA_PointGroup(*m_gdp));
        m_handlers_map.insert(std::make_pair<int, \
            GA_PointGroupPtr>(0, std::move(handle_group)));
    }
    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(m_rig, ptoff) { 
        const UT_Vector3 anchor_pos  = m_rig->getPos3(ptoff);
        const GA_Index  target_idx   = m_gdp_tree.findNearestIdx(anchor_pos);
        const GA_Offset target_ptoff = m_gdp->pointOffset(target_idx);

        int handle_id = 0;
        if (class_h.isValid()) {
            handle_id = class_h.get(ptoff); 
        }
        if(m_handlers_map.find(handle_id) == m_handlers_map.end()) {
            GA_PointGroupPtr handle_group(new GA_PointGroup(*m_gdp));
            m_handlers_map.insert(std::pair<int, \
                GA_PointGroupPtr>(handle_id, std::move(handle_group)));  
        }
        m_gq_detail->groupEdgePoints(target_ptoff, max_edges, partial_group);
        GA_PointGroup * handle_group = m_handlers_map[handle_id].get();
        handle_group->combine(&partial_group);
        partial_group.clear();
    }

    return m_handlers_map.size();
}


} // end of namespace facedeform