#pragma once 
#include <SOP/SOP_Node.h>

namespace facedeform {

#ifndef NDEBUG
#define DEBUG_PRINT(fmt, ...) fprintf(stderr, fmt, __VA_ARGS__)
#else
#define DEBUG_PRINT(fmt, ...) do {} while (0)
#endif

class ProximityCapture 
{
// public:
    // struct ProximityCaptureParms {};
    typedef std::unique_ptr<GA_PointGroup> GA_PointGroupPtr;
    typedef std::unique_ptr<GQ_Detail>     GQ_DetailPtr;
    typedef std::unordered_map<int, GA_PointGroupPtr> HandlerGroupMap;
public: 
        bool init(GU_Detail * mesh, const GU_Detail * rig);
        bool isInitialized() const { return m_init; }  
        bool capture( const int & max_edges, const float & radius, const int & dofalloff, const float & falloffrate);
        bool isCaptured() const { return m_capture; }  
private:
        int  findIslands(const int & max_edges);
    bool m_init       = false;
    bool m_capture    = false;
    int  init_counter    = 0;
    int  capture_counter = 0;
    GU_Detail       * m_gdp;
    const GU_Detail * m_rig;
    GQ_DetailPtr      m_gq_detail;
    GEO_PointTree     m_gdp_tree;
    GEO_PointTree     m_rest_tree;
    GU_RayIntersect   m_handler_ray_cache;
    HandlerGroupMap   m_handlers_map;
};





} // end of facedeform namespace