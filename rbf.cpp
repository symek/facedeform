#include "interpolation.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>


#include <GU/GU_Detail.h>
#include <UT/UT_Matrix.h>
#include <UT/UT_MatrixSolver.h>
#include <GEO/GEO_PointTree.h>


#define ALGLIB_MODEL_QNN 0
#define ALGLIB_MODEL_ML   1

#define ALGLIB_TERM_LINEAR 0
#define ALGLIB_TERM_CONST  1
#define ALGLIB_TERM_ZERO   2
// 
#define CUSTOM_RBF


class Timer 
{
    private:
        double begTime;
    public:
        void start()
        {
            begTime = clock();
        }

        double current() 
        {
        //int threads = UT_Thread::getNumProcessors();
            return (difftime(clock(), begTime) / CLOCKS_PER_SEC);// / (double) threads;
        }

        bool isTimeout(double seconds)
        {
            return seconds >= current();
        }
};

inline float norm(const UT_Vector3 a, const UT_Vector3 b)
{
    const UT_Vector3 s(b-a);
    return s.length();
}

 enum DistanceFunction
 {
    linear,
    cubic,
    quintic,
    thinPlate,
    qaussian,
    inverseMult,
    multiQuadratic
 };

inline float distance(const float r, const DistanceFunction type, const float epsilon)
{
    switch(type)
    {

        case linear:
            return r;
            break;
        case cubic:
            return SYSpow(r, 3);
            break;
        case quintic:
            return SYSpow(r, 5);
            break;
        case thinPlate:
            if (r == 0 ) return 0;
            return SYSpow(r, 2) * SYSlog(r);
            break;
        case qaussian:
            return SYSexp(-SYSpow(r/epsilon, 2));
            break;
        case inverseMult:
            return 1.0 / SYSsqrt(SYSpow(r / epsilon, 2) + 1.0);
            break;
        case multiQuadratic:
            return SYSsqrt(SYSpow(r/epsilon, 2) + 1);
            break;
    }
    return 0.0;
}

inline float interpolant(const GU_Detail *gdp, GEO_PointTree &tree, 
    const UT_VectorD &values, const UT_Vector3 &pos, const DistanceFunction type, 
    const UT_VectorD &epsilons) 
{
    const int npoints = gdp->getNumPoints();
    UT_VectorD dist(1, npoints); // 1-based index for LU comp;
    // GEO_PointTree::IdxArrayType  indexList;
    // UT_FloatArray                groupDist;

    for (int i=0; i<npoints; ++i)
    {
        const UT_Vector3 point = gdp->getPos3(GA_Offset(i));
        // int found = tree.findNearestGroupIdx(point, epsilon*2.0, 2, indexList, groupDist);
        // std::cout << groupDist(1) << ", ";
        dist(i+1) = distance(norm(point, pos), type, epsilons(i+1)); 
    }

    UT_VectorD sum(dist);
    sum *= values;
    float result = 0.0;

    for (int i=0; i<npoints; ++i)
        result += static_cast<float>(sum(i+1));
    
    return result;
}

class RbfDeformModel
{
public:
    int build(const GU_Detail &restGeom, const GU_Detail &deformed, 
        const DistanceFunction distance_func)
    {
        points  = &restGeom;
        npoints = restGeom.getNumPoints();
        M       = UT_Matrix(1, npoints, 1, npoints);
        P       = UT_Permutation(1, npoints);
        values  = UT_ValArray<UT_VectorD>(0, 2); 
        UT_VectorD xs = UT_VectorD(1, npoints);
        UT_VectorD ys = UT_VectorD(1, npoints);
        UT_VectorD zs = UT_VectorD(1, npoints);
        values(0) = xs; values(1) = ys; values(2) = zs;
        distances = UT_VectorD(1, npoints);
        func      = distance_func;

        for(int j=0; j<npoints; ++j)
        {   
            const UT_Vector3 b = restGeom.getPos3(GA_Offset(j));
            const UT_Vector3 c = deformed.getPos3(GA_Offset(j));
            const UT_Vector3 d(c - b);
            values(j+1)(0) = static_cast<double>(d.x());
            values(j+1)(1) = static_cast<double>(d.y());
            values(j+1)(2) = static_cast<double>(d.z());

            for(int i=0; i<npoints; ++i)
            {
                const UT_Vector3 a = restGeom.getPos3(GA_Offset(i));
                const float n      = norm(a, b);
                M.row(j+1)[i+1]    = static_cast<double>(n); // float to double coversion
                epsilon           += n;
            }
        }
        // Compute epsilon: 
        epsilon /= (npoints*npoints) - npoints; 
        // kdtree->build(samples);

        for(int j=0; j<npoints; ++j) {   
            for(int i=0; i<npoints; ++i) {
                const double m = M.row(j+1)[i+1]; 
                M.row(j+1)[i+1] = distance(m, func, epsilon);
            }
        }

        double determinant_sign = 0;
        if( solver.LUDecomp(M, P, determinant_sign) != 1 )
        {
            // Failed to create a LU-decomposition
            return 0;
        }
   
        solver.LUBackSub(M, P, values(0));
        solver.LUBackSub(M, P, values(1));
        solver.LUBackSub(M, P, values(2));


        return 1;
    }
    void interpolate(const UT_Vector3 &pos, UT_Vector3 &result)
    {
        for (int i=0; i<npoints; ++i)
        {
            const UT_Vector3 point = points->getPos3(GA_Offset(i));
            distances(i+1) = distance(norm(point, pos), func, epsilon); 
        }

        UT_VectorD xsum(distances);
        UT_VectorD ysum(distances);
        UT_VectorD zsum(distances);
        xsum *= values(0); 
        ysum *= values(1);
        zsum *= values(2);

       
        for (int i=0; i<npoints; ++i)
        {
            result.x() += static_cast<float>(xsum(i+1)); 
            result.y() += static_cast<float>(ysum(i+1)); 
            result.z() += static_cast<float>(zsum(i+1)); 
        }
    }

private:
    UT_ValArray\
    <UT_VectorD>    values;
    UT_VectorD      distances;
    UT_MatrixSolver solver;
    UT_Matrix       M;
    UT_Permutation  P;
    GEO_PointTree   kdtree;
    DistanceFunction func;
    const GU_Detail *points;
    float           epsilon;
    int             npoints;
};


int main(int argc, char **argv)
{
    Timer t;
    t.start();

    DistanceFunction type = linear; 
    float epsilon_mult    = 1.0;

    if (argc >= 2)
        type = DistanceFunction(atoi(argv[1]));
    if (argc >= 3)
        epsilon_mult = atof(argv[2]);
    


    GA_Offset ptoff;
    GU_Detail *samples = new GU_Detail();
    samples->load("./samples.bgeo");
    GU_Detail * values = new GU_Detail();
    values->load("./values.bgeo");
    int npoints = samples->getNumPoints();

    /* --------- Custom RBF part ---------------- */
    #ifdef CUSTOM_RBF


    #if 0

    RbfDeformModel model;
    model.build(*samples, *samples, type);

    #else



    float epsilon = 0.0f; 
    UT_Matrix M(1, npoints, 1, npoints); // double array
    UT_VectorD v(1, npoints);
    UT_VectorD nearest(1, npoints);
    GEO_PointTree kdtree;
    kdtree.build(samples);
    GEO_PointTree::IdxArrayType  indexList;
    UT_FloatArray                groupDist;


    for(int j=0; j<npoints; ++j)
    {   
        UT_Vector3 b = samples->getPos3(GA_Offset(j));
        int found    = kdtree.findNearestGroupIdx(b, 20.0, 2, indexList, groupDist);
        nearest(j+1) = groupDist(1) * epsilon_mult;
        v(j+1)       = static_cast<double>(b.y());
        for(int i=0; i<npoints; ++i)
        {
            UT_Vector3 a = samples->getPos3(GA_Offset(i));
            a.y() = 0.0f; // zero out as this are values, not coordintes
            b.y() = 0.0f; // as above.
            const float n   = norm(a, b);
            M.row(j+1)[i+1] = static_cast<double>(n); // float to double coversion
            epsilon        += n;
        }
    }

    // Compute epsilon: 
    epsilon /= (npoints*npoints) - npoints; 
    epsilon *= epsilon_mult;




    for(int j=0; j<npoints; ++j) 
    {   
        for(int i=0; i<npoints; ++i) 
        {
            // const UT_Vector3 point = samples->getPos3(GA_Offset(i));
            const double m = M.row(j+1)[i+1]; 
            M.row(j+1)[i+1] = distance(m, type, nearest(i+1)) ;
        }
    }

    UT_MatrixSolver solver;
    UT_Permutation P(1, npoints);
    double determinant_sign = 0;

    if( solver.LUDecomp(M, P, determinant_sign) != 1 )
    {
        // Failed to create a LU-decomposition
        std::cout << "Cant create LU_decomposition" << std::endl;
        delete samples;
        delete values;
        return 1;
        
    }
   
    std::cout << "Before Solve" << std::endl;
    for (int i=0; i<npoints;++i) std::cout << v(i+1) << ", "; std::cout << std::endl;
    
    solver.LUBackSub(M, P, v);

    std::cout << "after solve" << std::endl;
    for (int i=0; i<npoints;++i) std::cout << v(i+1) << ", "; std::cout << std::endl;

    std::cout << "Vector has: " << v.length() << std::endl;
    std::cout << "Npoints is: " << npoints   << std::endl;
    std::cout << "Epsilon is: "  << epsilon   << std::endl;


    /* Interpolation of target*/
    GA_FOR_ALL_PTOFF(values, ptoff)
    {
        const UT_Vector3 pos = values->getPos3(ptoff);
        const float result   = interpolant(samples, kdtree, v, pos, type, nearest);
        UT_Vector3 displace  = UT_Vector3(pos.x(), result, pos.z());
        values->setPos3(ptoff, displace);
    }

    #endif

    /* ------------ alglib part -------------- */
    #else

    alglib::real_2d_array rbf_data_model;
    rbf_data_model.setlength(npoints, 3);

    {   
        GA_FOR_ALL_PTOFF(samples, ptoff)
        {
            const UT_Vector3 pos   = samples->getPos3(ptoff);
            const double data[3]   = {pos.x(), pos.z(), pos.y()};
            for (int i=0; i<3; ++i)
                rbf_data_model[static_cast<uint>(ptoff)][i] = data[i];   
        }
    }

    std::string str_model;
    alglib::rbfmodel model;
    alglib::rbfreport report;
    try 
    {
        alglib::rbfcreate(2, 1, model);
    }
    catch (alglib::ap_error err)
    {
        std::cout << "Exception: " << err.msg << std::endl;
    }
    alglib::rbfsetpoints(model, rbf_data_model);
    alglib::rbfsetalgoqnn(model, 1.0, 5.0);
    alglib::rbfsetlinterm(model);

    // Finally build model:
    alglib::rbfbuildmodel(model, report);
    std::cout << static_cast<int>(report.terminationtype) << ", " <<  \
            static_cast<int>(report.iterationscount) << std::endl;

    // Execute storage:
    alglib::real_1d_array coord("[0,0]");
    alglib::real_1d_array result("[0]");

    GA_FOR_ALL_PTOFF(values, ptoff)
    {
        const UT_Vector3 pos = values->getPos3(ptoff);
        const double dp[2]   = {pos.x(), pos.z()};
        coord.setcontent(2, dp);
        alglib::rbfcalc(model, coord, result);
        UT_Vector3 displace = UT_Vector3(pos.x(), result[0], pos.z());
        values->setPos3(ptoff, displace);
        
    }

    #endif


    values->save("./result.bgeo", 0, 0);
    delete samples;
    delete values;
    std::cout << t.current() << std::endl;
    return 0;
}