#include "interpolation.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>

#include "rbf.hpp"




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


    #if 1

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
        nearest(j+1) = groupDist(1) * epsilon_mult; // FIXME: this is not always true.
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