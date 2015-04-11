


#ifndef __SOP_RBFDeform_h__
#define __SOP_RBFDeform_h__

#include <SOP/SOP_Node.h>
#include "interpolation.h"

namespace RBFDeform {
/// Run a sin() wave through geometry by deforming points
/// @see @ref HOM/SOP_HOMWave.C, SOP_HOMWave, SOP_CPPWave
class SOP_RBFDeform : public SOP_Node
{
public:
	     SOP_RBFDeform(OP_Network *net, const char *name, OP_Operator *op);
    virtual ~SOP_RBFDeform();

    static PRM_Template		 myTemplateList[];
    static OP_Node		*myConstructor(OP_Network*, const char *,
							    OP_Operator *);

    /// This method is created so that it can be called by handles.  It only
    /// cooks the input group of this SOP.  The geometry in this group is
    /// the only geometry manipulated by this SOP.
    virtual OP_ERROR		 cookInputGroups(OP_Context &context, 
						int alone = 0);

protected:
    /// Method to cook geometry for the SOP
    virtual OP_ERROR		 cookMySop(OP_Context &context);

private:
    void	getGroups(UT_String &str){ evalString(str, "group", 0, 0); }
    void    MODEL(UT_String &str)    { evalString(str, "model", 0, 0); }
    void    TERM(UT_String &str)     { evalString(str, "term", 0, 0); }
    fpreal	RADIUS(fpreal t)	{ return evalFloat("radius", 0, t); }
    int  	LAYERS(fpreal t)	{ return evalInt("layers", 0, t); }
    fpreal	LAMBDA(fpreal t)	{ return evalFloat("lambda", 0, t); }

    /// This is the group of geometry to be manipulated by this SOP and cooked
    /// by the method "cookInputGroups".
    const GA_PointGroup *myGroup;
};

class op_RBFDeform {
public:
    op_RBFDeform(const GA_RWAttributeRef &p; const alglib::rbfmodel &model)
        : myP(p), myModel(model) {}
            // Take a SplittableRange (not a GA_Range)
    void    operator()(const GA_SplittableRange &r) const
            {
                GA_RWPageHandleV3 p_ph(myP.getAttribute());
                // Execute storage:
                alglib::real_1d_array coord;
                alglib::real_1d_array result;
                coord.setlength(3);
                result.setlength(3);
                // Iterate over pages in the range
                for (GA_PageIterator pit = r.beginPages(); !pit.atEnd(); ++pit)
                {
                    GA_Offset start, end;
                    // iterate over the elements in the page.
                    for (GA_Iterator it(pit.begin()); it.blockAdvance(start, end); )
                    {
                        // Perform any per-page setup required, then
                        p_ph.setPage(start);
                        for (GA_Offset i = start; i < end; ++i)
                        {
                            UT_Vector3      pos = p_ph.get(i);
                            // N.normalize();
                            // v_ph.set(i, N);
                        }
                    }
                }
            }
    private:
        const GA_RWAttributeRef myP;
        const alglib::rbfmodel myModel;
};
void
rbfDeform(const GA_Range &range, const GA_RWAttributeRef &p, const alglib::rbfmodel &model)
{
    // Create a GA_SplittableRange from the original range
    UTparallelFor(GA_SplittableRange(range), op_RBFDeform(p, model));
}




} // End HDK_Sample namespace

#endif
