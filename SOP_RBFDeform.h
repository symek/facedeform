


#ifndef __SOP_RBFDeform_h__
#define __SOP_RBFDeform_h__

#include <SOP/SOP_Node.h>

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
} // End HDK_Sample namespace

#endif
