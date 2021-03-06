%module (package="nuto") ModuleNuToBase
%feature("autodoc","1");
%{
/* Put headers and other declarations here to be added in the wrapper files */
#include "base/Exception.h"
#include <string>
#include <iostream>
%}

%exceptionclass Exception;
%include "exception.i"

// avoid swig warning 
namespace std
{
class exception
{
};
}
%include "base/Exception.h"

// see Examples/python/exceptproxy in your swig installation for exception handling in swig
%exception
{
    try
    {
        $action
    }
    catch(NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl << std::endl;
        NuTo::Exception *ecopy = e.Clone();  //new NuToException(e);   
        PyObject *err = SWIG_NewPointerObj(&ecopy, SWIGTYPE_p_NuTo__Exception, 1);
        PyErr_SetObject(SWIG_Python_ExceptionType(SWIGTYPE_p_NuTo__Exception), err);
        SWIG_fail;
    }
    SWIG_CATCH_STDEXCEPT // catch std::exception
    catch (...)
    {
        SWIG_exception_fail(SWIG_UnknownError, "Unknown exception");
    }
}
