#include "hypergeometric.hh"

#include <cmath>
#include <stdlib.h>

#include <python2.7/Python.h>


typedef struct{
	double Rp, xmax;
} paraHG;


int intHypGeo(const int *ndim, const double x[], const int *ncomp, double f[], void *params){
    paraHG *pa = (paraHG*)params;
    
    double Rp = pa->Rp;
    double xmax = pa->xmax;
    
    double y = x[0]*xmax;
    
    f[0] = (pow(1.-y, Rp)-1.)/y*xmax;
    return(0);
}

double HypGeo_Py(double Rp, double xmax){
  //PyObject *pName, *pModule, *pDict, *pFunc;
  PyObject *pName, *pModule, *pFunc;
    PyObject *pArgs, *pValue;
    double result = 0.0;
    
    Py_Initialize();
    pName = PyString_FromString("mpmath");
    /* Error checking of pName left out */

    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, "hyp3f2");
        /* pFunc is a new reference */

        if (pFunc && PyCallable_Check(pFunc)) {
            pArgs = PyTuple_New(6);
            pValue = PyFloat_FromDouble(1.);
            if (!pValue) {
                Py_DECREF(pArgs);
                Py_DECREF(pModule);
                fprintf(stderr, "Cannot convert argument\n");
                return 1;
            }
            /* pValue reference stolen here: */
            PyTuple_SetItem(pArgs, 0, pValue);
            PyTuple_SetItem(pArgs, 1, pValue);
            
            pValue = PyFloat_FromDouble(1.-Rp);
            if (!pValue) {
                Py_DECREF(pArgs);
                Py_DECREF(pModule);
                fprintf(stderr, "Cannot convert argument\n");
                return 1;
            }
            PyTuple_SetItem(pArgs, 2, pValue);
            
            pValue = PyFloat_FromDouble(2.);
            if (!pValue) {
                Py_DECREF(pArgs);
                Py_DECREF(pModule);
                fprintf(stderr, "Cannot convert argument\n");
                return 1;
            }
            /* pValue reference stolen here: */
            PyTuple_SetItem(pArgs, 3, pValue);
            PyTuple_SetItem(pArgs, 4, pValue);
            
            pValue = PyFloat_FromDouble(xmax);
            if (!pValue) {
                Py_DECREF(pArgs);
                Py_DECREF(pModule);
                fprintf(stderr, "Cannot convert argument\n");
                return 1;
            }
            /* pValue reference stolen here: */
            PyTuple_SetItem(pArgs, 5, pValue);
            
            pValue = PyObject_CallObject(pFunc, pArgs);
            Py_DECREF(pArgs);
            if (pValue != NULL) {
                result = PyFloat_AsDouble(pValue);
                Py_DECREF(pValue);
            }
            else {
                Py_DECREF(pFunc);
                Py_DECREF(pModule);
                PyErr_Print();
                fprintf(stderr,"Call failed\n");
                return 1;
            }
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", "hyp3f2");
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", "mpmath");
        return 1;
    }
    Py_Finalize();
    return(-Rp*xmax*result);
}


