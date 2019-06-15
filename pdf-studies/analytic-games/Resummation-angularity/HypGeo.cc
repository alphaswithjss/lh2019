#include "Definitions.h"

// #include "cuba.h"
// #include <time.h>
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
    
//     cout << x[0] << endl;
    
    
    f[0] = (pow(1.-y, Rp)-1.)/y*xmax;
    
//     cout << Rp << "\t" << "\t" << y << "\t" << f[0] << endl;
    
    return(0);
}
/*
double HypGeo(double Rp, double xmax){
//     return(0);
    paraHG temp = {Rp, xmax};
    
    
    double result[1], error[1], prob[1];
    int neval, fail, nregions;

    int ndim = 1.;
    int ncomp = 1.;
    double relerr = 1.E-3;
    double abserr = 0.; 
    int flags = 0.;
    srand( time( NULL ) );
    int seed = rand() % 1000 + 1; // generates a random number between 1 to 1000 needed for monte carlo
// 	seed = 979;
//   	seed = 162;
    int mineval = 0.;
    int maxeval = 2.E7;
    char *statefile = NULL;
    int nvec = 1.;
    char *spin = NULL;

///Vegas variables///
    int nstart = 8.E3;
    int nincrease = 6E3;
    int nbatch = 500;
    int gridno = 1.;
    
//     cubacores(1,10000);
    
    
    Vegas(ndim, ncomp, intHypGeo, &temp, nvec, relerr, abserr, flags, seed, mineval, 10.*nstart/10., nstart/10., nincrease/30./100., nbatch, gridno, statefile, spin, &neval, &fail, result, error, prob);
    Vegas(ndim, ncomp, intHypGeo, &temp, nvec, relerr, abserr, flags, seed, mineval, maxeval, nstart, nincrease, nbatch, gridno, statefile, spin, &neval, &fail, result, error, prob);
    
//         cout << Rp << "\t" << "\t" << xmax << "\t" << result[0] << endl;
    
    return(result[0]);
}
*/
double HypGeo_Py(double Rp, double xmax){
    PyObject *pName, *pModule, *pDict, *pFunc;
    PyObject *pArgs, *pValue;
    double result;
    
//     char *cwd = wd.c_str();

//     Py_SetProgramName(strdup(wd.c_str()));
    
//     Py_SetProgramName(wd);
//     cout << wd << endl;
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
//                 printf("Result of call: %f\n", PyFloat_AsDouble(pValue));
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

double HypGeo_Grid(double Rp, double xmax){        
        int n1 = 0;
	int n2 = 0;
	
	int n1mx = 99.;
	int n2mx = 99.;
        
	while(double(n1)/99.<=xmax && n1 < n1mx) n1++;
	
	while(0.2*double(n2)/99.<=Rp && n2 < n2mx) n2++;
        
	if(n1==0) n1++;
	if(n2==0) n2++;

	double x1d = (xmax-(n1-1.)/99.)/(n1/99.-(n1-1.)/99.);
        double x2d = (Rp-0.2*(n2-1.)/99.)/(0.2*n2/99.-0.2*(n2-1.)/99.);
	
	double V00 = Grid[(n1-1)*(n2mx+1)+(n2-1)];
	double V10 = Grid[(n1)*(n2mx+1)+(n2-1)];
	double V01 = Grid[(n1-1)*(n2mx+1)+(n2)];
	double V11 = Grid[(n1)*(n2mx+1)+(n2)];
	
	double V0 = V00*(1.-x2d)+V01*x2d;
	double V1 = V10*(1.-x2d)+V11*x2d;

	double V = V0*(1.-x1d)+V1*x1d;
        
//         cout << V00 << "\t" << V10 << "\t" << V01 << "\t" << V11 << endl;
	/*
        double V_py = HypGeo_Py(Rp, xmax);
        
//         cout << V << "\t" << V_py << "\t" << (1.-V_py/V) << "\t";
        
        if(abs((1.-V_py/V))>1E-2){
            cout << Rp << "\t" << xmax << "\t" << V_py << endl;
            cout << V << "\t" << V_py << "\t" << (1.-V_py/V) << "\t";
            cout << n1 << "\t" << n2 << endl;
            cout << (n1-1)/99. << "\t" << (n1)/99. << "\t" << xmax << endl;
            cout << 0.2*(n2-1)/99. << "\t" << 0.2*(n2)/99. << "\t" << Rp << endl;
            cout << Grid[(n1-1)*(n2mx+1)+(n2-1)] << "\t" << Grid[(n1)*(n2mx+1)+(n2-1)] << endl;
            cout << Grid[(n1-1)*(n2mx+1)+(n2)] << "\t" << Grid[(n1)*(n2mx+1)+(n2)] << endl;
            cout << endl;
            exit(0);
        }*/
        
        
        return(V);
}

