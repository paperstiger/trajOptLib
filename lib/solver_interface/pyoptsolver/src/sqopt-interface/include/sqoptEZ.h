/* Author: Gao Tang
 * An easy-to-use interface for Sqopt
 * It will be a mix of sqoptProblem.hpp and other stuff.
 */
#ifndef SQOPTEZ_H
#define SQOPTEZ_H

#include "sqopt.h"
#include "assert.h"

class sqoptProblem{
public:
  // snopt part
  sqoptProblem();
  sqoptProblem(const char*name);
  ~sqoptProblem();

  void init2zero();

  char    Prob[30];

  int     initCalled, memCalled;

  int     leniw, lenrw;
  double *rw;
  int    *iw;

  int     lenru, leniu;
  double *ru;
  int    *iu;

  void allocI    (int leniw);
  void allocR    (int lenrw);
  void reallocI  (int leniw);
  void reallocR  (int lenrw);

  void setProbName    (const char *Prob);
  void setPrintFile   (const char *prtname);

  int getParameter    (const char *stroptin, char *stroptout);
  int getIntParameter (const char *stropt,   int    &opt);
  int getRealParameter(const char *stropt,   double &opt);
  int setParameter    (const char *stroptin);
  int setIntParameter (const char *stropt,   int     opt);
  int setRealParameter(const char *stropt,   double  opt);

  void setUserI       (int    *iu, int leniu);
  void setUserR       (double *ru, int lenru);
  void setUserspace   (int    *iu, int leniu, double *ru, int lenru);

  // sqopt itself variables
  isqLog sqLog;
  void initialize  (const char *prtfile, int summOn);
  int  setSpecsFile(const char *specname);
  void setLog      (isqLog sqLog);
  void setWorkspace(int m, int n, int neA, int ncObj, int nnH);

  int solve(int starttype, sqFunHx qpHx,
        int m, int n, int neA, int ncObj, int nnH,
        int iObj, double ObjAdd,
        double *A, int *indA, int *locA,
        double *bl, double *bu, double *cObj,
        int *eType, int *hs, double *x, double *pi, double *rc,
        int &nS, int &nInf, double &sInf, double &objective);
};

// copy from Bill code
class SqoptEZ{
public:
    sqoptProblem *qp;
    // n: number of variables
    // m: number of general inequalities, m = numberOfConstraints + 1
    int n, m;
    int inform;

    // H matrix(quadratic cost) in CSC format
    static double* Hv; // values
    static int* Hr; // row indices
    static int* Hc; // column indices

    // c vector(linear cost)
    double* c;
    int lenC; // length of linear cost. If lenC == 0, there is no linear cost

    // A matrix(linear constraint matrix) in CSC format
    double* Av; // values
    int* Ar; // row indices
    int* Ac; // column indices

    // lower bound
    double* lb;

    // upper bound
    double* ub;

    // primal variable and slacks
    // the first n elements are primal variables
    // the rest m elements are slacks
    // see Page 17: https://web.stanford.edu/group/SOL/guides/sqdoc7.pdf
    double* x;

    // elastic indicator
    int* hEtype;

    // sometimes contains a set of states for variables and slacks
    int* hs;

    // dual variables
    double* lambda;

    // vector of reduced cost
    double* rc;

    // objective
    double objective;

    bool initialized;

    static double infBnd;

    // https://stackoverflow.com/questions/2898316/using-a-member-function-pointer-within-a-class
    //typedef int (SqoptEZ::*HxPointer)(int *nnH, double x_[], double Hx[], int *nState,
    //                                 char   cu[], int   *lencu,
    //                                 int    iu[], int   *leniu,
    //                                 double ru[], int   *lenru);
    // HxPointer hxPointer;

    // allocate space for primal and dual variables
    // allocate space for sqoptProblem
    // we do not allocate space for H and A because they are sparse
    void initialize()
    {
        if (m <= 0 || n <= 0)
        {
            printf("Error: Initilization failed due to non positive m or n\n");
            return;
        }

        qp = new sqoptProblem();

        x = new double[n + m];
        hEtype = new int[n + m];
        hs = new int[n + m];
        lambda = new double[n + m];
        rc = new double[n + m];

        int mn = m + n;
        for (int i = 0; i < mn ; ++i)
        {
            x[i] = 0.0;
            hEtype[i] = 0;
            hs[i] = 0;
        }

        initialized = true;
    }

    // copy the first N-element of a vector to another
    void copyVector(const double* from, double** to, const int N)
    {
        if ((*to) != NULL)
            delete[] (*to);

        (*to) = new double[N];

        for (int i = 0; i < N; ++i)
        {
            (*to)[i] = from[i];
        }
    }

    // copy the first N-element of a vector to another
    void copyVector(const int* from, int** to, const int N)
    {
        if ((*to) != NULL)
            delete[] (*to);

        (*to) = new int[N];

        for (int i = 0; i < N; ++i)
        {
            (*to)[i] = from[i];
        }
    }

    // copy a lower bound or upper bound vector
    // will copy the first n elements from the source to destination,
    // then insert a +/-infBnd into the destination depending on the positive argument
    // then copy the rest (m-1) elements from source to destination
    void copyVectorBound(const double* from, double** to, const int n, const int m, const bool positive)
    {
        if ((*to) != NULL)
            delete[] (*to);
        (*to) = new double[n + m];

        int i = 0;
        for (; i < n; ++i)
        {
            (*to)[i] = from[i];
        }

        if (positive)
            (*to)[i] = +infBnd;
        else
            (*to)[i] = -infBnd;
        i++;

        for (; i < n + m; ++i)
        {
            (*to)[i] = from[i - 1];
        }
    }

    // copy a vector of length (m+n) while not copying one of the (n+1)th elements:
    // First will copy the first n elements from the source to destination,
    // then skip the (n+1)th element
    // then copy the rest (m-1) elements from source to destination
    template<typename T>void copyVectorDelete(const T* from, T** to, const int n, const int m)
    {
        // does not allocate space for to if to is not null
        //if ((*to) == NULL)
        //	(*to) = new double[n + m];

        int i = 0;
        for (; i < n; ++i)
        {
            (*to)[i] = from[i];
        }
        for (; i < n + m - 1; ++i)
        {
            (*to)[i] = from[i + 1];
        }
    }

    // copy a CSC matrix from "from" to "to", given the number of columns in the source matrix
    void copyMatrixCSC(const double* fromVal, const int* fromRowIndices, const int* fromColumnIndices,
                       double** toVal,   int** toRowIndices,   int** toColumnIndices,
                       const int numberOfColumns)
    {
        // number of constraints plus 1
        const int nCp1 = numberOfColumns + 1;
        // cout << "ncp1 " <<  nCp1 << endl;

        // number of non-zero elements in fromVal
        const int nnz  = fromColumnIndices[numberOfColumns];
        // cout << "nnz " <<  fromColumnIndices[numberOfColumns] << endl;
        // cout << "nnz " <<  fromColumnIndices[n] << endl;

        // clear the destination vector
        if ((*toVal) != NULL)
            delete [](*toVal);


        if ((*toRowIndices) != NULL)
            delete [](*toRowIndices);

        if ((*toColumnIndices) != NULL)
            delete [](*toColumnIndices);

        // allocate space for destination vector

        (*toVal) = new double[nnz];
        (*toRowIndices) = new int[nnz];
        (*toColumnIndices) = new int [nCp1];

        for (int i = 0; i < nnz; ++i)
        {
            (*toVal)[i] = fromVal[i];
            (*toRowIndices)[i] = fromRowIndices[i];
        }

        for (int i = 0; i < nCp1; ++i)
        {
            (*toColumnIndices)[i] = fromColumnIndices[i];
        }

        // cout << "nnz " <<  fromColumnIndices[numberOfColumns] << endl;
        // cout << "A/M nnz " <<  Hc[n] << endl;

    }

    // @todo may have to take a deeper look at this function
    // @todo
    static void Hx (int *nnH, double x_[], double Hx[], int *nState,
                   char   cu[], int   *lencu,
                   int    iu[], int   *leniu,
                   double ru[], int   *lenru)
    {
        // H is a n-by-n matrix
        for (int i = 0; i < *nnH; ++i)
        {
            Hx[i] = 0;
        }

        for (int j = 0; j < *nnH; ++j) {                      // Iterate over columns
            for (int ptr = Hc[j]; ptr < Hc[j + 1]; ++ptr) // Iterate over rows
            {
                int i = Hr[ptr];
                // H(i, j) is the (i, j) element of matrix H
                // H(i, j) is stored in Hv[ptr]

                // @todo
                if (1)
                {
                    Hx[j] += Hv[ptr] * x_[i];
                }
            }
        }
    }


    public:
    SqoptEZ():
        qp(NULL),
        n(0), m(0),
        inform(-1),
        // Hv(NULL), Hr(NULL), Hc(NULL),
        c(NULL),
        Av(NULL), Ar(NULL), Ac(NULL),
        lb(NULL), ub(NULL),
        x(NULL), hEtype(NULL), hs(NULL),
        lambda(NULL), rc(NULL),
        objective(0),
        initialized(false)
    {
        // hxPointer = &SqoptEZ::Hx;
    }

    /**
     * n_: number of variables
     * m_: number of constraints
     */
    SqoptEZ(int n_, int m_):
        SqoptEZ()
    {
        n = n_;

        // m = numberOfConstraints + 1, the extra one line is for linear objective, which is stacked
        // with linear constraints in Sqopt.
        m = m_ + 1;

        this->initialize();
    }

    // @todo need to implement this thing
    ~SqoptEZ()
    {
        delete[] c; // c vector(linear cost)

        // delete[] Hv;
        // delete[] Hr;
        // delete[] Hc;

        // A matrix(linear constraint matrix) in CSC format
        delete[] Av; // values
        delete[] Ar; // row indices
        delete[] Ac; // column indices

        // lower bound
        delete[] lb;

        // upper bound
        delete[] ub;

        // primal variable and slacks
        // the first n elements are primal variables
        // the rest m elements are slacks
        // see Page 17: https://web.stanford.edu/group/SOL/guides/sqdoc7.pdf
        delete[] x;

        // elastic indicator
        delete[] hEtype;

        // sometimes contains a set of states for variables and slacks
        delete[] hs;

        // dual variables
        delete[] lambda;

        // vector of reduced cost
        delete[] rc;

        delete qp;
    }

    // set matrix H(quadratic cost matrix) or A(linear constraint matrix) in CSC(Column Compressed) format
    bool setMatrixCSC(const char whichMatrix,
                      const double val[], const int rowIndices[], const int columnIndices[])
    {
        assert( initialized == 1 );

        // copy H matrix
        if (whichMatrix == 'H' || whichMatrix == 'h')
        {
            copyMatrixCSC(val, rowIndices, columnIndices, &Hv, &Hr, &Hc, n);
            // cout << "smc, Ac[n] " << Hc[n-1] << endl;
        }
        else if (whichMatrix == 'A' || whichMatrix == 'a')
        {
            copyMatrixCSC(val, rowIndices, columnIndices, &Av, &Ar, &Ac, n);
            // the the first row in matrix A is the linear objective vector,
            // so we push all the linear constraints down one row
            int nnz = Ac[n]; // number of non-zeros in the A
            for (int i = 0; i < nnz; ++i)
            {
                Ar[i]++;
            }
        }
        else
        {
            printf("Error: Unrecognized matrix: %c\n", whichMatrix);
            printf("====>: Only two kinds of matrices can be set\n");
            printf("====>: 'H': linear constraint matrix\n");
            printf("====>: 'A': linear constraint matrix\n");

            return false;
        }

        return true;
    }

    // set a vector to the values provided, the length of the vector will be determined automatically
    bool setVector(const char whichVector, const double vec[])
    {
        assert( initialized == 1 );

        // constant cost vector
        if (tolower(whichVector) == 'c')
        {
            copyVector(vec, &c, n);
            lenC = n;
        }
        // linear cost vector
        else if (tolower(whichVector) == 'i')
        {
            // linear cost vector will be put into the first row of matrix A
            // so matrix A not be null
            if (Av == NULL || Ar == NULL || Ac == NULL)
                return false;

            // @todo what is the difference between constant cost and linear cost vector?
            printf("Error: currently does not support adding linear cost.\n");
            return false;
        }
        // lower bound vector
        else if (tolower(whichVector) == 'l')
        {
            copyVectorBound(vec, &lb, n, m, false);
        }
        // upper bound vector
        else if (tolower(whichVector) == 'u')
        {
            copyVectorBound(vec, &ub, n, m, true);
        }
        // x
        else if (tolower(whichVector) == 'x')
        {
            copyVector(vec, &x, n);
        }
        else
        {
            printf("Error: Unrecognized vector: %c\n", whichVector);
            return false;
        }

        return true;

    }

    // set the linear cost vector to zero
    void setLinearCostToZero()
    {
        assert(initialized == 1);

        if (c != NULL)
            delete []c;

        // c can be any dimension if lenC == 0
        c = new double[1];
        lenC = 0;
    }

    int solve(int startType)
    {
        assert( initialized == 1 );
        // @todo perhapes we do not need to copy all these variables everytime we call solve
        int nS = 0, nInf=0, iObj=0;
        double sInf = 0, objAdd = 0;
        qp->initialize("", 0);
        inform = qp -> solve(startType, Hx, m, n, Ac[n], lenC, n, iObj,
                    objAdd, Av, Ar, Ac, lb, ub, c,
                    hEtype, hs, x, lambda, rc, nS, nInf, sInf, objective);
        return inform;
    }

    double getObjective()
    {
        assert( initialized == 1 );

        return objective;
    }

    void getPrimalSolution(double x_[])
    {
        assert( initialized == 1 );

        for (int i = 0; i < n; ++i)
        {
            x_[i] = x[i];
        }
    }

    void getDualSolution(double lambda_[])
    {
        assert( initialized == 1 );

        // @todo probably it is not rc
        copyVectorDelete(rc, &lambda_, n, m);

    }

    void getHsSolution(int hs_[])
    {
        copyVectorDelete(hs, &hs_, n, m);
    }

};
#endif
