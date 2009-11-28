// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef __HERMES1D_MATRIX_H
#define __HERMES1D_MATRIX_H

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <math.h>
#include <string.h>

#include "common.h"

/// Creates a new (full) matrix with m rows and n columns with entries of the type T.
/// The entries can be accessed by matrix[i][j]. To delete the matrix, just
/// do "delete matrix".
template<typename T>
T** new_matrix(int m, int n = 0)
{
  if (!n) n = m;
  T** vec = (T**) new char[sizeof(T*)*m + sizeof(T)*m*n];
  if (vec == NULL) error("Out of memory.");
  T* row = (T*) (vec + m);
  for (int i = 0; i < m; i++, row += n)
    vec[i] = row;
  return vec;
}


/// Transposes an m by n matrix. If m != n, the array matrix in fact has to be
/// a square matrix of the size max(m, n) in order for the transpose to fit inside it.
template<typename T>
void transpose(T** matrix, int m, int n)
{
  int min = std::min(m, n);
  for (int i = 0; i < min; i++)
    for (int j = i+1; j < min; j++)
      std::swap(matrix[i][j], matrix[j][i]);
    
  if (m < n)
    for (int i = 0; i < m; i++)
      for (int j = m; j < n; j++)
        matrix[j][i] = matrix[i][j];
  else if (n < m)
    for (int i = n; i < m; i++)
      for (int j = 0; j < n; j++)
        matrix[j][i] = matrix[i][j];
}


/// Changes the sign of a matrix
template<typename T>
void chsgn(T** matrix, int m, int n)
{
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      matrix[i][j] = -matrix[i][j];
}


/// Given a matrix a[n][n], this routine replaces it by the LU decomposition of a rowwise
/// permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
/// indx[n] is an output vector that records the row permutation effected by the partial
/// pivoting; d is output as +-1 depending on whether the number of row interchanges was even
/// or odd, respectively. This routine is used in combination with lubksb to solve linear equations
/// or invert a matrix.
void ludcmp(double** a, int n, int* indx, double* d);

/// Solves the set of n linear equations AX = B. Here a[n][n] is input, not as the matrix
/// A but rather as its LU decomposition, determined by the routine ludcmp. indx[n] is input
/// as the permutation vector returned by ludcmp. b[n] is input as the right-hand side vector
/// B, and returns with the solution vector X. a, n, and indx are not modified by this routine
/// and can be left in place for successive calls with different right-hand sides b. This routine takes
/// into account the possibility that b will begin with many zero elements, so it is efficient for use
/// in matrix inversion.
template<typename T>
void lubksb(double** a, int n, int* indx, T* b)
{
  int i, ip, j;
  T sum;

  for (i = 0; i < n; i++)
  {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    for (j = 0; j < i; j++) sum -= a[i][j]*b[j];
    b[i] = sum;
  }
  for (i = n-1; i >= 0; i--)
  {
    sum = b[i];
    for (j = i+1; j < n; j++) sum -= a[i][j]*b[j];
    b[i] = sum / a[i][i];
  }
}




/// Given a positive-definite symmetric matrix a[n][n], this routine constructs its Cholesky
/// decomposition, A = L*L^T . On input, only the upper triangle of a need be given; it is not
/// modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
/// elements which are returned in p[n].
void choldc(double **a, int n, double p[]);

/// Solves the set of n linear equations A*x = b, where a is a positive-definite symmetric matrix.
/// a[n][n] and p[n] are input as the output of the routine choldc. Only the lower
/// subdiagonal portion of a is accessed. b[n] is input as the right-hand side vector. The
/// solution vector is returned in x[n]. a, n, and p are not modified and can be left in place
/// for successive calls with different right-hand sides b. b is not modified unless you identify b and
/// x in the calling sequence, which is allowed. The right-hand side b can be complex, in which case
/// the solution x is also complex.
template<typename T>
void cholsl(double **a, int n, double p[], T b[], T x[])
{
  int i, k;
  T sum;

  for (i = 0; i < n; i++)
  {
    sum = b[i];
    k = i;
    while (--k >= 0)
      sum -= a[i][k] * x[k];
    x[i] = sum / p[i];
  }

  for (i = n-1; i >= 0; i--)
  {
    sum = x[i];
    k = i;
    while (++k < n)
      sum -= a[k][i] * x[k];
    x[i] = sum / p[i];
  }
}

class Matrix {
public:
    virtual ~Matrix() { }
    virtual int get_size() = 0;
    virtual void add(int m, int n, double v) = 0;
    virtual double get(int m, int n) = 0;
    virtual void zero() = 0;
    virtual void copy_into(Matrix *m) = 0;
    virtual void print() = 0;
};

class Triple {
    public:
        int i;
        int j;
        double v;
        Triple *next;

        Triple(int i, int j, double v) {
            this->i = i;
            this->j = j;
            this->v = v;
            this->next = NULL;
        }
};

class CooMatrix : public Matrix {
    public:
        CooMatrix(int size) {
            this->size = size;
            this->list = NULL;
            this->list_last = NULL;
            this->zero();
        }
        ~CooMatrix() {
            this->free_data();
        }
        void free_data() {
            Triple *t = this->list;
            while (t != NULL) {
                Triple *t_old = t;
                t = t->next;
                delete t_old;
            }
        }
        virtual void zero() {
            if (this->list != NULL)
                this->free_data();
            this->list = NULL;
            this->list_last = NULL;
        }
        virtual void add(int m, int n, double v) {
            Triple *t = new Triple(m, n, v);
            if (this->list == NULL) {
                this->list = t;
                this->list_last = this->list;
            } else {
                this->list_last->next = t;
                this->list_last = this->list_last->next;
            }
            if (m > this->size-1) {
              printf("m = %d, n = %d, this->size = %d\n", m, n, this->size);
              error("m is bigger than size");
            }
            if (n > this->size-1) {
              printf("m = %d, n = %d, this->size = %d\n", m, n, this->size);
              error("n is bigger than size");
            }
        }
        virtual double get(int m, int n) {
            double v=0;
            Triple *t = this->list;
            if (t != NULL) {
                while (t->next != NULL) {
                    if (m == t->i && n == t->j)
                        v += t->v;
                    t = t->next;
                }
            }
            return v;
        }

        virtual int get_size() {
            return this->size;
        }

        virtual void copy_into(Matrix *m) {
            m->zero();
            Triple *t = this->list;
            while (t != NULL) {
                m->add(t->i, t->j, t->v);
                t = t->next;
            }
        }

        virtual void print() {
            Triple *t = this->list;
            while (t != NULL) {
                printf("(%d, %d): %f\n", t->i, t->j, t->v);
                t = t->next;
            }
        }

    private:
        int size;
        /*
           We represent the COO matrix as a list of Triples (i, j, v), where
           (i, j) can be redundant (then the corresponding "v" have to be
           summed). We remember the pointers to the first (*list) and last
           (*list_last) element.
           */
        Triple *list;
        Triple *list_last;
};

class DenseMatrix : public Matrix {
    public:
        DenseMatrix(int size) {
            this->mat = new_matrix<double>(size, size);
            this->size = size;
            this->zero();
        }
        DenseMatrix(Matrix *m) {
            this->mat = new_matrix<double>(m->get_size(), m->get_size());
            //printf("%d %d", this->size, m->get_size());
            //exit(1);
            this->size = m->get_size();
            //this->size = size;
            m->copy_into(this);
        }
        ~DenseMatrix() {
            delete[] this->mat;
        }
        virtual void zero() {
            // erase matrix
            for(int i = 0; i < this->size; i++)
                for(int j = 0; j < this->size; j++)
                    this->mat[i][j] = 0;
        }
        virtual void add(int m, int n, double v) {
            this->mat[m][n] += v;
            //printf("calling add: %d %d %f\n", m, n, v);
        }
        virtual double get(int m, int n) {
            return this->mat[m][n];
        }

        virtual int get_size() {
            return this->size;
        }
        virtual void copy_into(Matrix *m) {
            m->zero();
            for (int i = 0; i < this->size; i++)
                for (int j = 0; j < this->size; j++) {
                    double v = this->get(i, j);
                    if (fabs(v) > 1e-12)
                        m->add(i, j, v);
                }
        }

        virtual void print() {
            for (int i = 0; i < this->size; i++) {
                for (int j = 0; j < this->size; j++) {
                    double v = this->get(i, j);
                    printf("%f ", v);
                    /*
                    if (fabs(v) > 1e-12)
                        printf("%f ", v);
                    else
                        printf("0 ");
                        */
                }
                printf("\n");
            }
        }

        // Return the internal matrix.
        double **get_mat() {
            return this->mat;
        }

    private:
        int size;
        double **mat;

};

class CSRMatrix : public Matrix {
    public:
        CSRMatrix(CooMatrix *m) {
            DenseMatrix *dmat = new DenseMatrix(m);
            this->copy_from_dense_matrix(dmat);
            delete dmat;
        }
        CSRMatrix(DenseMatrix *m) {
            this->copy_from_dense_matrix(m);
        }
        ~CSRMatrix() {
            delete[] this->A;
            delete[] this->IA;
            delete[] this->JA;
        }

        void copy_from_dense_matrix(DenseMatrix *m) {
            this->size = m->get_size();
            this->nnz = 0;
            for(int i = 0; i < this->size; i++)
                for(int j = 0; j < this->size; j++) {
                    double v = m->get(i, j);
                    if (fabs(v) > 1e-12)
                        this->nnz++;
                }
            this->A = new double[this->nnz];
            this->IA = new int[this->size+1];
            this->JA = new int[this->nnz];
            int count = 0;
            this->IA[0] = 0;
            for(int i = 0; i < this->size; i++) {
                for(int j = 0; j < this->size; j++) {
                    double v = m->get(i, j);
                    if (fabs(v) > 1e-12) {
                        this->A[count] = v;
                        this->JA[count] = j;
                        count++;
                    }
                }
                this->IA[i+1] = count;
            }
        }

        virtual void zero() {
            error("Not implemented.");
        }
        virtual void add(int m, int n, double v) {
            error("Not implemented.");
        }
        virtual double get(int m, int n) {
            error("Not implemented.");
        }

        virtual int get_size() {
            return this->size;
        }
        virtual void copy_into(Matrix *m) {
            error("Not implemented.");
        }

        virtual void print() {
            printf("(I, J): value:\n");
            for (int i = 0; i < this->nnz; i++) {
                printf("(%d, %d): %f\n", -1,
                        this->JA[i], this->A[i]);
            }
            printf("IA:\n");
            for (int i = 0; i < this->size+1; i++) {
                printf("%d ", this->IA[i]);
            }
            printf("\n");
        }

        int *get_IA() {
            return this->IA;
        }
        int *get_JA() {
            return this->JA;
        }
        double *get_A() {
            return this->A;
        }

    private:
        int size;
        int nnz;
        double *A;
        int *IA;
        int *JA;

};

// solve linear system
void solve_linear_system(Matrix *mat, double *res);
void solve_linear_system_dense(DenseMatrix *mat, double *res);

#endif
