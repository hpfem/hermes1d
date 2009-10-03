#include "discrete.h"

DiscreteProblem::DiscreteProblem(int neq, Mesh *mesh)
{
    this->neq = neq;
    this->mesh = mesh;
}

void DiscreteProblem::add_matrix_form(int i, int j, matrix_form fn)
{
}

void DiscreteProblem::add_vector_form(int i, int j, vector_form fn)
{
}
