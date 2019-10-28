#pragma once

#include "curves.hh"

// Leaves the first and last `fix' number of control points at their original positions.
// Assumes that the directions are unit vectors.
void optimize(BSplineCurve &curve, const VectorVector &directions, size_t n_samples, size_t fix);
