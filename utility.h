/*!
 * @file utility.h
 * @brief Contains auxillary functions mainly for Vector operations
 *
 *  Created on: May 4, 2022
 *      Author: Moritz Geßner
 */

#pragma once

#include "global.h"

/*!
 * Calculates the squared Euclidean norm of the Vec @p v
 */
T SqNorm(Vec v);

/*!
 * Calculates the scalar product of the Vecs @p v and @p w
 * the size of @p v must be smaller than that of @p w
 */
T ScalarProd(const Vec& v, const Vec& w);

/*
 * Scale the Vec @p v linearly with @p a
 */
void Scale(Vec& v, T a);

/*!
 * Returns the element-wise sum of two Vecs @p v and @p w
 */
Vec Sum(const Vec& v, const Vec& w);

/*!
 * @brief Adds the components of a vector @p w to @p v elemt-wise and in-place
 * Add the components of @p w to @p v elementwise
 */
void VecAdd(Vec& v, const Vec& w);

/*!
 * Appends @p b to @p a
 */
void AppendVector(Vec& a, const Vec& b);

/*!
 * Calculates the velocity squared, which is proportional to the conserved
 * classical kinetic energy of a particle with the state vector @p q
 * @param q is a state vector \f$ q=(x, y, z, v_x, v_y, v_z)^T \f$
 * \return The classical velocity squared \f$ \vec{v}^2 =  v_x^2 + v_y^2 + v_z^2 \f$
 */
T VSquared(const Vec& q);

/*!
 * Prints the time in format [hh:mm:ss:] without newline
 */
void PrintTime();

/*!
 * Constructs a Unit vector from a given @p zeta with \f$ \zeta \in [-1,1] \f$ and @p phi with \f$ \phi \in [0,2 \pi ) \f$
 */
Vec UnitVec(const T zeta, const T phi);
