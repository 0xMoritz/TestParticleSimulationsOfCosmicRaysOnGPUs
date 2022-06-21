/*
 * utility.h
 *
 *  Created on: May 4, 2022
 *      Author: moritz
 */

#pragma once

#include "global.h"


template <typename type>
std::string To_string_nDigits(const type a_value, const int n = 6);

/*!
 * Calculates the squared Euclidean norm of the Vec @v
 */
T SqNorm(Vec v);

/*!
 * calculates the scalar product of the Vecs @v and @w
 * the size of @v must be smaller than that of @w
 */
T ScalarProd(const Vec& v, const Vec& w);

/*
 * Scale the Vec @v linearly with @a
 */
void Scale(Vec& v, T a);

/*!
 * Return the element-wise sum of two Vecs @v and @w
 */
Vec Sum(const Vec& v, const Vec& w);

/*!
 * Add the components of @w to @v elementwise
 */
void VecAdd(Vec& v, const Vec& w);

/*!
 * Appends @b to @a
 */
//template <typename Type>
//void AppendVector(Vec<Type>& a, const Vec<Type>& b);
void AppendVector(Vec& a, const Vec& b);

/*!
* Calculates the velocity squared, which is proportional to the conserved
* classical kinetic energy of a particle with the state vector @q
*/
T VSquared(const Vec& q);

/*!
 * Prints the time in format [hh:mm:ss:] without newline
 */
void PrintTime();

/*!
 * Constructs a Unit vector from a given @eta in [-1,1] and @phi in [0,2\f$ \pi\f$ )
 */
Vec UnitVec(const T eta, const T phi);

