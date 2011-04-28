/*
 * deformation.h
 * Copyright (C) 2011 Till Theato <root@ttill.de>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#ifndef ROTO_DEFORMATION_H
#define ROTO_DEFORMATION_H

typedef struct PointF
{
    double x;
    double y;
} PointF;

inline PointF vVec( double x, double y );
inline PointF *vVecD( double x, double y, PointF *v );

/**
 * Moves \param v by applying a rigid deformation using Moving Least Squares.
 * Algorithm explained and developed in
 * "Image Deformation using Moving Least Squares," Proceedings of ACM SIGGRAPH, pp. 533-540, 2006, Schaefer S., McPhail T., Warren J.
 * http://faculty.cs.tamu.edu/schaefer/research/mls.pdf
 * 
 * \param p List of handles in original position
 * \param q List of handles in new position
 * \param count Number of handles
 * \param v Point on which deformation will be applied on
 * \param alpha Additional parameter controlling interpolation
 * \param result Calculated new point
 */
void deform( PointF *p, PointF *q, int count, PointF v, double alpha, PointF *result );

#endif
