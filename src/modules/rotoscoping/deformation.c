/*
 * deformation.c
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

#include "deformation.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


inline PointF vVec( double x, double y )
{
    PointF v;
    v.x = x;
    v.y = y;
    return v;
}

inline PointF *vVecD( double x, double y, PointF *v )
{
    v->x = x;
    v->y = y;
    return v;
}

inline int vIsEqual( PointF *v1, PointF *v2 )
{
    return v1->x == v2->x && v1->y == v2->y;
}

static PointF *vSub( PointF *v1, PointF *v2, PointF *v )
{
    v->x = v1->x - v2->x;
    v->y = v1->y - v2->y;
    return v;
}

static PointF *vAdd( PointF *v1, PointF *v2, PointF *v )
{
    v->x = v1->x + v2->x;
    v->y = v1->y + v2->y;
    return v;
}

static void vAddD( PointF *v, PointF *v1 )
{
    v->x += v1->x;
    v->y += v1->y;
}

static PointF *vFact( PointF *v1, double f, PointF *v )
{
    v->x = v1->x * f;
    v->y = v1->y * f;
    return v;
}

static void vFactD( PointF *v, double f )
{
    v->x *= f;
    v->y *= f;
}

static double vLength( PointF *v )
{
    return sqrt( v->x * v->x + v->y * v->y );
}

static double vLengthSquared( PointF *v )
{
    return v->x * v->x + v->y * v->y;
}

void deform( PointF *p, PointF *q, int count, PointF v, double alpha, PointF *result )
{
    if ( !count )
    {
        *result = v;
        return;
    }

    int i;
    double *w = malloc( count * sizeof( double ) );
    double wSum = 0;
    PointF pS = vVec( 0, 0 );                              // p*
    PointF qS = vVec( 0, 0 );                              // q*
    PointF vSubPS;                                         // v - p*
    PointF tmp;
    

    for ( i = 0; i < count; ++i )
    {
        if ( vIsEqual( &p[i], &v ) )
        {
            *result = q[i];
            free( w );
            return;
        }

        w[i] = 1 / pow( vLengthSquared( vSub( &p[i], &v, &tmp ) ), alpha );
        wSum += w[i];
        vAddD( &pS, vFact( &p[i], w[i], &tmp ) );
        vAddD( &qS, vFact( &q[i], w[i], &tmp ) );
    }

    vFactD( &pS, 1 / wSum );
    vFactD( &qS, 1 / wSum );

    vSub( &v, &pS, &vSubPS );

    PointF fr = vVec( 0, 0 );                              // fr(v)
    for ( i = 0; i < count; ++i )
    {
        PointF pSR, qSR;
        vSub( &p[i], &pS, &pSR );
        vSub( &q[i], &qS, &qSR );
        vFactD( &qSR, w[i] );

        vAddD( &fr, vVecD( qSR.x * pSR.x + qSR.y * pSR.y, qSR.x * pSR.y - qSR.y * pSR.x, &tmp ) );
    }

    vFactD( &fr, 1 / vLength( &fr ) );
    vAdd( vVecD( fr.x * vSubPS.x + fr.y * vSubPS.y, fr.x * vSubPS.y - fr.y * vSubPS.x, &tmp ) , &qS, result );

    free( w );
}
