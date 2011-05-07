/*
 * spline_handling.c - ...
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


#include <framework/mlt_pool.h>

#include "spline_handling.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


#ifndef MAX
#define MAX( x, y ) ( ( x ) > ( y ) ? ( x ) : ( y ) )
#endif
#ifndef MIN
#define MIN( x, y ) ( ( x ) < ( y ) ? ( x ) : ( y ) )
#endif
#define SQR( x ) ( x ) * ( x )


/** Linear interp */
static inline void lerp( const PointF *a, const PointF *b, PointF *result, double t )
{
    result->x = a->x + ( b->x - a->x ) * t;
    result->y = a->y + ( b->y - a->y ) * t;
}

/** Linear interp. with t = 0.5 */
static inline void lerpHalf( const PointF *a, const PointF *b, PointF *result )
{
    result->x = ( a->x + b->x ) * .5;
    result->y = ( a->y + b->y ) * .5;
}

/** Turns a json array with two children into a point (x, y tuple). */
void jsonGetPoint( cJSON *json, PointF *point )
{
    if ( cJSON_GetArraySize( json ) == 2 )
    {
        point->x = json->child->valuedouble;
        point->y = json->child->next->valuedouble;
    }
}

void jsonSetPoint( cJSON *json, PointF *point )
{
    if ( cJSON_GetArraySize( json ) == 2 )
    {
        json->child->valuedouble       = point->x;
        json->child->next->valuedouble = point->y;
    }
}

void normalizePoints( BPointF *points, int count, int width, int height )
{
    int i;
    for ( i = 0; i < count; i++ )
    {
        points[i].h1.x /= width;
        points[i].p.x  /= width;
        points[i].h2.x /= width;
        points[i].h1.y /= height;
        points[i].p.y  /= height;
        points[i].h2.y /= height;
    }
}

void denormalizePoints( BPointF *points, int count, int width, int height )
{
    int i;
    for ( i = 0; i < count; i++ )
    {
        points[i].h1.x *= width;
        points[i].p.x  *= width;
        points[i].h2.x *= width;
        points[i].h1.y *= height;
        points[i].p.y  *= height;
        points[i].h2.y *= height;
    }
}

/**
 * Turns the array of json elements into an array of Bézier points.
 * \param array cJSON array. values have to be Bézier points: handle 1, point , handl2
 *              ( [ [ [h1x, h1y], [px, py], [h2x, h2y] ], ... ] )
 * \param points pointer to array of points. Will be allocated and filled with the points in \param array
 * \return number of points
 */
int json2BCurves( cJSON *array, BPointF **points, int *tracked )
{
    int count = cJSON_GetArraySize( array );
    cJSON *child = array->child;
    *points = mlt_pool_alloc( count * sizeof( BPointF ) );
    if ( tracked )
    {
        *tracked = 0;
        if ( child && child->type == cJSON_String )
        {
            // keyframe created by tracker ?
            if ( strcmp( child->valuestring, "t" ) == 0 )
                *tracked = 1;
            child = child->next;
        }
    }

    int i = 0;
    do
    {
        if ( child && cJSON_GetArraySize( child ) == 3 )
        {
            jsonGetPoint( child->child , &(*points)[i].h1 );
            jsonGetPoint( child->child->next, &(*points)[i].p );
            jsonGetPoint( child->child->next->next, &(*points)[i].h2 );
            i++;
        }
    } while ( child && ( child = child->next ) );

    if ( i < count )
        *points = mlt_pool_realloc( *points, i * sizeof( BPointF ) );

    return i;
}

/** Determines the point in the middle of the Bézier curve (t = 0.5) defined by \param p1 and \param p2
 * using De Casteljau's algorithm.
 */
static void deCasteljau( BPointF *p1, BPointF *p2, BPointF *mid )
{
    struct PointF ab, bc, cd;

    lerpHalf( &(p1->p), &(p1->h2), &ab );
    lerpHalf( &(p1->h2), &(p2->h1), &bc );
    lerpHalf( &(p2->h1), &(p2->p), &cd );
    lerpHalf( &ab, &bc, &(mid->h1) ); // mid->h1 = abbc
    lerpHalf( &bc, &cd, &(mid->h2) ); // mid->h2 = bccd
    lerpHalf( &(mid->h1), &(mid->h2), &(mid->p) );

    p1->h2 = ab;
    p2->h1 = cd;
}

void curvePoints( BPointF p1, BPointF p2, PointF **points, int *count, int *size )
{
    double errorSqr = SQR( p1.p.x - p2.p.x ) + SQR( p1.p.y - p2.p.y );

    if ( *size + 1 >= *count )
    {
        *size += (int)sqrt( errorSqr / 2 );
        *points = mlt_pool_realloc( *points, *size * sizeof ( struct PointF ) );
    }
    
    (*points)[(*count)++] = p1.p;

    if ( errorSqr <= 2 )
        return;

    BPointF mid;
    deCasteljau( &p1, &p2, &mid );

    curvePoints( p1, mid, points, count, size );

    curvePoints( mid, p2, points, count, size );

    (*points)[*(count)++] = p2.p;
}

int getSplineAt( cJSON *root, mlt_position time, mlt_position in, BPointF **points, int *count, int *isOriginalKeyframe )
{
    *isOriginalKeyframe = 0;
    if ( root == NULL )
        return 0;

    if ( root->type == cJSON_Array )
    {
        /*
         * constant
         */
        *count = json2BCurves( root, points, NULL );
        if ( time == in )
            *isOriginalKeyframe = 1;
    }
    else if ( root->type == cJSON_Object )
    {
        /*
         * keyframes
         */

        int i;
        mlt_position pos1, pos2;

        cJSON *keyframe = root->child;
        cJSON *keyframeOld = keyframe;

        if ( !keyframe )
            return 0;

        while ( atoi( keyframe->string ) < time && keyframe->next )
        {
            keyframeOld = keyframe;
            keyframe = keyframe->next;
        }

        pos1 = atoi( keyframeOld->string );
        pos2 = atoi( keyframe->string );

        if ( pos1 >= pos2 || time >= pos2 )
        {
            // keyframes in wrong order or before/at first / after/at last keyframe
            int tracked = 0;
            *count = json2BCurves( keyframe, points, &tracked );

            if ( !tracked && ( time == pos1 || time == pos2 ) )
                *isOriginalKeyframe = 1;
        }
        else
        {
            /*
             * pos1 < time < pos2
             */

            BPointF *p1, *p2;
            int c1 = json2BCurves( keyframeOld, &p1, NULL );
            int c2 = json2BCurves( keyframe, &p2, NULL );

            // range 0-1
            double position = ( time - pos1 ) / (double)( pos2 - pos1 + 1 );

            *count = MIN( c1, c2 );  // additional points are ignored
            *points = mlt_pool_alloc( *count * sizeof( BPointF ) );
            for ( i = 0; i < *count; i++ )
            {
                lerp( &(p1[i].h1), &(p2[i].h1), &((*points)[i].h1), position );
                lerp( &(p1[i].p), &(p2[i].p), &((*points)[i].p), position );
                lerp( &(p1[i].h2), &(p2[i].h2), &((*points)[i].h2), position );
            }

            mlt_pool_release( p1 );
            mlt_pool_release( p2 );
        }
    }
    else
    {
        return 0;
    }

    return 1;
}

void jsonCreatePoint( cJSON *parent, PointF *point )
{
    cJSON *p = cJSON_CreateArray();
    cJSON *x = cJSON_CreateNumber( point->x );
    cJSON *y = cJSON_CreateNumber( point->y );
    cJSON_AddItemToArray( p, x );
    cJSON_AddItemToArray( p, y );
    cJSON_AddItemToArray( parent, p );
}

void setKeyframeAt( cJSON* root, mlt_position time, BPointF* points, int count, int isOriginalKeyframe, int keyWidth )
{
    if ( !root )
        return;

    if ( root->type == cJSON_Array )
    {
        // TODO: insert conversion to cJSON_Object here
        return;
    }

    if ( !keyWidth )
        keyWidth = strlen( root->child->string );

    int i = 0;
    cJSON *newKeyframe = cJSON_CreateArray();
    if ( !isOriginalKeyframe )
    {
        cJSON *flag = cJSON_CreateString( "t" );
        cJSON_AddItemToArray( newKeyframe, flag );
    }
    for ( ; i < count; ++i )
    {
        cJSON *bpoint = cJSON_CreateArray();
        jsonCreatePoint( bpoint, &points[i].h1 );
        jsonCreatePoint( bpoint, &points[i].p );
        jsonCreatePoint( bpoint, &points[i].h2 );
        cJSON_AddItemToArray( newKeyframe, bpoint );
    }

    char *key = malloc( ( keyWidth + 1 ) * sizeof( char ) );
    for ( i = 0; i < keyWidth; ++i )
        key[i] = '0';
    
    sprintf( key + keyWidth - ( time == 0 ? 0 : (int)( log10( time ) + 1 ) ), "%d", time );
    key[keyWidth] = '\0';

    cJSON *nextKeyframe = root->child;

    while ( atoi ( nextKeyframe->string ) < time && nextKeyframe->next )
    {
        nextKeyframe = nextKeyframe->next;
    }

    if ( atoi( nextKeyframe->string ) == time )
    {
        cJSON_ReplaceItemInObject( root, nextKeyframe->string, newKeyframe );
    }
    else if ( atoi( nextKeyframe->string ) < time )
    {
        cJSON_AddItemToObject( root, key, newKeyframe );
    }
    else
    {
        newKeyframe->string = strdup( key );
        if (nextKeyframe->prev)
        {
            newKeyframe->prev = nextKeyframe->prev;
            nextKeyframe->prev->next = newKeyframe;
        }
        else
        {
            newKeyframe->prev = NULL;
            root->child = newKeyframe;
        }
        newKeyframe->next = nextKeyframe;
        nextKeyframe->prev = newKeyframe;
    }

    free( key );
}


void setStatic( cJSON *root, BPointF *points )
{
    if ( !root || !root->type == cJSON_Array )
        return;

    int count = cJSON_GetArraySize( root );
    cJSON *child = root->child;
    int i = 0;

    do
    {
        jsonSetPoint( child->child , &points[i].h1 );
        jsonSetPoint( child->child->next, &points[i].p );
        jsonSetPoint( child->child->next->next, &points[i].h2 );
        i++;
    } while ( child && ( child = child->next ) );
}
