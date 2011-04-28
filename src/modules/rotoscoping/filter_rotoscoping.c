/*
 * rotoscoping.c -- keyframable vector based rotoscoping
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

#include <framework/mlt_filter.h>
#include <framework/mlt_frame.h>

#include "deformation.h"
#include "spline_handling.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <opencv/cv.h>

#ifndef MAX
#define MAX( x, y ) ( ( x ) > ( y ) ? ( x ) : ( y ) )
#endif
#ifndef MIN
#define MIN( x, y ) ( ( x ) < ( y ) ? ( x ) : ( y ) )
#endif


#define DEBUG_TRACK

typedef struct ocvData
{
    IplImage *last;
    int points_count;
    CvPoint2D32f *points;
} ocvData;

enum MODES { MODE_ALPHA, MODE_LUMA };
const char *MODESTR[2] = { "alpha", "luma" };

enum ALPHAOPERATIONS { ALPHA_CLEAR, ALPHA_MAX, ALPHA_MIN, ALPHA_ADD, ALPHA_SUB };
const char *ALPHAOPERATIONSTR[5] = { "clear", "max", "min", "add", "sub" };


void freeOcvData( ocvData *data )
{
    if ( data->last )
        cvReleaseImage( &data->last );
    free( data->points );
    free( data );
}

/** Returns the index of \param string in \param stringList.
 * Useful for assigning string parameters to enums. */
int stringValue( const char *string, const char **stringList, int max )
{
    int i;
    for ( i = 0; i < max; i++ )
        if ( strcmp( stringList[i], string ) == 0 )
            return i;
    return 0;
}

/** Sets "spline_is_dirty" to 1 if property "spline" was changed.
 * We then know when to parse the json stored in "spline" */
static void rotoPropertyChanged( mlt_service owner, mlt_filter filter, char *name )
{
    if ( !strcmp( name, "spline" ) )
        mlt_properties_set_int( MLT_FILTER_PROPERTIES( filter ), "_spline_is_dirty", 1 );
}

/** Helper for using qsort with an array of integers. */
int ncompare( const void *a, const void *b )
{
    return *(const int*)a - *(const int*)b;
}

/** Blurs \param src horizontally. \See funtion blur. */
void blurHorizontal( uint8_t *src, uint8_t *dst, int width, int height, int radius)
{
    int x, y, kx, yOff, total, amount, amountInit;
    amountInit = radius * 2 + 1;
    for (y = 0; y < height; ++y)
    {
        total = 0;
        yOff = y * width;
        // Process entire window for first pixel
        int size = MIN(radius + 1, width);
        for ( kx = 0; kx < size; ++kx )
            total += src[yOff + kx];
        dst[yOff] = total / ( radius + 1 );
        // Subsequent pixels just update window total
        for ( x = 1; x < width; ++x )
        {
            amount = amountInit;
            // Subtract pixel leaving window
            if ( x - radius - 1 >= 0 )
                total -= src[yOff + x - radius - 1];
            else
                amount -= radius - x;
            // Add pixel entering window
            if ( x + radius < width )
                total += src[yOff + x + radius];
            else
                amount -= radius - width + x;
            dst[yOff + x] = total / amount;
        }
    }
}

/** Blurs \param src vertically. \See funtion blur. */
void blurVertical( uint8_t *src, uint8_t *dst, int width, int height, int radius)
{
    int x, y, ky, total, amount, amountInit;
    amountInit = radius * 2 + 1;
    for (x = 0; x < width; ++x)
    {
        total = 0;
        int size = MIN(radius + 1, height);
        for ( ky = 0; ky < size; ++ky )
            total += src[x + ky * width];
        dst[x] = total / ( radius + 1 );
        for ( y = 1; y < height; ++y )
        {
            amount = amountInit;
            if ( y - radius - 1 >= 0 )
                total -= src[( y - radius - 1 ) * width + x];
            else
                amount -= radius - y;
            if ( y + radius < height )
                total += src[( y + radius ) * width + x];
            else
                amount -= radius - height + y;
            dst[y * width + x] = total / amount;
        }
    }
}

/**
 * Blurs the \param map using a simple "average" blur.
 * \param map Will be blured; 1bpp
 * \param width x dimension of channel stored in \param map
 * \param height y dimension of channel stored in \param map
 * \param radius blur radius
 * \param passes blur passes
 */
void blur( uint8_t *map, int width, int height, int radius, int passes )
{
    uint8_t *src = mlt_pool_alloc( width * height );
    uint8_t *tmp = mlt_pool_alloc( width * height );

    int i;
    for ( i = 0; i < passes; ++i )
    {
        memcpy( src, map, width * height );
        blurHorizontal( src, tmp, width, height, radius );
        blurVertical( tmp, map, width, height, radius );
    }

    mlt_pool_release(src);
    mlt_pool_release(tmp);
}

/**
 * Determines which points are located in the polygon and sets their value in \param map to \param value
 * \param vertices points defining the polygon
 * \param count number of vertices
 * \param with x range
 * \param height y range
 * \param value value identifying points in the polygon
 * \param map array of integers of the dimension width * height.
 *            The map entries belonging to the points in the polygon will be set to \param set * 255 the others to !set * 255.
 */
void fillMap( PointF *vertices, int count, int width, int height, int invert, uint8_t *map )
{
    int nodes, nodeX[1024], pixelY, i, j, value;

    value = !invert * 255;
    memset( map, invert * 255, width * height );

    // Loop through the rows of the image
    for ( pixelY = 0; pixelY < height; pixelY++ )
    {
        /*
         * Build a list of nodes.
         * nodes are located at the borders of the polygon
         * and therefore indicate a move from in to out or vice versa
         */
        nodes = 0;
        for ( i = 0, j = count - 1; i < count; j = i++ )
            if ( (vertices[i].y > (double)pixelY) != (vertices[j].y > (double)pixelY) )
                nodeX[nodes++] = (int)(vertices[i].x + (pixelY - vertices[i].y) / (vertices[j].y - vertices[i].y) * (vertices[j].x - vertices[i].x) );

        qsort( nodeX, nodes, sizeof( int ), ncompare );

        // Set map values for points between the node pairs to 1
        for ( i = 0; i < nodes; i += 2 )
        {
            if ( nodeX[i] >= width )
                break;

            if ( nodeX[i+1] > 0 )
            {
                nodeX[i] = MAX( 0, nodeX[i] );
                nodeX[i+1] = MIN( nodeX[i+1], width );
                memset( map + width * pixelY + nodeX[i], value, nodeX[i+1] - nodeX[i] );
            }
        }
    }
}

/** Do it :-).
*/
static int filter_get_image( mlt_frame frame, uint8_t **image, mlt_image_format *format, int *width, int *height, int writable )
{
    mlt_filter filter = mlt_frame_pop_service( frame );
    mlt_properties filter_properties = MLT_FILTER_PROPERTIES( filter );
    mlt_properties frame_properties = mlt_frame_pop_service( frame );

    int mode = mlt_properties_get_int( frame_properties, "mode" );
    int doTrack = mlt_properties_get_int( filter_properties, "track" );

    // Get the image
    if ( doTrack )
        *format = mlt_image_rgb24a;
    int error = mlt_frame_get_image( frame, image, format, width, height, writable );

    // Only process if we have no error and a valid colour space
    if ( !error )
    {
        BPointF *bpoints;
        int bcount, length, count, size, i, j;

        cJSON *root = mlt_properties_get_data( filter_properties, "_spline_parsed", NULL );
        if ( !splineAt( root, mlt_frame_get_position( frame ), &bpoints, &bcount ) )
            return error;

        length = *width * *height;

        if ( doTrack )
        {
            mlt_service_lock( MLT_FILTER_SERVICE( filter ) );

            CvSize cSize = cvSize( *width, *height );

            IplImage *cOrig = cvCreateImage( cSize, IPL_DEPTH_8U, 4 );
            IplImage *cImg = cvCreateImage( cSize, IPL_DEPTH_8U, 1 );

            memcpy( cOrig->imageData, *image, length * 4 );
            cvCvtColor( cOrig, cImg, CV_RGBA2GRAY );

            ocvData *data = mlt_properties_get_data( filter_properties, "_roto_tracking", NULL );
            if ( !data || !data->points_count )
            {
                if (data)
                    freeOcvData( data );
                data = calloc( 1, sizeof( ocvData ) );
                data->points_count = 1000;
                data->points = malloc( data->points_count * sizeof( CvPoint2D32f ) );
                IplImage *cEig = cvCreateImage( cSize, IPL_DEPTH_32F, 1 );
                IplImage *cTmp = cvCreateImage( cSize, IPL_DEPTH_32F, 1 );
                cvGoodFeaturesToTrack( cImg, cEig, cTmp, data->points, &data->points_count, 0.01, 5, NULL, 3, 0, 0.04 );
                data->last = cImg;
                cvReleaseImage( &cEig );
                cvReleaseImage( &cTmp );

                mlt_properties_set_data( filter_properties, "_roto_tracking", data, 0, (mlt_destructor)freeOcvData, NULL );
            }
            else
            {
                CvPoint2D32f *pNew = calloc( data->points_count, sizeof( CvPoint2D32f ) );
                CvPoint2D32f *pOld = malloc( data->points_count * sizeof( CvPoint2D32f ) );
                char status[data->points_count];
                cvCalcOpticalFlowPyrLK( data->last, cImg, NULL, NULL, data->points, pNew, data->points_count, cvSize(20, 20), 5, status, NULL,
                                        cvTermCriteria(CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, 20, 0.01), 0 );

                cvReleaseImage( &data->last );
                data->last = cImg;

                int points_count = 0;
                for (i = 0; i < data->points_count; ++i )
                {
                    if ( status[i] )
                    {
                        pOld[points_count] = data->points[i];
                        pNew[points_count++] = pNew[i];
                    }
                }

                pOld = realloc( pOld, points_count * sizeof ( CvPoint2D32f ) );
                pNew = realloc( pNew, points_count * sizeof ( CvPoint2D32f ) );

#ifdef DEBUG_TRACK
                for ( i = 0; i < points_count; ++i )
                {
                    cvLine( cOrig, cvPoint( (int)(pOld[i].x + .5), (int)(pOld[i].y + .5) ),
                                   cvPoint( (int)(pNew[i].x + .5), (int)(pNew[i].y + .5) ),
                            cvScalar( 255, 0, 0, 0 ), 1, 8, 0 );
                }
#endif

                if ( mlt_properties_get_int( filter_properties, "deform" ) )
                {
                    cJSON *root = mlt_properties_get_data( filter_properties, "_spline_parsed", NULL );
                    if ( !root || !root->type == cJSON_Array )
                        ;
                    PointF *p = malloc( points_count * sizeof( PointF ) );
                    PointF *q = malloc( points_count * sizeof( PointF ) );
                    for ( i = 0; i < points_count; ++i )
                    {
                        vVecD( pOld[i].x / *width, pOld[i].y / *height, &p[i] );
                        vVecD( pNew[i].x / *width, pNew[i].y / *height, &q[i] );
                    }

                    PointF r;
                    BPointF r1;
                    double a = 2;
                    for ( i = 0; i < bcount; ++i )
                    {
                        deform( p, q, points_count, vVec( bpoints[i].h1.x, bpoints[i].h1.y ), a, &r );
                        r1.h1.x = r.x;
                        r1.h1.y = r.y;
                        deform( p, q, points_count, vVec( bpoints[i].p.x, bpoints[i].p.y ), a, &r );
                        r1.p.x = r.x;
                        r1.p.y = r.y;
                        deform( p, q, points_count, vVec( bpoints[i].h2.x, bpoints[i].h2.y ), a, &r );
                        r1.h2.x = r.x;
                        r1.h2.y = r.y;
                        bpoints[i] = r1;
                    }
                    free( p );
                    free( q );

                    setStatic( root, bpoints );

                }

                free( data->points );
                free( pOld );
                data->points = pNew;
                data->points_count = points_count;

                memcpy( *image, cOrig->imageData, length * 4);
            }

            cvReleaseImage( &cOrig );

            mlt_service_unlock( MLT_FILTER_SERVICE( filter ) );
        }
        else
        {
            mlt_properties_set_data( filter_properties, "_roto_tracking", NULL, 0, NULL, NULL );
        }

        denormalizePoints( bpoints, bcount, *width, *height );

        count = 0;
        size = 1;
        struct PointF *points;
        points = mlt_pool_alloc( size * sizeof( struct PointF ) );
        for ( i = 0; i < bcount; i++ )
        {
            j = (i + 1) % bcount;
            curvePoints( bpoints[i], bpoints[j], &points, &count, &size );
        }

        if ( count )
        {
            uint8_t *map = mlt_pool_alloc( length );
            int invert = mlt_properties_get_int( filter_properties, "invert" );
            fillMap( points, count, *width, *height, invert, map );

            int feather = mlt_properties_get_int( filter_properties, "feather" );
            if ( feather )
                blur( map, *width, *height, feather, mlt_properties_get_int( filter_properties, "feather_passes" ) );

            int bpp;
            size = mlt_image_format_size( *format, *width, *height, &bpp );
            uint8_t *p = *image;
            uint8_t *q = *image + size;

            i = 0;
            uint8_t *alpha;

            switch ( mode )
            {
            case MODE_LUMA:
                switch ( *format )
                {
                    case mlt_image_rgb24:
                    case mlt_image_rgb24a:
                    case mlt_image_opengl:
                        while ( p != q )
                        {
                            p[0] = p[1] = p[2] = map[i++];
                            p += bpp;
                        }
                        break;
                    case mlt_image_yuv422:
                        while ( p != q )
                        {
                            p[0] = map[i++];
                            p[1] = 128;
                            p += 2;
                        }
                        break;
                    case mlt_image_yuv420p:
                        memcpy( p, map, length );
                        memset( p + length, 128, length / 2 );
                        break;
                    default:
                        break;
                }
                break;
            case MODE_ALPHA:
                switch ( *format )
                {
                case mlt_image_rgb24a:
                case mlt_image_opengl:
                    switch ( mlt_properties_get_int( frame_properties, "alpha_operation" ) )
                    {
                    case ALPHA_CLEAR:
                        while ( p != q )
                        {
                            p[3] = map[i++];
                            p += 4;
                        }
                        break;
                    case ALPHA_MAX:
                        while ( p != q )
                        {
                            p[3] = MAX( p[3], map[i] );
                            p += 4;
                            i++;
                        }
                        break;
                    case ALPHA_MIN:
                        while ( p != q )
                        {
                            p[3] = MIN( p[3], map[i] );
                            p += 4;
                            i++;
                        }
                        break;
                    case ALPHA_ADD:
                        while ( p != q )
                        {
                            p[3] = MIN( p[3] + map[i], 255 );
                            p += 4;
                            i++;
                        }
                        break;
                    case ALPHA_SUB:
                        while ( p != q )
                        {
                            p[3] = MAX( p[3] - map[i], 0 );
                            p += 4;
                            i++;
                        }
                        break;
                    }
                    break;
                default:
                    alpha = mlt_frame_get_alpha_mask( frame );
                    switch ( mlt_properties_get_int( frame_properties, "alpha_operation" ) )
                    {
                    case ALPHA_CLEAR:
                        memcpy( alpha, map, length );
                        break;
                    case ALPHA_MAX:
                        for ( ; i < length; i++, alpha++ )
                            *alpha = MAX( map[i], *alpha );
                        break;
                    case ALPHA_MIN:
                        for ( ; i < length; i++, alpha++ )
                            *alpha = MIN( map[i], *alpha );
                        break;
                    case ALPHA_ADD:
                        for ( ; i < length; i++, alpha++ )
                            *alpha = MIN( *alpha + map[i], 255 );
                        break;
                    case ALPHA_SUB:
                        for ( ; i < length; i++, alpha++ )
                            *alpha = MAX( *alpha - map[i], 0 );
                        break;
                    }
                    break;
                }
                break;
            }

            mlt_pool_release( map );
        }

        mlt_pool_release( points );
    }

    return error;
}

/** Filter processing.
*/
static mlt_frame filter_process( mlt_filter filter, mlt_frame frame )
{
    mlt_properties properties = MLT_FILTER_PROPERTIES( filter );
    int splineIsDirty = mlt_properties_get_int( properties, "_spline_is_dirty" );
    char *modeStr = mlt_properties_get( properties, "mode" );
    cJSON *root = mlt_properties_get_data( properties, "_spline_parsed", NULL );

    if ( splineIsDirty || root == NULL )
    {
        // we need to (re-)parse
        char *spline = mlt_properties_get( properties, "spline" );
        root = cJSON_Parse( spline );
        mlt_properties_set_data( properties, "_spline_parsed", root, 0, (mlt_destructor)cJSON_Delete, NULL );
        mlt_properties_set_int( properties, "_spline_is_dirty", 0 );
    }

    mlt_properties unique = mlt_frame_unique_properties( frame, MLT_FILTER_SERVICE( filter ) );
    mlt_properties_set_int( unique, "mode", stringValue( modeStr, MODESTR, 2 ) );
    mlt_properties_set_int( unique, "alpha_operation", stringValue( mlt_properties_get( properties, "alpha_operation" ), ALPHAOPERATIONSTR, 5 ) );
    mlt_frame_push_service( frame, unique );
    mlt_frame_push_service( frame, filter );
    mlt_frame_push_get_image( frame, filter_get_image );

    return frame;
}

/** Constructor for the filter.
*/
mlt_filter filter_rotoscoping_init( mlt_profile profile, mlt_service_type type, const char *id, char *arg )
{
        mlt_filter filter = mlt_filter_new( );
        if ( filter )
        {
                filter->process = filter_process;
                mlt_properties properties = MLT_FILTER_PROPERTIES( filter );
                mlt_properties_set( properties, "mode", "alpha" );
                mlt_properties_set( properties, "alpha_operation", "clear" );
                mlt_properties_set_int( properties, "invert", 0 );
                mlt_properties_set_int( properties, "feather", 0 );
                mlt_properties_set_int( properties, "feather_passes", 1 );
                if ( arg )
                    mlt_properties_set( properties, "spline", arg );

                mlt_events_listen( properties, filter, "property-changed", (mlt_listener)rotoPropertyChanged );
        }
        return filter;
}
