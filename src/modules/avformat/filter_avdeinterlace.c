/*
 * filter_avdeinterlace.c -- deinterlace filter
 * Copyright (C) 2003-2004 Ushodaya Enterprises Limited
 * Author: Charles Yates <charles.yates@pandora.be>
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

#include "filter_avdeinterlace.h"

#include <framework/mlt_frame.h>

#include <string.h>
#include <stdlib.h>

// ffmpeg Header files
#include <ffmpeg/avformat.h>

#ifdef USE_MMX
#include "mmx.h"
#endif

#ifdef USE_MMX
#define DEINT_INPLACE_LINE_LUM \
                    movd_m2r(lum_m4[0],mm0);\
                    movd_m2r(lum_m3[0],mm1);\
                    movd_m2r(lum_m2[0],mm2);\
                    movd_m2r(lum_m1[0],mm3);\
                    movd_m2r(lum[0],mm4);\
                    punpcklbw_r2r(mm7,mm0);\
                    movd_r2m(mm2,lum_m4[0]);\
                    punpcklbw_r2r(mm7,mm1);\
                    punpcklbw_r2r(mm7,mm2);\
                    punpcklbw_r2r(mm7,mm3);\
                    punpcklbw_r2r(mm7,mm4);\
                    paddw_r2r(mm3,mm1);\
                    psllw_i2r(1,mm2);\
                    paddw_r2r(mm4,mm0);\
                    psllw_i2r(2,mm1);\
                    paddw_r2r(mm6,mm2);\
                    paddw_r2r(mm2,mm1);\
                    psubusw_r2r(mm0,mm1);\
                    psrlw_i2r(3,mm1);\
                    packuswb_r2r(mm7,mm1);\
                    movd_r2m(mm1,lum_m2[0]);

#define DEINT_LINE_LUM \
                    movd_m2r(lum_m4[0],mm0);\
                    movd_m2r(lum_m3[0],mm1);\
                    movd_m2r(lum_m2[0],mm2);\
                    movd_m2r(lum_m1[0],mm3);\
                    movd_m2r(lum[0],mm4);\
                    punpcklbw_r2r(mm7,mm0);\
                    punpcklbw_r2r(mm7,mm1);\
                    punpcklbw_r2r(mm7,mm2);\
                    punpcklbw_r2r(mm7,mm3);\
                    punpcklbw_r2r(mm7,mm4);\
                    paddw_r2r(mm3,mm1);\
                    psllw_i2r(1,mm2);\
                    paddw_r2r(mm4,mm0);\
                    psllw_i2r(2,mm1);\
                    paddw_r2r(mm6,mm2);\
                    paddw_r2r(mm2,mm1);\
                    psubusw_r2r(mm0,mm1);\
                    psrlw_i2r(3,mm1);\
                    packuswb_r2r(mm7,mm1);\
                    movd_r2m(mm1,dst[0]);
#endif

/* filter parameters: [-1 4 2 4 -1] // 8 */
static void deinterlace_line(uint8_t *dst, 
			     const uint8_t *lum_m4, const uint8_t *lum_m3, 
			     const uint8_t *lum_m2, const uint8_t *lum_m1, 
			     const uint8_t *lum,
			     int size)
{
#ifndef USE_MMX
    uint8_t *cm = cropTbl + MAX_NEG_CROP;
    int sum;

    for(;size > 0;size--) {
        sum = -lum_m4[0];
        sum += lum_m3[0] << 2;
        sum += lum_m2[0] << 1;
        sum += lum_m1[0] << 2;
        sum += -lum[0];
        dst[0] = cm[(sum + 4) >> 3];
        lum_m4++;
        lum_m3++;
        lum_m2++;
        lum_m1++;
        lum++;
        dst++;
    }
#else

    {
        mmx_t rounder;
        rounder.uw[0]=4;
        rounder.uw[1]=4;
        rounder.uw[2]=4;
        rounder.uw[3]=4;
        pxor_r2r(mm7,mm7);
        movq_m2r(rounder,mm6);
    }
    for (;size > 3; size-=4) {
        DEINT_LINE_LUM
        lum_m4+=4;
        lum_m3+=4;
        lum_m2+=4;
        lum_m1+=4;
        lum+=4;
        dst+=4;
    }
#endif
}
static void deinterlace_line_inplace(uint8_t *lum_m4, uint8_t *lum_m3, uint8_t *lum_m2, uint8_t *lum_m1, uint8_t *lum,
                             int size)
{
#ifndef USE_MMX
    uint8_t *cm = cropTbl + MAX_NEG_CROP;
    int sum;

    for(;size > 0;size--) {
        sum = -lum_m4[0];
        sum += lum_m3[0] << 2;
        sum += lum_m2[0] << 1;
        lum_m4[0]=lum_m2[0];
        sum += lum_m1[0] << 2;
        sum += -lum[0];
        lum_m2[0] = cm[(sum + 4) >> 3];
        lum_m4++;
        lum_m3++;
        lum_m2++;
        lum_m1++;
        lum++;
    }
#else

    {
        mmx_t rounder;
        rounder.uw[0]=4;
        rounder.uw[1]=4;
        rounder.uw[2]=4;
        rounder.uw[3]=4;
        pxor_r2r(mm7,mm7);
        movq_m2r(rounder,mm6);
    }
    for (;size > 3; size-=4) {
        DEINT_INPLACE_LINE_LUM
        lum_m4+=4;
        lum_m3+=4;
        lum_m2+=4;
        lum_m1+=4;
        lum+=4;
    }
#endif
}

/* deinterlacing : 2 temporal taps, 3 spatial taps linear filter. The
   top field is copied as is, but the bottom field is deinterlaced
   against the top field. */
static void deinterlace_bottom_field(uint8_t *dst, int dst_wrap,
                                    const uint8_t *src1, int src_wrap,
                                    int width, int height)
{
    const uint8_t *src_m2, *src_m1, *src_0, *src_p1, *src_p2;
    int y;

    src_m2 = src1;
    src_m1 = src1;
    src_0=&src_m1[src_wrap];
    src_p1=&src_0[src_wrap];
    src_p2=&src_p1[src_wrap];
    for(y=0;y<(height-2);y+=2) {
        memcpy(dst,src_m1,width);
        dst += dst_wrap;
        deinterlace_line(dst,src_m2,src_m1,src_0,src_p1,src_p2,width);
        src_m2 = src_0;
        src_m1 = src_p1;
        src_0 = src_p2;
        src_p1 += 2*src_wrap;
        src_p2 += 2*src_wrap;
        dst += dst_wrap;
    }
    memcpy(dst,src_m1,width);
    dst += dst_wrap;
    /* do last line */
    deinterlace_line(dst,src_m2,src_m1,src_0,src_0,src_0,width);
}

static void deinterlace_bottom_field_inplace(uint8_t *src1, int src_wrap,
					     int width, int height)
{
    uint8_t *src_m1, *src_0, *src_p1, *src_p2;
    int y;
    uint8_t *buf;
    buf = (uint8_t*)av_malloc(width);

    src_m1 = src1;
    memcpy(buf,src_m1,width);
    src_0=&src_m1[src_wrap];
    src_p1=&src_0[src_wrap];
    src_p2=&src_p1[src_wrap];
    for(y=0;y<(height-2);y+=2) {
        deinterlace_line_inplace(buf,src_m1,src_0,src_p1,src_p2,width);
        src_m1 = src_p1;
        src_0 = src_p2;
        src_p1 += 2*src_wrap;
        src_p2 += 2*src_wrap;
    }
    /* do last line */
    deinterlace_line_inplace(buf,src_m1,src_0,src_0,src_0,width);
    av_free(buf);
}


/* deinterlace - if not supported return -1 */
static int mlt_avpicture_deinterlace(AVPicture *dst, const AVPicture *src,
                          int pix_fmt, int width, int height)
{
    int i;

    if (pix_fmt != PIX_FMT_YUV420P &&
        pix_fmt != PIX_FMT_YUV422P &&
        pix_fmt != PIX_FMT_YUV422 &&
        pix_fmt != PIX_FMT_YUV444P &&
	pix_fmt != PIX_FMT_YUV411P)
        return -1;
    if ((width & 3) != 0 || (height & 3) != 0)
        return -1;

	if ( pix_fmt != PIX_FMT_YUV422 )
	{
      for(i=0;i<3;i++) {
          if (i == 1) {
              switch(pix_fmt) {
              case PIX_FMT_YUV420P:
                  width >>= 1;
                  height >>= 1;
                  break;
              case PIX_FMT_YUV422P:
                  width >>= 1;
                  break;
              case PIX_FMT_YUV411P:
                  width >>= 2;
                  break;
              default:
                  break;
              }
          }
          if (src == dst) {
              deinterlace_bottom_field_inplace(dst->data[i], dst->linesize[i],
                                   width, height);
          } else {
              deinterlace_bottom_field(dst->data[i],dst->linesize[i],
                                          src->data[i], src->linesize[i],
                                          width, height);
          }
	  }
    }
	else {
      if (src == dst) {
          deinterlace_bottom_field_inplace(dst->data[0], dst->linesize[0],
                               width, height);
      } else {
          deinterlace_bottom_field(dst->data[0],dst->linesize[0],
                                      src->data[0], src->linesize[0],
                                      width, height);
      }
	}

#ifdef USE_MMX
    emms();
#endif
    return 0;
}

/** Do it :-).
*/

static int filter_get_image( mlt_frame this, uint8_t **image, mlt_image_format *format, int *width, int *height, int writable )
{
	int error = 0;
	
	// Check that we want progressive and we aren't already progressive
	if ( *format == mlt_image_yuv422 &&
		 !mlt_properties_get_int( mlt_frame_properties( this ), "progressive" ) &&
		 mlt_properties_get_int( mlt_frame_properties( this ), "consumer_deinterlace" ) )
	{
		// Create a picture
		AVPicture *output = mlt_pool_alloc( sizeof( AVPicture ) );

		// Get the input image
		error = mlt_frame_get_image( this, image, format, width, height, 1 );

		// Fill the picture
		if ( *format == mlt_image_yuv422 )
		{
			avpicture_fill( output, *image, PIX_FMT_YUV422, *width, *height );
			mlt_avpicture_deinterlace( output, output, PIX_FMT_YUV422, *width, *height );
		}

		// Free the picture
		mlt_pool_release( output );

		// Make sure that others know the frame is deinterlaced
		mlt_properties_set_int( mlt_frame_properties( this ), "progressive", 1 );
	}
	else
	{
		// Get the input image
		error = mlt_frame_get_image( this, image, format, width, height, writable );
	}

	return error;
}

/** Deinterlace filter processing - this should be lazy evaluation here...
*/

static mlt_frame deinterlace_process( mlt_filter this, mlt_frame frame )
{
	// Push the get_image method on to the stack
	mlt_frame_push_get_image( frame, filter_get_image );
	
	return frame;
}

/** Constructor for the filter.
*/

mlt_filter filter_avdeinterlace_init( void *arg )
{
	mlt_filter this = mlt_filter_new( );
	if ( this != NULL )
		this->process = deinterlace_process;
	return this;
}
