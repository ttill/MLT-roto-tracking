/*
 * filter_sox.c -- apply any number of SOX effects using libst
 * Copyright (C) 2003-2004 Ushodaya Enterprises Limited
 * Author: Dan Dennedy <dan@dennedy.org>
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

#include "filter_sox.h"

#include <framework/mlt_frame.h>
#include <framework/mlt_tokeniser.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sox/st.h>

#define BUFFER_LEN 8192
#define AMPLITUDE_NORM 0.2511886431509580 /* -12dBFS */
#define AMPLITUDE_MIN 0.00001

/** Compute the mean of a set of doubles skipping unset values flagged as -1
*/
static inline double mean( double *buf, int count )
{
	double mean = 0;
	int i;
	int j = 0;
	
	for ( i = 0; i < count; i++ )
	{
		if ( buf[ i ] != -1.0 )
		{
			mean += buf[ i ];
			j ++;
		}
	}
	if ( j > 0 )
		mean /= j;
	
	return mean;
}

/** Create an effect state instance for a channels
*/
static int create_effect( mlt_filter this, char *value, int count, int channel, int frequency )
{
	mlt_tokeniser tokeniser = mlt_tokeniser_init();
	eff_t eff = mlt_pool_alloc( sizeof( struct st_effect ) );
	char id[ 256 ];
	int error = 1;

	// Tokenise the effect specification
	mlt_tokeniser_parse_new( tokeniser, value, " " );

	// Locate the effect
	int opt_count = st_geteffect_opt( eff, tokeniser->count, tokeniser->tokens );
	
	// If valid effect
	if ( opt_count != ST_EOF )
	{
		// Supply the effect parameters
		if ( ( * eff->h->getopts )( eff, opt_count, &tokeniser->tokens[ tokeniser->count - opt_count ] ) == ST_SUCCESS )
		{
			// Set the sox signal parameters
			eff->ininfo.rate = frequency;
			eff->outinfo.rate = frequency;
			eff->ininfo.channels = 1;
			eff->outinfo.channels = 1;
			
			// Start the effect
			if ( ( * eff->h->start )( eff ) == ST_SUCCESS )
			{
				// Construct id
				sprintf( id, "_effect_%d_%d", count, channel );

				// Save the effect state
				mlt_properties_set_data( mlt_filter_properties( this ), id, eff, 0, mlt_pool_release, NULL );
				error = 0;
			}
		}
	}
	// Some error occurred so delete the temp effect state
	if ( error == 1 )
		mlt_pool_release( eff );
	
	mlt_tokeniser_close( tokeniser );
	
	return error;
}

/** Get the audio.
*/

static int filter_get_audio( mlt_frame frame, int16_t **buffer, mlt_audio_format *format, int *frequency, int *channels, int *samples )
{
	// Get the properties of the frame
	mlt_properties properties = mlt_frame_properties( frame );

	// Get the filter service
	mlt_filter filter = mlt_frame_pop_audio( frame );

	// Get the filter properties
	mlt_properties filter_properties = mlt_filter_properties( filter );

	// Get the properties
	st_sample_t *input_buffer = mlt_properties_get_data( filter_properties, "input_buffer", NULL );
	st_sample_t *output_buffer = mlt_properties_get_data( filter_properties, "output_buffer", NULL );
	int channels_avail = *channels;
	int i; // channel
	int count = mlt_properties_get_int( filter_properties, "effect_count" );

	// Restore the original get_audio
	frame->get_audio = mlt_frame_pop_audio( frame );

	// Get the producer's audio
	mlt_frame_get_audio( frame, buffer, format, frequency, &channels_avail, samples );

	// Duplicate channels as necessary
	if ( channels_avail < *channels )
	{
		int size = *channels * *samples * sizeof( int16_t );
		int16_t *new_buffer = mlt_pool_alloc( size );
		int j, k = 0;
		
		// Duplicate the existing channels
		for ( i = 0; i < *samples; i++ )
		{
			for ( j = 0; j < *channels; j++ )
			{
				new_buffer[ ( i * *channels ) + j ] = (*buffer)[ ( i * channels_avail ) + k ];
				k = ( k + 1 ) % channels_avail;
			}
		}
		
		// Update the audio buffer now - destroys the old
		mlt_properties_set_data( properties, "audio", new_buffer, size, ( mlt_destructor )mlt_pool_release, NULL );
		
		*buffer = new_buffer;
	}
	else if ( channels_avail == 6 && *channels == 2 )
	{
		// Nasty hack for ac3 5.1 audio - may be a cause of failure?
		int size = *channels * *samples * sizeof( int16_t );
		int16_t *new_buffer = mlt_pool_alloc( size );
		
		// Drop all but the first *channels
		for ( i = 0; i < *samples; i++ )
		{
			new_buffer[ ( i * *channels ) + 0 ] = (*buffer)[ ( i * channels_avail ) + 2 ];
			new_buffer[ ( i * *channels ) + 1 ] = (*buffer)[ ( i * channels_avail ) + 3 ];
		}

		// Update the audio buffer now - destroys the old
		mlt_properties_set_data( properties, "audio", new_buffer, size, ( mlt_destructor )mlt_pool_release, NULL );
		
		*buffer = new_buffer;
	}

	// Even though some effects are multi-channel aware, it is not reliable
	// We must maintain a separate effect state for each channel
	for ( i = 0; i < *channels; i++ )
	{
		char id[ 256 ];
		sprintf( id, "_effect_0_%d", i );
		
		// Get an existing effect state
		eff_t e = mlt_properties_get_data( filter_properties, id, NULL );
		
		// Validate the existing effect state
		if ( e != NULL && ( e->ininfo.rate != *frequency || 
							e->outinfo.rate != *frequency ) )
			e = NULL;
		
		// (Re)Create the effect state
		if ( e == NULL )
		{
			int j = 0;
			
			// Reset the count
			count = 0;
	
			// Loop over all properties
			for ( j = 0; j < mlt_properties_count( filter_properties ); j ++ )
			{
				// Get the name of this property
				char *name = mlt_properties_get_name( filter_properties, j );
	
				// If the name does not contain a . and matches effect
				if ( !strncmp( name, "effect", 6 ) )
				{
					// Get the effect specification
					char *value = mlt_properties_get( filter_properties, name );
	
					// Create an instance
					if ( create_effect( filter, value, count, i, *frequency ) == 0 )
						count ++;
				}
			}
			
			// Save the number of filters
			mlt_properties_set_int( filter_properties, "effect_count", count );
			
		}
		if ( *samples > 0 && count > 0 )
		{
			st_sample_t *p = input_buffer;
			st_sample_t *end = p + *samples;
			int16_t *q = *buffer + i;
			st_size_t isamp = *samples;
			st_size_t osamp = *samples;
			double rms = 0;
			int j;
			char *normalise = mlt_properties_get( filter_properties, "normalise" );
			double normalised_gain = 1.0;
			
			// Convert to sox encoding
			while( p != end )
			{
				*p = ST_SIGNED_WORD_TO_SAMPLE( *q );
				
				// Compute rms amplitude while we are accessing each sample
				rms += ( double )*p * ( double )*p;
				
				p ++;
				q += *channels;
			}
			
			// Compute final rms amplitude
			rms = sqrt( rms / *samples / ST_SSIZE_MIN / ST_SSIZE_MIN );
			
			if ( normalise )
			{
				int window = mlt_properties_get_int( filter_properties, "window" );
				double *smooth_buffer = mlt_properties_get_data( filter_properties, "smooth_buffer", NULL );
				double max_gain = mlt_properties_get_double( filter_properties, "max_gain" );
				
				// Default the maximum gain factor to 20dBFS
				if ( max_gain == 0 )
					max_gain = 10.0;
				
				// The smoothing buffer prevents radical shifts in the gain level
				if ( window > 0 && smooth_buffer != NULL )
				{
					int smooth_index = mlt_properties_get_int( filter_properties, "_smooth_index" );
					smooth_buffer[ smooth_index ] = rms;
					
					// Ignore very small values that adversely affect the mean
					if ( rms > AMPLITUDE_MIN )
						mlt_properties_set_int( filter_properties, "_smooth_index", ( smooth_index + 1 ) % window );
					
					// Smoothing is really just a mean over the past N values
					normalised_gain = AMPLITUDE_NORM / mean( smooth_buffer, window );
				}
				else if ( rms > 0 )
				{
					// Determine gain to apply as current amplitude
					normalised_gain = AMPLITUDE_NORM / rms;
				}
					
				//printf("filter_sox: rms %.3f gain %.3f\n", rms, normalised_gain );
				
				// Govern the maximum gain
				if ( normalised_gain > max_gain )
					normalised_gain = max_gain;
			}
			
			// For each effect
			for ( j = 0; j < count; j++ )
			{
				sprintf( id, "_effect_%d_%d", j, i );
				e = mlt_properties_get_data( filter_properties, id, NULL );
				
				// We better have this guy
				if ( e != NULL )
				{
					float saved_gain = 1.0;
					
					// XXX: hack to apply the normalised gain level to the vol effect
					if ( normalise && strcmp( e->name, "vol" ) == 0 )
					{
						float *f = ( float * )( e->priv );
						saved_gain = *f;
						*f = saved_gain * normalised_gain;
					}
					
					// Apply the effect
					if ( ( * e->h->flow )( e, input_buffer, output_buffer, &isamp, &osamp ) == ST_SUCCESS )
					{
						// Swap input and output buffer pointers for subsequent effects
						p = input_buffer;
						input_buffer = output_buffer;
						output_buffer = p;
					}
					
					// XXX: hack to restore the original vol gain to prevent accumulation
					if ( normalise && strcmp( e->name, "vol" ) == 0 )
					{
						float *f = ( float * )( e->priv );
						*f = saved_gain;
					}
				}
			}
			
			// Convert back to signed 16bit
			p = input_buffer;
			q = *buffer + i;
			end = p + *samples;
			while ( p != end )
			{
				*q = ST_SAMPLE_TO_SIGNED_WORD( *p ++ );
				q += *channels;
			}
		}
	}

	return 0;
}

/** Filter processing.
*/

static mlt_frame filter_process( mlt_filter this, mlt_frame frame )
{
	if ( frame->get_audio != NULL )
	{
		// Add the filter to the frame
		mlt_frame_push_audio( frame, frame->get_audio );
		mlt_frame_push_audio( frame, this );
		frame->get_audio = filter_get_audio;
		
		// Parse the window property and allocate smoothing buffer if needed
		mlt_properties properties = mlt_filter_properties( this );
		int window = mlt_properties_get_int( properties, "window" );
		if ( mlt_properties_get( properties, "smooth_buffer" ) == NULL && window > 1 )
		{
			// Create a smoothing buffer for the calculated "max power" of frame of audio used in normalisation
			double *smooth_buffer = (double*) calloc( window, sizeof( double ) );
			int i;
			for ( i = 0; i < window; i++ )
				smooth_buffer[ i ] = -1.0;
			mlt_properties_set_data( properties, "smooth_buffer", smooth_buffer, 0, free, NULL );
		}
	}

	return frame;
}

/** Constructor for the filter.
*/

mlt_filter filter_sox_init( char *arg )
{
	mlt_filter this = mlt_filter_new( );
	if ( this != NULL )
	{
		void *input_buffer = mlt_pool_alloc( BUFFER_LEN );
		void *output_buffer = mlt_pool_alloc( BUFFER_LEN );
		mlt_properties properties = mlt_filter_properties( this );
		
		this->process = filter_process;
		
		if ( arg != NULL )
			mlt_properties_set( properties, "effect", arg );
		mlt_properties_set_data( properties, "input_buffer", input_buffer, BUFFER_LEN, mlt_pool_release, NULL );
		mlt_properties_set_data( properties, "output_buffer", output_buffer, BUFFER_LEN, mlt_pool_release, NULL );
		mlt_properties_set_int( properties, "window", 75 );
	}
	return this;
}

// What to do when a libst internal failure occurs
void cleanup(void){}

// Is there a build problem with my sox-devel package?
#ifndef gsm_create
void gsm_create(void){}
#endif
#ifndef gsm_decode
void gsm_decode(void){}
#endif
#ifndef gdm_encode
void gsm_encode(void){}
#endif
#ifndef gsm_destroy
void gsm_destroy(void){}
#endif
#ifndef gsm_option
void gsm_option(void){}
#endif