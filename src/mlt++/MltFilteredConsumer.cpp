/**
 * MltFilteredConsumer.cpp - MLT Wrapper
 * Copyright (C) 2004-2005 Charles Yates
 * Author: Charles Yates <charles.yates@pandora.be>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "MltFilteredConsumer.h"
using namespace Mlt;

FilteredConsumer::FilteredConsumer( Profile& profile, char *id, char *arg ) :
	Consumer( profile, id, arg )
{
	// Create a reference to the first service
	first = new Service( *this );
}

FilteredConsumer::FilteredConsumer( Consumer &consumer ) :
	Consumer( consumer )
{
	// Create a reference to the first service
	first = new Service( *this );
}

FilteredConsumer::~FilteredConsumer( )
{
	// Delete the reference to the first service
	delete first;
}

int FilteredConsumer::connect( Service &service )
{
	// All producers must connect to the first service, hence the use of the virtual here
	return first->connect_producer( service );
}

int FilteredConsumer::attach( Filter &filter )
{
	int error = 0;
	if ( filter.is_valid( ) )
	{
		Service *producer = first->producer( );
		error = filter.connect( *producer );
		if ( error == 0 )
		{
			first->connect_producer( filter );
			delete first;
			first = new Service( filter );
		}
		delete producer;
	}
	else
	{
		error = 1;
	}
	return error;
}

int FilteredConsumer::last( Filter &filter )
{
	int error = 0;
	if ( filter.is_valid( ) )
	{
		Service *producer = this->producer( );
		error = filter.connect( *producer );
		if ( error == 0 )
			connect_producer( filter );
		delete producer;
	}
	else
	{
		error = 1;
	}
	return error;
}

int FilteredConsumer::detach( Filter &filter )
{
	if ( filter.is_valid( ) )
	{
		Service *it = new Service( *first );
		while ( it->is_valid( ) && it->get_service( ) != filter.get_service( ) )
		{
			Service *consumer = it->consumer( );
			delete it;
			it = consumer;
		}
		if ( it->get_service( ) == filter.get_service( ) )
		{
			Service *producer = it->producer( );
			Service *consumer = it->consumer( );
			consumer->connect_producer( *producer );
			Service dummy( NULL );
			it->connect_producer( dummy );
			if ( first->get_service( ) == it->get_service( ) )
			{
				delete first;
				first = new Service( *consumer );
			}
		}
		delete it;
	}
	return 0;
}

