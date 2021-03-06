/**
 * MltFilteredConsumer.h - MLT Wrapper
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

#ifndef _MLTPP_FILTERED_CONSUMER_H_
#define _MLTPP_FILTERED_CONSUMER_H_

#include "config.h"

#include "MltConsumer.h"
#include "MltFilter.h"
#include "MltService.h"

namespace Mlt
{
	class Consumer;
	class Service;
	class Filter;
	class Profile;

	class MLTPP_DECLSPEC FilteredConsumer : public Consumer
	{
		private:
			Service *first;
		public:
			FilteredConsumer( Profile& profile, char *id, char *arg = NULL );
			FilteredConsumer( Consumer &consumer );
			virtual ~FilteredConsumer( );
			int connect( Service &service );
			int attach( Filter &filter );
			int last( Filter &filter );
			int detach( Filter &filter );
	};
}

#endif

