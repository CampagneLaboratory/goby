/**
 * Copyright (C) 2010 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#ifndef GOBY_COMMON_H
#define GOBY_COMMON_H

#include <fcntl.h>

namespace goby {

#if defined(_MSC_VER)
  #ifdef LIBGOBY_EXPORTS
    #define LIBGOBY_EXPORT __declspec(dllexport)
    #define LIBGOBY_EXPIMP_TEMPLATE
  #else
    #define LIBGOBY_EXPORT __declspec(dllimport)
    #define LIBGOBY_EXPIMP_TEMPLATE extern
  #endif
#else
  #define LIBGOBY_EXPORT
#endif

#ifndef O_BINARY
#ifdef _O_BINARY
#define O_BINARY _O_BINARY
#else
#define O_BINARY 0     // If this isn't defined, the platform doesn't need it.
#endif
#endif
}

#endif // GOBY_COMMON_H
