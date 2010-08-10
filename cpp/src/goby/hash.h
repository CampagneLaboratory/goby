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

#ifndef GOBY_HASH_H
#define GOBY_HASH_H

#if HAVE_CONFIG_H
#include <config.h>
#endif

// modern microsoft compilers define hash_map
#ifdef _MSC_VER
#include <hash_map>
#define LIBGOBY_HASH_MAP ::std::hash_map

// if there is support for C++ Technical Report 1
#elif defined(HAVE_TR1_UNORDERED_MAP)
#include <tr1/unordered_map>
#define LIBGOBY_HASH_MAP ::std::tr1::unordered_map

// otherwise try for the hash map extension
#elif defined(HAVE_EXT_HASH_MAP)
#include <ext/hash_map>
namespace std { using namespace __gnu_cxx; }
#define LIBGOBY_HASH_MAP ::std::hash_map

// otherwise just use plain old map
#else
#include <map>
#define LIBGOBY_HASH_MAP ::std::map
#endif

// modern microsoft compilers define hash_set
#ifdef _MSC_VER
#include <hash_set>
#define LIBGOBY_HASH_SET ::std::hash_set

// if there is support for C++ Technical Report 1
#elif defined(HAVE_TR1_UNORDERED_SET)
#include <tr1/unordered_set>
#define LIBGOBY_HASH_SET ::std::tr1::unordered_set

// otherwise try for the hash set extension
#elif defined(HAVE_EXT_HASH_SET)
#include <ext/hash_map>
namespace std { using namespace __gnu_cxx; }
#define LIBGOBY_HASH_SET ::std::hash_set

// otherwise just use plain old set
#else
#include <map>
#define LIBGOBY_HASH_SET ::std::set
#endif

#endif // GOBY_HASH_H
