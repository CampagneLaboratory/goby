//
// Copyright (C) 2009-2012 Institute for Computational Biomedicine,
//                         Weill Medical College of Cornell University
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

/**
 * C API / functions to directly support GSNAP alignment / parsing.
 */

#ifndef C_GSNAP_H_
#define C_GSNAP_H_

#include "goby/C_CompactHelpers.h"

#ifdef __cplusplus
extern "C" {
#endif
    void gobyGsnap_startAlignment(CAlignmentsWriterHelper *writerHelper);
    void gobyGsnap_parse(CAlignmentsWriterHelper *writerHelper, char *alignment);
    void gobyGsnap_test_registerTargets(CAlignmentsWriterHelper *writerHelper, char *targets);
    void gobyGsnap_destoryAlignment(CAlignmentsWriterHelper *writerHelper);
#ifdef __cplusplus
}
#endif

#endif /* C_GSNAP_H_ */
