/**
 * C API / functions to directly support GSNAP alignment / parsing.
 */

#ifndef C_GSNAP_H_
#define C_GSNAP_H_

#include "goby/C_CompactHelpers.h"

#ifdef __cplusplus
extern "C" {
#endif
    void gobyGsnap_parse(CAlignmentsWriterHelper *writerHelper, char *alignment);
    void gobyGsnap_test_registerTargets(CAlignmentsWriterHelper *writerHelper, char *targets);
    void gobyGsnap_destoryAlignment(CAlignmentsWriterHelper *writerHelper);
#ifdef __cplusplus
}
#endif

#endif /* C_GSNAP_H_ */
