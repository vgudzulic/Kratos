//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// Project includes
#include "testing/testing.h"
#include "utilities/zlib_utilities.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ZlibUtilities, KratosCoreFastSuite)
{
    const std::string str = "Lorem ipsum dolor sit amet, consectetur adipiscing elit";
    const std::string compressed_str = ZlibUtilities::CompressString(str);
    const std::string decompressed_str = ZlibUtilities::DecompressString(compressed_str);
    KRATOS_CHECK_STRING_EQUAL(str, decompressed_str);
    KRATOS_CHECK_LESS(compressed_str.size(), str.size());
}

}   // namespace Testing
}  // namespace Kratos.
