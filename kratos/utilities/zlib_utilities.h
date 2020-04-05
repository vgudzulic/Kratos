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

#if !defined(KRATOS_ZLIB_UTILITIES)
#define KRATOS_ZLIB_UTILITIES

// System includes

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
    
/**
 * @namespace ZlibUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities using the library triangle
 * @author Vicente Mataix Ferrandiz
 */
namespace ZlibUtilities
{
    /**
     * @brief This method compresses the given string
     * @param rInputString The input string
     * @return The compressed string
     */
    std::string KRATOS_API(KRATOS_CORE) CompressString(const std::string& rInputString);

    /**
     * @brief This method decompresses the given string
     * @param rInputString The input string
     * @return The decompressed string
     */
    std::string KRATOS_API(KRATOS_CORE) DecompressString(const std::string& rInputString);

}; // namespace ZlibUtilities
}  // namespace Kratos
#endif /* KRATOS_ZLIB_UTILITIES defined */
