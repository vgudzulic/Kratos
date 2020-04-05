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

// System includes
#include <string.h>

// External includes
#include "zlib.h"

// Project includes
#include "utilities/zlib_utilities.h"

namespace Kratos
{
namespace ZlibUtilities
{

std::string CompressString(const std::string& rInputString)
{
    std::string str;

    // zlib struct
    const std::size_t data_size = rInputString.size()+1;
    std::size_t compressed_data_size = data_size;
    unsigned char aux_char[data_size];
    strcpy(reinterpret_cast<char*>(aux_char), rInputString.c_str());
    unsigned char* p_compressed_data = new unsigned char[compressed_data_size];

    const int result = compress2(p_compressed_data, &compressed_data_size, aux_char, data_size, Z_BEST_COMPRESSION);

    if (result == Z_OK) {
       str = (const char*)p_compressed_data;
    } else {
        KRATOS_ERROR << "Compression failed" << std::endl;
    }

    delete [] p_compressed_data;

    return str;
}

/***********************************************************************************/
/***********************************************************************************/

std::string DecompressString(const std::string& rInputString)
{
    std::string str;

    // zlib struct
    const std::size_t data_size = rInputString.size()+1;
    std::size_t decompressed_data_size = data_size;
    unsigned char aux_char[data_size];
    strcpy(reinterpret_cast<char*>(aux_char), rInputString.c_str());
    unsigned char* p_decompressed_data = new unsigned char[decompressed_data_size];

    const int result = uncompress(p_decompressed_data, &decompressed_data_size, aux_char, data_size);

    if (result == Z_OK) {
        str = (const char*)p_decompressed_data;
    } else {
        KRATOS_ERROR << "Decompression failed" << std::endl;
    }

    delete [] p_decompressed_data;

    return str;
}

} // namespace ZlibUtilities
} // namespace Kratos

#undef REAL
