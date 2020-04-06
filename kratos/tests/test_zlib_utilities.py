from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestZlibUtilities(KratosUnittest.TestCase):

    def test_slib_utilities(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        lorem_ipsum = "Lorem ipsum dolor sit amet, consectetur adipiscing elit"

        lorem_ipsum_compressed = KratosMultiphysics.ZlibUtilities.CompressString(lorem_ipsum)
        lorem_ipsum_decompressed = KratosMultiphysics.ZlibUtilities.DecompressString(lorem_ipsum_compressed)

        self.assertEqual(lorem_ipsum, lorem_ipsum_decompressed)
        self.assertLess(len(lorem_ipsum_compressed), len(lorem_ipsum))

if __name__ == '__main__':
    KratosUnittest.main()

