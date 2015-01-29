# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org


import unittest

##
# Test base classes
class TestBase(unittest.TestCase):

    ##
    # General function for testing bases
    def baseTest(self, b, points):
        for t in points:
            val = b.eval(t[0], t[1], t[2])
            self.failUnlessAlmostEqual(val, t[3], msg = ("%f != %f => (%d, %d) @ %f"%(val, t[3], t[0], t[1], t[2])))

    def testLinear(self):
        from pysgpp.base import SLinearBase
        
        b = SLinearBase()
        
        points = [(1, 1, 0.5, 1.0),
                  (1, 1, 0.25, 0.5),
                  (2, 1, 0.25, 1.0),
                  (2, 1, 0.125, 0.5),
                  
                  ]

        self.baseTest(b, points)
        
    def testLinearBoundary(self):
        from pysgpp.base import SLinearBoundaryBase
        
        b = SLinearBoundaryBase()
        
        points = [(0, 0, 0.0, 1.0),
                  (0, 0, 1.0, 0.0),
                  (0, 0, 0.25, 0.75),
                  (0, 0, 0.75, 0.25),
                  (0, 1, 0.0, 0.0),
                  (0, 1, 1.0, 1.0),
                  (0, 1, 0.75, 0.75),
                  (0, 1, 0.25, 0.25),
                  (1, 1, 0.5, 1.0),
                  (1, 1, 0.25, 0.5),
                  (2, 1, 0.25, 1.0),
                  (2, 1, 0.125, 0.5),
                  
                  ]

        self.baseTest(b, points)             
        
    def testLinearModified(self):
        from pysgpp.base import SLinearModifiedBase
        
        b = SLinearModifiedBase()
        
        points = [(1, 1, 0.5, 1.0),
                  (1, 1, 0.25, 1.0),
                  
                  (2, 1, 0.25, 1.0),
                  (2, 1, 0.125, 1.5),
                  (2, 1, 0.375, 0.5),
                  
                  (2, 3, 0.75, 1.0),
                  (2, 3, 0.75 + 0.125, 1.5),
                  (2, 3, 0.75 - 0.125, 0.5),
                  
                  (3, 3, 0.375 + 0.0625, 0.5)
                  ]
        
        self.baseTest(b, points)

    def testLinearStretched(self):
        from pysgpp.base import SLinearStretchedBase
        
        b = SLinearStretchedBase()
        
# Point data format: C, L, R, value 
# (1,1) and (2,1) tried
        points = [(1.312631659, 0.5, 7, 0.125020255230769),
		  (0.96716821, 0.5, 7, 0.071872032307692)
                  ]
        
        self.baseTest(b, points)
        
        

#    def testPoly(self):
#        from pysgpp.base import SPolyBase
#        
#        self.failUnlessRaises(Exception, SPolyBase, 0)
#        
#        b = SPolyBase(2)
#        
#        points = [(1, 1, 0.5, 1.0),
#                  (1, 1, 0.25, 0.75),
#                  (1, 1, 0.75, 0.75),
#                  (2, 1, 0.25, 1.0),
#                  (2, 1, 0.125, 0.75),
#                  (2, 1, 0.25+0.125, 0.75),
#                  ]
# 
#        self.baseTest(b, points)
#        
#        b = SPolyBase(3)
#
#        points = [(1, 1, 0.5, 1.0),
#                  (1, 1, 0.25, 0.75),
#                  (1, 1, 0.75, 0.75),
#                  
#                  (2, 1, 0.25, 1.0),
#                  (2, 1, 0.125, 0.875),
#                  (2, 1, 0.25+0.125, 0.625),
#                  
#                  (3, 1, 0.0625, 0.875),
#                  (3, 1, 0.125+0.0625, 0.625),
#                  
#                  (3, 3, 0.375-0.0625, 0.625),
#                  (3, 3, 0.375+0.0625, 0.875),
#                  ]
# 
#        self.baseTest(b, points)
#
#    def testModPoly(self):
#        from pysgpp.base import SPolyModifiedBase
#        
#        self.failUnlessRaises(Exception, SPolyModifiedBase, -1)
#        
#        b = SPolyModifiedBase(0)
#        points = [(1, 1, 0.5, 1.0),
#                  (1, 1, 0.25, 1.0),
#                  (2, 1, 0.25, 1.0),
#                  (2, 1, 0.125, 1.0),
#                  (2, 3, 0.75, 1.0),
#                  (3, 1, 0.125, 1.0),
#                  ]
#
#        self.baseTest(b, points)
#
#        b = SPolyModifiedBase(1)
#        points = [(1, 1, 0.5, 1.0),
#                  (1, 1, 0.25, 1.0),
#                  
#                  (2, 1, 0.25, 1.0),
#                  (2, 1, 0.0, 2.0),
#                  (2, 1, 0.5, 0.0),
#                  
#                  (2, 3, 0.75, 1.0),
#                  (2, 3, 1.0, 2.0),
#                  (2, 3, 0.5, 0.0),
#
#                  (3, 1, 0.125, 1.0),
#                  (3, 1, 0.0, 2.0),
#                  (3, 1, 0.25, 0.0),
#                  
#                  (3, 3, 0.25 + 0.125, 1.0),
#                  (3, 3, 0.5, 2.0),
#                  (3, 3, 0.25, 0.0),
#                  ]
#
#        self.baseTest(b, points)
#
#        b = SPolyModifiedBase(2)
#        points = [(1, 1, 0.5, 1.0),
#                  (1, 1, 0.25, 1.0),
#                  
#                  (2, 1, 0.25, 1.0),
#                  (2, 1, 0.0, 2.0),
#                  (2, 1, 0.5, 0.0),
#                  
#                  (2, 3, 0.75, 1.0),
#                  (2, 3, 1.0, 2.0),
#                  (2, 3, 0.5, 0.0),
#                  
#                  (3, 1, 0.125, 1.0),
#                  (3, 1, 0.0, 2.0 + 2.0/3.0),
#                  (3, 1, 0.25, 0.0),
#
#                  (3, 3, 0.375, 1.0),
#                  (3, 3, 0.25, 0.0),
#                  (3, 3, 0.5, 0.0),
#                  (3, 3, (0.25+0.375)/2, 0.75)
#                  
#                  ]
#
#        self.baseTest(b, points)

        
class TestFunctions(unittest.TestCase):
    def testGetAffected(self):
        from pysgpp.base import HashGridIndex, HashGridStorage, SLinearBase
        from pysgpp.base import SGetAffectedBasisFunctions
        
        i = HashGridIndex(1)
        s = HashGridStorage(1)
        
        b = SLinearBase()
        
        i.set(0,1,1)
        s.insert(i)
        
        ga = SGetAffectedBasisFunctions(s)
        
        x = ga(b, [0.25])
        
        self.failUnlessEqual(x[0][0], 0)
        self.failUnlessEqual(x[0][1], 0.5)
        
        
    def testGetAffectedBoundaries(self):
        from pysgpp.base import HashGridIndex, HashGridStorage, SLinearBoundaryBase
        from pysgpp.base import SGetAffectedBasisFunctionsBoundaries
        
        i = HashGridIndex(1)
        s = HashGridStorage(1)
        
        b = SLinearBoundaryBase()

        i.set(0,0,0)
        s.insert(i)
        i.set(0,0,1)
        s.insert(i)        
        i.set(0,1,1)
        s.insert(i)
        
        ga = SGetAffectedBasisFunctionsBoundaries(s)
        
        x = ga(b, [0.5])
        
        self.failUnlessEqual(x[0][0], 0)
        self.failUnlessEqual(x[0][1], 0.5)
        self.failUnlessEqual(x[1][0], 1)
        self.failUnlessEqual(x[1][1], 0.5)
        self.failUnlessEqual(x[2][0], 2)
        self.failUnlessEqual(x[2][1], 1.0)   

    def testGetAffectedLinearStretched(self):
        from pysgpp.base import HashGridIndex, HashGridStorage, SLinearStretchedBoundaryBase
        from pysgpp.base import SGetAffectedBasisFunctionsLinearStretchedBoundaries
	from pysgpp.base import Stretching, Stretching1D, DimensionBoundary
        
	str1d = Stretching1D()
	str1d.type='log'
	str1d.x_0=1
	str1d.xsi=10
	dimBound = DimensionBoundary() 
 	dimBound.leftBoundary=0.5
	dimBound.rightBoundary=7
	stretch=Stretching(1,dimBound,str1d)

        i = HashGridIndex(1)
        s = HashGridStorage(1)
	s.setStretching(stretch)
        
        b = SLinearStretchedBoundaryBase()
        
        i.set(0,0,0)
        s.insert(i)
        i.set(0,0,1)
        s.insert(i)        
        i.set(0,1,1)
        s.insert(i)
        
        ga = SGetAffectedBasisFunctionsLinearStretchedBoundaries(s)
        
        x = ga(b, [0.25])
        
        self.failUnlessEqual(x[0][0], 0)
        self.failUnlessEqual(x[0][1], 1.0384615384615385)     

# Run tests for this file if executed as application 
if __name__=='__main__':
    unittest.main()
    
