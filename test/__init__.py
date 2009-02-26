
import unittest, traceback, os, shutil


def fequal(f1, f2, rel=.0001):
    """assert whether two floats are approximately equal"""
    
    if f1 == f2:
        return
    
    if f2 == 0:
        err = f1
    else:
        err = abs(f1 - f2) / abs(f2)
    x = (err < rel)

    assert x, "%f != %f  [err=%f]" % (f1, f2, err)



def prep_dir(path):
    try:
        if os.path.exists(path):
            shutil.rmtree(path)
    except:
        pass
    else:
        os.makedirs(path)
    


class TestResult (unittest.TestResult):
    """Custom test result for writing result to terminal"""
    
    def __init__(self):
        unittest.TestResult.__init__(self)
        self.successes = 0

    def startTest(self, test):
        print "=" * 78
        print test.shortDescription()

    def addError(self, test, err):
        unittest.TestResult.addError(self, test, err)
        print "[ ERROR ]", test.shortDescription()
        traceback.print_exception(*err) #type(error), error, tracebk)
        
           
    def addFailure(self, test, err):
        unittest.TestResult.addFailure(self, test, err)        
        print "[ FAIL ]", test.shortDescription()
        traceback.print_exception(*err) #type(error), error, tracebk)
        
   
    def addSuccess(self, test):
        unittest.TestResult.addSuccess(self, test)
        print "ok."
        self.successes += 1
   
    
        
class TestRunner (unittest.TextTestRunner):
    """Custom test runner"""

    def __init__(self):
        unittest.TextTestRunner.__init__(self)

    def run(self, test):
        result = TestResult()
        test.run(result)

        print "=" * 78
        print "successes:", result.successes
        print "errors:   ", len(result.errors)
        print "failures: ", len(result.failures)

        return result

