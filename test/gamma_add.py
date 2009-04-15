
import sys, os

import pygsl
import pygsl.sf

while "python" not in os.listdir("."):
    os.chdir("..")

sys.path.append("python")
import spidir

from rasmus.common import *
from test import *

if os.system("xpdf") != 0:
    rplot_set_viewer("display")


def sample_xs(alphas, betas):
    return [random.gammavariate(a, 1/b) for a, b in zip(alphas, betas)]

def prod(lst):
    p = 1.0
    for i in lst:
        p *= i
    return p

def gamma_add_pdf2(z, params):
    a, b = params
    n = len(a)

    coeff = 1.0
    for i in xrange(n-1):
        coeff *= b[i]**a[i] / gamma(a[i])
    
    def func(xs):
        y = sum(xs[:-1])
        if y >= z:
            return 0.0
        return gammaPdf(z - y, (a[-1], b[-1])) * \
               exp(-sum(x*u for x,u in zip(xs, b[:-1]))) * \
               prod(x**(u-1) for x,u in zip(xs, a[:-1]))

    step = .02
    xpts = list(frange(0.01, 3.0, step))
    m = len(xpts)

    def integrate(xs):
        if len(xs) < n - 1:
            return sum(
                integrate(xs + [xpts[i]])
                for i in xrange(m))
        else:
            #print xs
            return step**(n-1) * func(xs)
        
    p = coeff * integrate([])
    print z, p
    return p


def gamma_add_pdf(y, params, nterms=None, tol=.02):

    if y <= 0.0:
        return 0.0

    if nterms is None:
        nterms = 100

    a, b = params
    b = [1.0 / i for i in b] # convert definition of beta
    b1 = min(b)


    C = prod((b1/v)**u for u, v in zip(a,b))
    p = sum(a)

    prob = 0.0
    gammas = []    # originally 1-indexed
    deltas = [1.0] # zero-indexed
    for i in xrange(nterms+1):
        k = i + 1
        gammas.append(sum(
            u*(1.0-b1/v)**k/k
            for u, v in zip(a,b)))

        k = i
        deltas.append(1.0 / (k+1) *
                      sum(j * gammas[j-1] * deltas[k+1-j]
                          for j in xrange(1,k+2)))
        
        prob2 = (deltas[k] * y**(p+k-1) * exp(-y/b1) /
                 (gamma(p+k) * b1**(p+k)))

        prob += prob2
        if prob2 / prob < tol:
            print i
            break        
    return C * prob

    
    #return C * sum(deltas[k] * y**(p+k-1) * exp(-y/b1) /
    #               (gamma(p+k) * b1**(p+k))
    #               for k in xrange(nterms))
    

class TestGammaAdd (unittest.TestCase):

    def test_add(self):
        """Test adding of independent gammas"""

        prep_dir("test/output/gamma_add")

        def make_plot(name, data, params, nterms):
            
            x, y = distrib(data, 40)            
            y2 = [[gamma_add_pdf(i, params, n, .02) for i in x]
                  for n in nterms]
            

            rplot_start("test/output/gamma_add/%s.pdf" % name)
            rplot("plot", x, y, t="l")
            for i, n in enumerate(nterms):
                k = i/float(len(nterms))
                rp.lines(x, y2[i], col=rp.rgb(1-k, k, 0))
            rplot_end(True)

        def make_plot2(name, data, params, tol=.01):
            
            x, y = distrib(data, 40)
            n = len(params[0])
            y2 = [[spidir.gammaSumPdf(i, n,
                                      params[0],
                                      params[1], tol)
                   for i in x]]
            

            rplot_start("test/output/gamma_add/%s.pdf" % name)
            rplot("plot", x, y, t="l")
            for i in xrange(len(y2)):
                k = i/float(len(y2))
                rp.lines(x, y2[i], col=rp.rgb(1-k, k, 0))
            rplot_end(True)

        '''
        alphas = [1.0, 2.0]
        betas = [1.2, 2.2]
        params = (alphas, betas)
        data = [sum(sample_xs(alphas, betas))
                for i in xrange(10000)]
        make_plot("2", data, params, [None])#, 5, 7, 10, 20])

        alphas = [1.0, 5.0, 3.0]
        betas = [1.0, 2.0, 3.0]
        params = (alphas, betas)
        data = [sum(sample_xs(alphas, betas))
                for i in xrange(10000)]
        make_plot("3", data, params, [None])#, 5, 7, 10, 20])
        '''
        
        alphas = [1.0, 5.0, 3.0, 16.0]
        betas = [1.0, 2.0, 3.0, 7.0]
        params = (alphas, betas)
        data = [sum(sample_xs(alphas, betas))
                for i in xrange(10000)]
        make_plot("4", data, params, [None])#, 5, 7, 10, 20, 40, 50])
        

        alphas = [1.0, 5.0, 3.0, 17.0]
        betas = [1.0, 2.0, 3.0, 8.0]
        params = (alphas, betas)
        data = [sum(sample_xs(alphas, betas))
                for i in xrange(10000)]
        make_plot2("c", data, params, .02)

        
if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())
