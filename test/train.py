
import sys, os

import pygsl
import pygsl.roots
import pygsl.sf

while "python" not in os.listdir("."):
    os.chdir("..")

sys.path.append("python")
import spidir

from rasmus.common import *
from test import *

if os.system("xpdf 2&> /dev/null") != 0:
    rplot_set_viewer("display")


def gamma_deriv_x(x, a, b):
    return b**a / gamma(a) * exp(-b*x) * pow(x, a-2) * (a - b*x - 1)


def gamma_deriv_a(x, a, b):
    return exp(-b*x) * pow(x,a-1) * pow(b, a) / gamma(a) * (
        log(b) + log(x) - pygsl.sf.psi_n(0, a)[0])


def gamma_deriv_b(x, a, b):
    return pow(x, a-1) / gamma(a) * exp(-b*x) * pow(b, a-1) * (a - b*x)


def invgammaPdf(x, params):
    a, b = params
    return exp(log(b)*a + log(x)*(-a - 1) + (-b / x) - gammaln(a))


class TestTrain (unittest.TestCase):


    def _test_matrix(self):
        m = spidir.c_matrix(spidir.c_int, [[1,2,3],[4,5,6]])
        print m[1][1]

    def test_train(self):

        # read data
        stree = read_tree("test/data/flies.norm.stree")
        print "stree"
        #dat = read_delim("test/data/flies-one2one.fastleaf.lens")

        dat = read_delim("test/data/flies-one2one.uniform.lens")
        
        species = dat[0]
        lens = map2(float, dat[1:])
        times = [node.dist for node in stree.postorder() if node.parent]
        gene_sizes = [1000] * len(lens)
        
        ntrees = len(lens)
        nrates = 20
        
        #params = spidir.train_params(gene_sizes, lens, times, species,
        #                             nrates=10, max_iter=10)

        print "alloc"
        em = spidir.alloc_rates_em(gene_sizes, lens, times, species, nrates)

        spidir.RatesEM_Init(em)
        print "estep"
        spidir.RatesEM_EStep(em)
        logl = spidir.RatesEM_likelihood(em)
        print "logl:", logl
        
        for i in xrange(4):
            spidir.RatesEM_EStep(em)
            spidir.RatesEM_MStep(em)
            print spidir.RatesEM_likelihood(em)

        print "get params"
        params = spidir.rates_em_get_params(em, species)
        pd(params)

        spidir.free_rates_em(em)

        #prep_dir("test/output/train.fastleaf/")
        #spidir.write_params("test/output/train.fastleaf/flies.fastleaf.out.params", params)
        
    
    def _test_gamma_deriv(self):
        """gamma derivatives"""

        dx = .001
        a = 2.0
        b = 3.0
        m = 2.0
        x = list(frange(.1, 10, .01))

        def gammav(x, v):
            return gammaPdf(x, (1.0/v, 1.0/v))

        def dgammav(x, v, dx):
            return (gammav(x, v+dx) - gammav(x, v)) / dx


        '''
        dy = [(gammaPdf(i+dx, (a,b)) - gammaPdf(i, (a,b))) / dx for i in x]
        dy2 = [gamma_deriv_x(i, a, b) for i in x]

        prep_dir("test/output/train/")
        
        rplot_start("test/output/train/gamma_dx.pdf")
        rplot("plot", x, dy, t="l")
        rp.lines(x, dy2, t="l", col="red")
        rplot_end(True)

        dy = [(gammaPdf(m, (i+dx,b)) - gammaPdf(m, (i,b))) / dx for i in x]
        dy2 = [gamma_deriv_a(m, i, b) for i in x]

        
        rplot_start("test/output/train/gamma_da.pdf")
        rplot("plot", x, dy, t="l")
        rp.lines(x, dy2, t="l", col="red")
        rplot_end(True)


        dy = [(gammaPdf(m, (a,i+dx)) - gammaPdf(m, (a,i))) / dx for i in x]
        dy2 = [gamma_deriv_b(m, a, i) for i in x]

        
        rplot_start("test/output/train/gamma_db.pdf")
        rplot("plot", x, dy, t="l")
        rp.lines(x, dy2, t="l", col="red")
        rplot_end(True)
        '''

        x = list(frange(.01, 1.5, .01))

        dy = [(gammav(m, i+dx) - gammav(m, i)) / dx for i in x]
        dy2 = [spidir.gammaDerivV(m, i) for i in x]
        
        rplot_start("test/output/train/gamma_dv.pdf")
        rplot("plot", x, dy, t="l")
        rp.lines(x, dy2, t="l", col="red")
        rplot_end(True)
        

        #dy = [(dgammav(m, i+dx, dx/2.0) - dgammav(m, i, dx/2.0))/dx
        #      for i in x]
        dy = [(spidir.gammaDerivV(m, i+dx) - spidir.gammaDerivV(m, i)) / dx
              for i in x]
        dy2 = [spidir.gammaDerivV2(m, i) for i in x]
        
        rplot_start("test/output/train/gamma_dv2.pdf")
        rplot("plot", x, dy, t="l")
        rp.lines(x, dy2, t="l", col="red")
        rplot_end(True)



    def _test_invgamma(self):
        """Test inverse gamma quantile function"""
        a = 2
        b = 3 #34

        x = list(frange(.001, 5, .01))
        y = [invgammaPdf(i, (a, b)) for i in x]
        y2 = [spidir.invgamma(i, a, b) for i in x]
        y3 = [spidir.invgammaCdf(i, a, b) for i in x]
        y4 = [spidir.quantInvgamma(i, a, b) for i in x]

        #print spidir.quantInvgamma(.5, a, b)
        #print spidir.gammaDerivA(145.166061, 1.000000, 1.000000), \
        #      spidir.gammaPdf(145.166061, 1.000000, 1.000000), \
        #      spidir.gammalog(145.166061, 1.000000, 1.000000)

        prep_dir("test/output/train-invgamma/")

        rplot_start("test/output/train-invgamma/invgamma.pdf")
        rplot("plot", x, y, t="l") #, ylim=[min(x), max(x)])
        rp.lines(x, y2, t="l", col="red")
        rp.lines(x, y3, t="l", col="blue")
        rp.lines(y4, x, t="l", col="green") # show overlay CDF
        rplot_end(True)



    def _test_invgamma_pdf(self):
        """Test inverse gamma PDF/CDF functions"""

        prep_dir("test/output/train-invgamma-pdf/")
        
        #a = 23.0
        #b = 53.0
        a = 3
        b = .5
        n = 10000

        data = [1 / random.gammavariate(a, 1.0/b) for i in xrange(n)]

        x, y = distrib(data, 200)
        y2 = [spidir.invgamma(i, a, b) for i in x]
        
        rplot_start("test/output/train-invgamma-pdf/invgamma-pdf.pdf")
        rplot("plot", x, y, t="l")
        rp.lines(x, y2, t="l", col="red")
        rplot_end(True)


        x, y = cdf(data)
        y2 = [spidir.invgammaCdf(i, a, b) for i in x]

        rplot_start("test/output/train-invgamma-pdf/invgamma-cdf.pdf")
        rplot("plot", x, y, t="l")
        rp.lines(x, y2, t="l", col="red")
        rplot_end(True)


        s = 2
        x = list(frange(0, 10, .1))
        y = [spidir.incompleteGammaC(s, i) for i in x]

        rplot_start("test/output/train-invgamma-pdf/incompletegammac.pdf")
        rplot("plot", x, y, t="l")
        #rp.lines(x, y2, t="l", col="red")
        rplot_end(True)



    def _test_estep(self):


        def gene_post(g, lens, times, params):
            (gene_a, gene_b), sp_params = params

            gene_len = 10
            changes = [gene_len * l for l in lens]

            A = prod((g*t*gene_len + b)**(-a-c)
                    for t, a, b, c in zip(times,
                                          util.cget(sp_params, 0),
                                          util.cget(sp_params, 1),
                                          changes))

            return g**(gene_a-1+sum(changes)) * exp(-gene_b*g) * A
        
        
        def deriv_log_gene_post(g, lens, times, params):
            (gene_a, gene_b), sp_params = params

            gene_len = 10
            changes = [gene_len * l for l in lens]

            A = sum((-a-c) * (t * gene_len) / (g * t * gene_len + b)
                    for t, a, b, c in zip(times,
                                          util.cget(sp_params, 0),
                                          util.cget(sp_params, 1),
                                          changes))

            return (gene_a - 1 + sum(changes)) / g - gene_b + A
            
        

        def gene_post_mode(lens, times, params):
            (gene_a, gene_b), sp_params = params
            
            def f(x, params2=None):
                print x
                return deriv_log_gene_post(x, lens, times, params)

            mu = gene_a / gene_b

            func = pygsl.gsl_function.gsl_function(f, [])
            sol = pygsl.roots.bisection(func)
            sol.set(mu*.1, mu*10)

            for i in range(20):
                sol.iterate()
            return sol.root()


        def find_upper_g(lens, times, params, tol1=.05, tol2=.01):
            
            m = gene_post_mode(lens, times, params)
            fm = gene_post(m, lens, times, params)
            top = 2*m
            bot = m

            # extent top
            while True:
                ftop = gene_post(top, lens, times, params)
                print top, ftop/fm
                if ftop/fm <= tol2:
                    break
                top *= 2

            # binary search
            while True:
                u = (top + bot) / 2.0
                fu = gene_post(u, lens, times, params)

                if fu / fm > tol1:
                    bot = u
                elif fu / fm < tol2:
                    top = u
                else:
                    break

                print bot, u, top

            return u


        def find_lower_g(lens, times, params, tol1=.05, tol2=.01):
            
            m = gene_post_mode(lens, times, params)
            fm = gene_post(m, lens, times, params)
            top = m
            bot = 0
            
            # binary search
            while True:
                u = (top + bot) / 2.0
                fu = gene_post(u, lens, times, params)

                if fu / fm > tol1:
                    top = u
                elif fu / fm < tol2:
                    bot = u
                else:
                    break

                print bot, u, top

            return u


        
        prep_dir("test/output/train-estep/")

        nrates = 20  
        lens = [.1, .2, .3]
        times = [1.0, 0.5, 0.1]
        params = ((90.0, .1), ((2.0, 3.0),
                               (2.0, 3.0),
                               (2.0, 3.0)))

        m = gene_post_mode(lens, times, params)
        top = find_upper_g(lens, times, params)
        ftop = gene_post(top, lens, times, params)
        bot = find_lower_g(lens, times, params)
        fbot = gene_post(bot, lens, times, params)
        print "range:", bot, m, top
        

        x = list(frange(.01, 2*m, 2*m/60.0))
        y = [gene_post(g, lens, times, params) for g in x]
        #y = map(log, y)
        #y2 = [deriv_log_gene_post(g, lens, times, params) for g in x]

        step1 = (m - bot) / (nrates / 2)
        step2 = (top - m) / (nrates / 2)
        box = list(frange(bot, m+step1/2, step1)) + \
              list(frange(m, top+step2/2, step2))
        fbox = [gene_post(g, lens, times, params) for g in box]
        
        rplot_start("test/output/train-estep/gene_post.pdf")
        rplot("plot", x, y, t="l")
        #rplot("lines", x, y2, t="l", col="blue")
        rplot("lines", [bot, bot], [0, fbot], col="green")
        rplot("lines", [top, top], [0, ftop], col="green")
        rplot("lines", box, fbox, col="blue")
        rplot("lines", box, fbox, col="blue", t="h")
        rplot("lines", [m, m], [0, max(y)], col="red")
        rplot_end(True)



    def _test_estep_old(self):


        def gene_post(g, lens, times, params):
            (gene_a, gene_b), sp_params = params

            A = gene_a - 1 - sum(util.cget(sp_params, 0))
            B = sum(b * l / t for b, l, t in zip(util.cget(sp_params, 1),
                                                 lens,
                                                 times))
            
            return g**A * exp(-gene_b * g - B/g)

        def deriv_gene_post(g, lens, times, params):
            (gene_a, gene_b), sp_params = params

            A = gene_a - 1 - sum(util.cget(sp_params, 0))
            B = sum(b * l / t for b, l, t in zip(util.cget(sp_params, 1),
                                                 lens,
                                                 times))
            ex = exp(-gene_b * g - B/g)
            
            return A*g**(A-1) * ex + g**A * ex * (-gene_b + B/g/g)
                   

        def deriv_log_gene_post(g, lens, times, params):
            (gene_a, gene_b), sp_params = params


            A = gene_a - 1 - sum(util.cget(sp_params, 0))
            B = sum(b * l / t for b, l, t in zip(util.cget(sp_params, 1),
                                                 lens,
                                                 times))
            
            return A/g - gene_b + B / (g*g)
        

        def gene_post_mode(lens, times, params):
            (gene_a, gene_b), sp_params = params

            A = - gene_b
            B = (gene_a - 1 - sum(util.cget(sp_params, 0)))
            C = sum(b * l / t for b, l, t in zip(util.cget(sp_params, 1),
                                                   lens,
                                                   times))
            print A, B, C, B*B - 4*A*C

            return (-B - sqrt(B*B - 4*A*C)) / (2*A)


        def find_upper_g(lens, times, params, tol1=.05, tol2=.01):
            
            m = gene_post_mode(lens, times, params)
            fm = gene_post(m, lens, times, params)
            top = 2*m
            bot = m

            # extent top
            while True:
                ftop = gene_post(top, lens, times, params)
                print top, ftop/fm
                if ftop/fm <= tol2:
                    break
                top *= 2

            # binary search
            while True:
                u = (top + bot) / 2.0
                fu = gene_post(u, lens, times, params)

                if fu / fm > tol1:
                    bot = u
                elif fu / fm < tol2:
                    top = u
                else:
                    break

                print bot, u, top

            return u


        def find_lower_g(lens, times, params, tol1=.05, tol2=.01):
            
            m = gene_post_mode(lens, times, params)
            fm = gene_post(m, lens, times, params)
            top = m
            bot = 0
            
            # binary search
            while True:
                u = (top + bot) / 2.0
                fu = gene_post(u, lens, times, params)

                if fu / fm > tol1:
                    top = u
                elif fu / fm < tol2:
                    bot = u
                else:
                    break

                print bot, u, top

            return u


        
        prep_dir("test/output/train-estep/")

        nrates = 10  
        lens = [.1, .2, .3]
        times = [1.0, 0.5, 0.1]
        params = ((90.0, .1), ((2.0, 3.0),
                               (2.0, 3.0),
                               (2.0, 3.0)))

        m = gene_post_mode(lens, times, params)
        top = find_upper_g(lens, times, params)
        ftop = gene_post(top, lens, times, params)
        bot = find_lower_g(lens, times, params)
        fbot = gene_post(bot, lens, times, params)
        print "range:", bot, m, top
        

        x = list(frange(.01, 5*m, 5*m/30.0))
        y = [gene_post(g, lens, times, params) for g in x]
        #y = map(log, y)
        y2 = [deriv_gene_post(g, lens, times, params) for g in x]

        step1 = (m - bot) / (nrates / 2)
        step2 = (top - m) / (nrates / 2)
        box = list(frange(bot, m+step1/2, step1)) + \
              list(frange(m, top+step2/2, step2))
        fbox = [gene_post(g, lens, times, params) for g in box]
        
        rplot_start("test/output/train-estep/gene_post.pdf")
        rplot("plot", x, y, t="l")
        #rplot("lines", x, y2, t="l", col="blue")
        rplot("lines", [m, m], [0, max(y)], col="red")
        rplot("lines", [bot, bot], [0, fbot], col="green")
        rplot("lines", [top, top], [0, ftop], col="green")
        rplot("lines", box, fbox, col="blue")
        rplot("lines", box, fbox, col="blue", t="h")
        rplot_end(True)

        

        
if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())
