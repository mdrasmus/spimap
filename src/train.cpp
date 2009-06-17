/*=============================================================================

    SPIDIR    
    Parameter Estimation (Training)
    
    train.cpp
    started: Fri Mar 27 14:22:58 EDT 2009

=============================================================================*/

// c++ headers
#include <math.h>
#include <time.h>

// 3rd party
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multimin.h>

// spidir headers
#include "common.h"
#include "Matrix.h"
#include "seq_likelihood.h"
#include "parsimony.h"
#include "spidir.h"
#include "Tree.h"


namespace spidir {



class RatesEM
{
public:

    RatesEM(int ntrees, int nspecies, int nrates, 
            int *gene_sizes,
            float **lengths, float *times,
            float *sp_alpha, float *sp_beta, 
            float gene_alpha, float gene_beta) : 
        ntrees(ntrees),
        nspecies(nspecies), 
        gene_sizes(gene_sizes),
        lengths(lengths),
        times(times),
        sp_alpha(sp_alpha),
        sp_beta(sp_beta),
        gene_alpha(gene_alpha),
        gene_beta(gene_beta),
        nrates(nrates),
        gtab(ntrees, nrates),
        pgtab(ntrees, nrates)
    {
        // allocate gene rate optizer
	sol_gene_rat = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

        // setup optimizer for gene rates
	opt_gene_rate.function = &gene_rate_f;
        opt_gene_rate.params = this; 

        // setup optimizer for species rates
        const int ndim = 2;
        sol_sp_rate = gsl_multimin_fdfminimizer_alloc(
            gsl_multimin_fdfminimizer_vector_bfgs2, ndim);
        opt_sp_rate.f = &sp_rate_f;
        opt_sp_rate.df = &sp_rate_df;
        opt_sp_rate.fdf = &sp_rate_fdf;
        opt_sp_rate.n = ndim;
        opt_sp_rate.params = this;
    }


    ~RatesEM()
    {
	//gsl_root_fdfsolver_free(sol_gene_rat);
	gsl_root_fsolver_free(sol_gene_rat);
        gsl_multimin_fdfminimizer_free(sol_sp_rate);
    }

    // gene rates function and derivative
  
    static double gene_rate_f(double x, void *params)
    {
        RatesEM *em = (RatesEM*) params;
        double gene_nu = x;
	double q = 1.0 / gene_nu;

        // clamp gamma params
        if (gene_nu < .001)
            gene_nu = .001;

        double sum = 0.0;
        for (int j=0; j<em->ntrees; j++) {
            for (int k=0; k<em->nrates; k++) {
                double g = em->gtab[j][k];
                double lngamma = gammalog(g, q, q);
                double gamma = exp(lngamma);

                if (!isnan(gamma) && gamma != 0.0) {
                    sum += em->pgtab[j][k] * 
                        gammaDerivV(g, gene_nu) / gamma;
                }
            }
        }

        return sum;
    }

    static double gene_rate_df(double x, void *params)
    {
        RatesEM *em = (RatesEM*) params;
        double gene_nu = x;
	double q = 1.0 / gene_nu;

        // clamp gamma params
        if (gene_nu < .001)
            gene_nu = .001;

        double sum = 0.0;
        for (int j=0; j<em->ntrees; j++) {
            for (int k=0; k<em->nrates; k++) {
                double g = em->gtab[j][k];
                double lngamma = gammalog(g, q, q);
                double gamma = exp(lngamma);
		double dgamma = gammaDerivV(g, gene_nu);

                if (!isnan(gamma) && gamma != 0.0) {
                    sum += em->pgtab[j][k] * 
                        (gamma * gammaDerivV2(g, gene_nu) - dgamma*dgamma)
			 / (gamma*gamma);
                }
            }
        }

        return sum;
    }


    static void gene_rate_fdf(double x, void *params, 
                              double *f, double *df)
    {
        RatesEM *em = (RatesEM*) params;
        double gene_nu = x;
        double q = 1.0 / gene_nu;
   
        
        // clamp gamma params
        if (gene_nu < .001)
            gene_nu = .001;

        //printf("try %f %f\n", gene_alpha, gene_beta);

        double sum = 0.0;
        double sum2 = 0.0;
        for (int j=0; j<em->ntrees; j++) {
            for (int k=0; k<em->nrates; k++) {
                double g = em->gtab[j][k];
                double lngamma = gammalog(g, q, q);
                double gamma = exp(lngamma);
		double dgamma = gammaDerivV(g, gene_nu);
                
                // TODO: substitute a better test
                if (gamma != 0.0) {
		    sum += em->pgtab[j][k] * 
                        gammaDerivV(g, gene_nu) / gamma;
		    sum2 += em->pgtab[j][k] * 
                        (gamma * gammaDerivV2(g, gene_nu) - dgamma*dgamma)
			 / (gamma*gamma);
                }
            }
        }

        // set return
        *f = sum;

        *df = sum2;

        //printf("gene f = %f; fdf = (%f, %f); a=%f, b=%f\n", 
        //       sum,
        //       alpha_sum / em->nrates, 
        //       beta_sum / em->nrates,
        //       gene_alpha, gene_beta);
    }
    
    
    
    //======================================================
    // species rates function and derivative

    // f(a_i, b_i) = - sum_j sum_k pgtab_jk 
    //                       log(NB(c_ij, a_i, b_i  / (g_jk d_j t_i + b_i)))
    static double sp_rate_f(const gsl_vector *x, void *params)
    {
        RatesEM *em = (RatesEM*) params;
        double sp_alpha_i = gsl_vector_get(x, 0);
        double sp_beta_i = gsl_vector_get(x, 1);
        int i = em->cur_species;

        double sum = 0.0;
        for (int j=0; j<em->ntrees; j++) {
            for (int k=0; k<em->nrates; k++) {
                double lngamma = gammalog(em->lengths[j][i], sp_alpha_i, 
                                sp_beta_i / (em->gtab[j][k] *
                                             em->times[i]));
                if (!isnan(lngamma))
                    sum += em->pgtab[j][k] * lngamma;
            }
        }

        return -sum;
    }
    
    // d/d a_i f(a_i, b_i) = - sum_j sum_k pgtab_jk 
    //                (NB'_r(c_ij, a_i, b_i  / (g_jk d_j t_i + b_i)) /
    //                    NB(c_ij, a_i, b_i  / (g_jk d_j t_i + b_i)))
    static void sp_rate_df(const gsl_vector *x, void *params, gsl_vector *df)
    {
        RatesEM *em = (RatesEM*) params;
        double sp_alpha_i = gsl_vector_get(x, 0);
        double sp_beta_i = gsl_vector_get(x, 1);
        int i = em->cur_species;

        double alpha_sum = 0.0;
        double beta_sum = 0.0;
        for (int j=0; j<em->ntrees; j++) {
            double l = em->lengths[j][i];
            for (int k=0; k<em->nrates; k++) {
                double g = em->gtab[j][k];
                double bgt = sp_beta_i / (g * em->times[i]);
                double lngamma = gammalog(l, sp_alpha_i, bgt);
                double gamma = exp(lngamma);

                if (gamma != 0.0) {
                    alpha_sum += em->pgtab[j][k] * 
                        gammaDerivA(l, sp_alpha_i, bgt) / gamma;
                    beta_sum += em->pgtab[j][k] *
                        gammaDerivB(l, sp_alpha_i, bgt) /
                        (gamma * g * em->times[i]);
                }
            }
        }

        gsl_vector_set(df, 0, -alpha_sum);
        gsl_vector_set(df, 1, -beta_sum);
    }

    // f(a_i, b_i) = - sum_j sum_k pgtab_jk 
    //                       log(NB(c_ij, a_i, b_i  / (g_jk d_j t_i + b_i)))
    // d/d a_i f(a_i, b_i) = - sum_j sum_k pgtab_jk 
    //                (NB'_r(c_ij, a_i, b_i  / (g_jk d_j t_i + b_i)) /
    //                    NB(c_ij, a_i, b_i  / (g_jk d_j t_i + b_i)))
    static void sp_rate_fdf(const gsl_vector *x, void *params, 
                              double *f, gsl_vector *df)
    {
        RatesEM *em = (RatesEM*) params;
        double sp_alpha_i = gsl_vector_get(x, 0);
        double sp_beta_i = gsl_vector_get(x, 1);  
        int i = em->cur_species;
   
        double sum = 0.0;
        double alpha_sum = 0.0;
        double beta_sum = 0.0;
        for (int j=0; j<em->ntrees; j++) {
            double l = em->lengths[j][i];
            for (int k=0; k<em->nrates; k++) {
                double g = em->gtab[j][k];
                double bgt = sp_beta_i / (g * em->times[i]);
                double lngamma = gammalog(l, sp_alpha_i, bgt);
                double gamma = exp(lngamma);

                sum += em->pgtab[j][k] * lngamma;
                
                if (gamma != 0.0) {
                    alpha_sum += em->pgtab[j][k] *
                        gammaDerivA(l, sp_alpha_i, bgt) / gamma;
                    beta_sum += em->pgtab[j][k] * 
                        gammaDerivB(l, sp_alpha_i, bgt) /
                        (gamma * g * em->times[i]);
                }
            }
        }

        // set return
        *f = -sum;

        //printf("sp %d, sum=%f, %f, %f\n", i, sum, alpha_sum, beta_sum);

        gsl_vector_set(df, 0, -alpha_sum);
        gsl_vector_set(df, 1, -beta_sum);
    }


    float likelihood()
    {
        double logl = 0.0;
	
        for (int j=0; j<ntrees; j++) {
            double sum = 0.0;
            for (int k=0; k<nrates; k++) {
                double prod = 0.0;
                for (int i=0; i<nspecies; i++)
                    prod += gammalog(lengths[j][i], sp_alpha[i], 
				     sp_beta[i] / (times[i] * gtab[j][k]));
                sum += pgtab[j][k] * exp(prod);
            }
            logl += log(sum);
        }

        return logl;
    }

    //====================

    void init_params()
    {
        for (int i=0; i<nspecies; i++) {
            sp_alpha[i] = 1.0;
            sp_beta[i] = 1.0;
        }

        gene_nu = 1.0; 
    }

    inline double gene_post(double g, double A, double B, double C)
    {
        return pow(g, B) * exp(A*g - C / g);
    }
    

    double find_upper_g(double m, double A, double B, double C,
                        double tol1=.05, double tol2=.01)
    {        
        double fm = gene_post(m, A, B, C);
        double top = 2*m;
        double bot = m;

        // extent top
        while (true) {
            double ftop = gene_post(top, A, B, C);
            if (ftop/fm <= tol2)
                break;
            top *= 2;
        }

        // binary search
        while (true) {
            double u = (top + bot) / 2.0;
            double fu = gene_post(u, A, B, C);

            if (fu / fm > tol1)
                bot = u;
            else if (fu / fm < tol2)
                top = u;
            else
                return u;
        }
    }


    double find_lower_g(double m, double A, double B, double C,
                       double tol1=.05, double tol2=.01)
    {
        double fm = gene_post(m, A, B, C);
        double top = m;
        double bot = 0;       

        // binary search
        while (true) {
            double u = (top + bot) / 2.0;
            double fu = gene_post(u, A, B, C);

            if (fu / fm > tol1)
                top = u;
            else if (fu / fm < tol2)
                bot = u;
            else
                return u;
        }
    }


    // populate gene rate posteriors
    void EStep()
    {
        // temp variables for PDF of posterior gene rate
        float x[nrates+1];
        double y[nrates+1];
	float q = 1.0 / gene_nu;
	double gene_alpha = q, gene_beta = q;

        // determine commonly used coefficients
        double A = - gene_beta;
        double B = gene_alpha - 1.0;
        for (int i=0; i<nspecies; i++)
            B -= sp_alpha[i];

        for (int j=0; j<ntrees; j++) {

            // determine commonly used coefficients
            double C = 0.0;
            for (int i=0; i<nspecies; i++)
                C += sp_beta[i] * lengths[j][i] / times[i];

            // find main range of gene rates
            double mid = (-B - sqrt(B*B - 4*A*C)) / (2*A);
            double top = find_upper_g(mid, A, B, C);
            double bot = find_lower_g(mid, A, B, C);

            int half_nrates = (nrates + 1) / 2;
            double step1 = (mid - bot) / half_nrates;
            double step2 = (top - mid) / (nrates + 1 - half_nrates);

            // compute x, y for posterior gene rate PDF
            for (int k=0; k<half_nrates; k++) {
                x[k] = bot + step1 * k;
                y[k] = gene_post(x[k], A, B, C);
            }
            for (int k=half_nrates; k<nrates+1; k++) {
                x[k] = mid + step2 * (k-half_nrates);
                y[k] = gene_post(x[k], A, B, C);
            }

            // compute gtab and pgtab
            double total = 0.0;
            for (int k=0; k<nrates; k++) {
                gtab[j][k] = (x[k] + x[k+1]) / 2.0;
                pgtab[j][k] = (y[k] + y[k+1]) * (x[k+1] - x[k]) / 2.0;
                total += pgtab[j][k];
            }

            // normalize pgtab[j]
            for (int k=0; k<nrates; k++) {
                pgtab[j][k] /= total;
            }
        }
    }
    

    // maximize each model parameter given the hidden data estimated from
    // last iteration
    void MStep()
    {
        // optimization config
        double step_size = .01;
        double tol = .1;
        const double epsabs = .01;
        gsl_vector *init_x = gsl_vector_alloc(2);
        int status;
        //double r = 1.0, r0 = 1.0;
	double low = .0001, high = gene_nu * 5;
	int iter = 0;
	const int maxiter = 10;
        
        // optimize gene rate parameters
	gsl_root_fsolver_set(sol_gene_rat, &opt_gene_rate, low, high);

	status = GSL_CONTINUE;
	for (iter=0; iter<maxiter && status==GSL_CONTINUE; iter++) {
            // do one iteration
	    status = gsl_root_fsolver_iterate(sol_gene_rat);
            if (status)
                break;

            // check convergence
	    low = gsl_root_fsolver_x_lower(sol_gene_rat);
	    high = gsl_root_fsolver_x_upper(sol_gene_rat);
	    status = gsl_root_test_interval(low, high, 0, 0.01);
	}
	gene_nu = gsl_root_fsolver_root(sol_gene_rat);
        fprintf(stderr, "nu = %f (iter=%d)\n", gene_nu, iter);
        
        
        // optimize each species rate parmater set
        for (int i=0; i<nspecies; i++) {
            cur_species = i;
            gsl_vector_set(init_x, 0, sp_alpha[i]);
            gsl_vector_set(init_x, 1, sp_beta[i]);
            gsl_multimin_fdfminimizer_set(sol_sp_rate, &opt_sp_rate, init_x, 
                                          step_size, tol);            
            status = GSL_CONTINUE;

            
	    for (iter=0; iter<maxiter && status==GSL_CONTINUE; iter++) {
                // do one iteration
                status = gsl_multimin_fdfminimizer_iterate(sol_sp_rate);
                if (status)
                    break;        
                // get gradient
                status = gsl_multimin_test_gradient(sol_sp_rate->gradient, epsabs);
            }

            fprintf(stderr, "calc lk\n");
            double lk = likelihood();
            fprintf(stderr, "species %d %d %f\n", i, iter, lk);

            sp_alpha[i] = gsl_vector_get(sol_sp_rate->x, 0);
            sp_beta[i] = gsl_vector_get(sol_sp_rate->x, 1);

            //printf("sp[%d] = (%f, %f)\n", i, sp_alpha[i], sp_beta[i]);
        }

        gsl_vector_free(init_x);
    }

    // data
    int ntrees;
    int nspecies;
    int *gene_sizes;
    float **lengths;

    // given fixed parameters
    float *times;

    // model parameters
    float *sp_alpha;
    float *sp_beta;
    float gene_alpha;
    float gene_beta;
    float gene_nu;
    
    //protected:

    // hidden data
    int nrates;
    Matrix<float> gtab;
    Matrix<double> pgtab;

    // optimizers
    gsl_multimin_fdfminimizer *sol_sp_rate;
    gsl_root_fsolver *sol_gene_rat;
    gsl_function opt_gene_rate;
    gsl_multimin_function_fdf opt_sp_rate;
    int cur_species;
};


extern "C" {

void train(int ntrees, int nspecies, int *gene_sizes,
           float **lengths, float *times,
           float *sp_alpha, float *sp_beta, 
           float *gene_alpha, float *gene_beta,
           int nrates, int max_iter)
{
    RatesEM em(ntrees, nspecies, nrates, gene_sizes, lengths, times,
               sp_alpha, sp_beta, *gene_alpha, *gene_beta);
    
    
    // make initial guess for model parameters
    em.init_params();
    
    // iterate until convergence
    for (int iter=0; iter<max_iter; iter++) {
        em.EStep();
        em.MStep();
    }
      
}


RatesEM *allocRatesEM(int ntrees, int nspecies, int nrates,
                      int *gene_sizes,
                      float **lengths, float *times,
                      float *sp_alpha, float *sp_beta, 
                      float gene_alpha, float gene_beta)
{
    int *gene_sizes2 = new int [ntrees];
    for (int j=0; j<ntrees; j++)
        gene_sizes2[j] = gene_sizes[j];

    // copy lengths
    float **lengths2 = new float* [ntrees];
    for (int j=0; j<ntrees; j++) {
        lengths2[j] = new float [nspecies];
        for (int i=0; i<nspecies; i++)
            lengths2[j][i] = lengths[j][i];
    }

    float *times2 = new float [nspecies];
    float *sp_alpha2 = new float [nspecies];
    float *sp_beta2 = new float [nspecies];

    for (int i=0; i<nspecies; i++) {
        times2[i] = times[i];
        sp_alpha2[i] = sp_alpha[i];
        sp_beta2[i] = sp_beta[i];
    }

    return new RatesEM(ntrees, nspecies, nrates, gene_sizes2,
                       lengths2, times2,
		       sp_alpha2, sp_beta2, gene_alpha, gene_beta);
}


void freeRatesEM(RatesEM *em)
{
    delete [] em->gene_sizes;

    for (int j=0; j<em->ntrees; j++)
        delete [] em->lengths[j];
    delete [] em->lengths;

    delete [] em->times;
    delete [] em->sp_alpha;
    delete [] em->sp_beta;
    
    delete em;
}


void RatesEM_Init(RatesEM *em)
{
    em->init_params();
}


void RatesEM_EStep(RatesEM* em)
{
    em->EStep();
}

void RatesEM_MStep(RatesEM* em)
{
    em->MStep();
}

float RatesEM_likelihood(RatesEM *em)
{
    return em->likelihood();
}


void RatesEM_getParams(RatesEM *em, float *params)
{
    params[0] = 1.0 / em->gene_nu; //em->gene_alpha;
    params[1] = 1.0 / em->gene_nu; //em->gene_beta;

    for (int i=0; i<em->nspecies; i++) {
        params[2+2*i] = em->sp_alpha[i];
        params[2+2*i+1] = em->sp_beta[i];
    }
}



} // extern "C"



} // spidir
