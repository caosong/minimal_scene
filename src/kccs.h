#include <string>
#include <iostream>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::RowMajor;
using Eigen::ColMajor;
typedef Eigen::Triplet<double> T;
typedef Eigen::Triplet<int> TI;
typedef Eigen::SparseMatrix<int,ColMajor> SpColIntMat;
typedef Eigen::SparseMatrix<double,ColMajor> SpColMat;
typedef Eigen::SparseMatrix<double,RowMajor> SpRowMat;

double getSurvFuncDirect(int kc, int cov_count, int pt_idx, int img_idx,
                         float p, SpRowMat& covered_mat, bool adding)
{
    int N = cov_count;
    if (adding)
        N++;

    if (N < kc)
        return 0.0f;
    else if (N == kc){
        double prod = 1.0;
        for (SpRowMat::InnerIterator it(covered_mat, img_idx); it; ++it) {
            prod *= it.value();
        }
        if (adding) {
            prod *= p;
        }
        return prod;
    }

    double factor = 1.0 / (double)(N+1.0);
    double Q = 1.0 - (double)kc * factor;

    std::complex<double> sum(0.0, 0.0);
    for (int i = 1; i <= N; i++){
        std::complex<double> result(0.0, 0.0);
        std::complex<double> exp_1(0.0, -2.0*i*kc*M_PI/(double)(N+1));
        std::complex<double> exp_2(0.0, -2.0*i*M_PI/(double)(N+1));
        std::complex<double> exp_3(0.0, 2.0*i*M_PI/(double)(N+1));

        result = (1.0 - std::exp(exp_1)) / (1.0 - std::exp(exp_2));

        for (SpRowMat::InnerIterator it(covered_mat, img_idx); it; ++it){
            result *= it.value() * std::exp(exp_3) + (1.0 - it.value());
        }

        if (adding){
            result *= (double)p * std::exp(exp_3) + (double)(1.0 - p);
        }
        sum += result;
    }
    sum *= factor;
    //printf("sum real: %f, imag: %f\n", sum.real(), sum.imag());
    Q -= sum.real();
    return Q;
}

// adding update_prob (false true): initial round computing all P(m)
// adding update_prob (true false): incremental rounds evaluation
// adding update_prob (true true): incremental rounds addition
double getSurvFunc(int kc, int cov_count, int pt_idx, int img_idx,
                   float p, SpRowMat& covered_mat, SpRowMat& prob_mat,
                   bool adding, bool update_prob)
{
    int N = cov_count;
    if (adding)
        N++;

    // Survival function value
    double Q = 0.0; 

    // initial round: compute P(m) from scratch
    if (!adding && update_prob) {
        double factor = 1.0 / (double)(N+1.0);
        double new_prob_val = 0.0;
        
        for (int m = 0; m <= N; m++) {
            std::complex<double> sum(0.0, 0.0);
            for (int n = 0; n <= N; n++){
                double exp_val1 = -2.0 * M_PI * n * m * factor;
                double exp_val2 = 2.0 * M_PI * n * factor;
                
                //std::complex<double> exp_1(0.0, -2.0 * M_PI * n * m * factor);
                //std::complex<double> exp_2(0.0, 2.0 * M_PI * n * factor);
                //result = std::exp(exp_1);
                //exp_factor = std::exp(exp_2);
                std::complex<double> result(cos(exp_val1), sin(exp_val1));
                std::complex<double> exp_factor(cos(exp_val2), sin(exp_val2));

                for (SpRowMat::InnerIterator it(covered_mat, img_idx); it; ++it){
                    result *= it.value() * exp_factor + (1.0 - it.value());
                }

                //if (adding) {
                //    result *= (double)p * exp_factor + (1.0 - (double)p);
                //}
                sum += result;
            }
            sum *= factor;
            new_prob_val = sum.real();
            if (new_prob_val < 0.0){
                new_prob_val = 0.0;
            }

            prob_mat.coeffRef(img_idx, m) = new_prob_val;
            if (m >= kc){
                Q += new_prob_val;
            }
        }
    } 

    // incremental adding round
    else if (adding && update_prob) {
        //std::vector<double> new_prob_vec;
        //new_prob_vec.reserve(N+1);
        //for (int i = 0; i <= N; i++){
        //    new_prob_vec.push_back(0.0);
        //}
        double prob = 0.0, last_prob = 0.0;

        //printf("previous distribution:");
        for (SpRowMat::InnerIterator it(prob_mat, img_idx); it; ++it){
            //printf(" %d:%.5e", it.col(), it.value());
            if (it.col() == 0){
                last_prob = it.value();
                prob_mat.coeffRef(img_idx, 0) = (1 - p) * it.value();
                //new_prob_vec[0] = (1 - p) * it.value();
                continue;
            }
            //prob = p * last_prob + (1 - p) * it.value();
            prob = it.value();
            prob_mat.coeffRef(img_idx, it.col()) = p * last_prob + (1 - p) * prob;
            last_prob = prob;
            //new_prob_vec[it.col()] = prob;

            if (it.col() >= kc)
                Q += it.value();
        }
        //printf("\np:%.5e\n", p);
        prob_mat.coeffRef(img_idx, N) = p * last_prob;
        //new_prob_vec[N] = prob;
        if (N >= kc)
            Q += prob;

        //printf("current distribution:");        
        //for (unsigned i = 0; i < new_prob_vec.size(); i++){
            //printf(" %d:%.5e", i, new_prob_vec.at(i));
            //prob_mat.coeffRef(img_idx, i) = new_prob_vec.at(i);
        //}
        //printf("\nQ:%.5e\n", Q);
        //new_prob_vec.clear();
    }

    /*
    // incremental testing round
    else if (adding && !update_prob) {
        double prob = 0.0;
        
        for (int m = kc; m <= N-1; m++) {
            prob = p * prob_mat.coeffRef(img_idx, m-1);
            prob += (1 - p) * prob_mat.coeffRef(img_idx, m);
            Q += prob;
        }
        Q += p * prob_mat.coeffRef(img_idx, N-1);
    }
    */
    return Q;
}


bool printSparseMat(SpRowMat& ip_mat, int num_show)
{
    int num_shown = 0;
    for (int k = 0; k < ip_mat.outerSize(); k++)
        for (SpRowMat::InnerIterator it(ip_mat,k); it; ++it)
            {
                printf(" (pt %d, img %d, %.5f)", it.row(), it.col(), it.value());
                num_shown++;
                if (num_shown >= num_show){
                    printf("\n");
                    return true;
                }
            }    
    return false;
}
