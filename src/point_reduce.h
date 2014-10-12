#include <assert.h>
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

float get_thresh_weight(float thresh, float dist)
{
    assert(thresh > 0);
    
    if (dist > thresh)
        return 1.0f;
    else
        return dist / thresh;
}

float get_cdf_weight(std::vector<double>& cdf_vec, int max_cdf_key, float dist)
{
    float half_dist = dist * 0.5f;
    int cdf_key = ceil(half_dist);
    //double intpart = 0.0;
    //cdf_key = (modf(half_min_dist, &intpart)>=0.5) ? ceil(half_min_dist) : floor(half_min_dist);

    if (cdf_key >= max_cdf_key)
        return 1.0f;
    else {
        float cdf_val = 0.5 + 0.5 * cdf_vec[cdf_key];
        if (false) {
            printf("cdf_key: %d, raw cdf val: %.5f, cdf_val: %.5f, ", cdf_key, cdf_vec[cdf_key], cdf_val);
        }
        return cdf_val;
    }
}

bool add_point(int max_pt_id, int kc, 
               std::vector<bool>& chosen_pts, std::vector<int>& selected_ids, int& num_pts_chosen, 
               std::vector<int>& covered_times, std::vector<float>& covered_mass, std::vector<Point>& points,
               std::vector<bool>& covered_imgs, int& num_img_covered, const int num_img)
{
    chosen_pts[max_pt_id] = true;
    selected_ids.push_back(max_pt_id);
    num_pts_chosen++;

    int img_id = 0;
    float ip_prob = 0.0f;

    for (unsigned int i = 0; i < points[max_pt_id].m_imgs.size(); i++) {
        img_id = points[max_pt_id].m_imgs[i];
        ip_prob = points[max_pt_id].m_probs[i];

        covered_times[img_id]++;
        covered_mass[img_id] += ip_prob;

        if (!covered_imgs[img_id] && covered_mass[img_id] >= kc){
            covered_imgs[img_id] = true;
            num_img_covered++;
            if (num_img_covered % 500 == 0){
                printf("Covered %d images (%f%%)...\n", 
                       num_img_covered, (float)num_img_covered * 100.0 / num_img);
                fflush(stdout);
            }
        }
    }
    return true;
}

bool add_point_pb(int max_pt_id, int kc, float min_prob, float max_pt_prob, std::vector<double>& covered_prob,
                  SpRowMat& prob_mat, SpRowMat& covered_mat,
                  std::vector<bool>& chosen_pts, std::vector<int>& selected_ids, int& num_pts_chosen, 
                  std::vector<int>& covered_times, std::vector<float>& covered_mass, std::vector<Point>& points,
                  std::vector<bool>& covered_imgs, int& num_img_covered, const int num_img)
{
    chosen_pts[max_pt_id] = true;
    selected_ids.push_back(max_pt_id);
    num_pts_chosen++;

    int img_id = 0;
    float ip_prob = 0.0f;

    for (unsigned int i = 0; i < points[max_pt_id].m_imgs.size(); i++) {
        img_id = points[max_pt_id].m_imgs[i];

        points[max_pt_id].m_probs[i] = points[max_pt_id].m_probs[i] * max_pt_prob;
        ip_prob = points[max_pt_id].m_probs[i];
                
        covered_prob[img_id] += prob_mat.coeffRef(img_id, kc - 1) * ip_prob;
        if (!covered_imgs[img_id]){
            getSurvFunc(kc, covered_times[img_id], max_pt_id, img_id,
                        ip_prob, covered_mat, prob_mat,
                        true, true);
            /*
            //double max_pb_val = 0.0;
            for (SpRowMat::InnerIterator it (prob_mat, i); it; ++it){
            dirty_vec[it.col()] = true;
            max_gains[it.col()] -= ip_prob;
            //if (it.value() > max_pb_val)
            //    max_pb_val = it.value();
            }
            //max_pb_vals[img_id] = max_pb_val;
            */
        }
        covered_mat.coeffRef(img_id, max_pt_id) = ip_prob;
        covered_times[img_id]++;
        covered_mass[img_id] += points[max_pt_id].m_probs[i];

        if (!covered_imgs[img_id] && covered_prob[img_id] >= min_prob){
            covered_imgs[img_id] = true;
            num_img_covered++;

            printf("Covered image %d with prob %f, %d out of %d covered so far (%f%%)...\n", 
                   img_id, covered_prob[img_id], num_img_covered, num_img, (float)num_img_covered*100.0f/num_img);
        }
    }

    return true;
}

bool reduce_redundancy(int kc, std::vector<bool>& chosen_pts, std::vector<int>& selected_ids, int& num_pts_chosen,
                       std::vector<int>& covered_times, std::vector<float>& covered_mass, std::vector<Point>& points)
{
    for (int i = 0; i < (int)selected_ids.size(); i++){
        bool redundant = true;
        int exam_pt_id = selected_ids[i];
        int img_id = 0;
        float ip_prob = 0.0f;
        
        for (unsigned int j = 0; j < points[exam_pt_id].m_imgs.size(); j++){
            img_id = points[exam_pt_id].m_imgs[j];
            ip_prob = points[exam_pt_id].m_probs[j];

            if (covered_mass[img_id] - ip_prob < kc){
                redundant = false;
                break;
            }
        }

        if (redundant) {
            for (unsigned int j = 0; j < points[exam_pt_id].m_imgs.size(); j++){
                img_id = points[exam_pt_id].m_imgs[j];
                ip_prob = points[exam_pt_id].m_probs[j];

                covered_times[img_id]--;
                covered_mass[img_id] -= ip_prob;
            }
            printf("removing point %d for redundancy.\n", exam_pt_id);
            num_pts_chosen--;
            chosen_pts[exam_pt_id] = false;
            selected_ids.erase(selected_ids.begin()+i);
            i--;
        }
    }
    return true;
}

                       
