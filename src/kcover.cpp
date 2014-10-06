/* kccs.cpp */
#include "reduce.h"
#include "kccs.h"

int main(int argc, char **argv)
{
    if (argc < 8 || argc > 9) {
        printf("Usage: %s <input bundle file> <output bundle file> "
               "<point idx file> <k cover> <percentage> "
               "<min prob (<0 for additive)>  <img pt edge weight> "
               "[weighted img_pt_mat file]\n", argv[0]);
        return 1;
    }

    const char *input_bundle_file = argv[1];
    const char *output_bundle_file = argv[2];
    const char *pt_idx_path = argv[3];
    //int num_img = atoi(argv[3]);
    //int num_pts = atoi(argv[4]);
    int kc = atoi(argv[4]);
    float pct = atof(argv[5]);

    float min_prob = atof(argv[6]);
    bool only_add = (min_prob < 0);

    float ip_weight =  atof(argv[7]);

    const char *img_pt_matw_file = NULL;
    bool use_mat = false;
    if (argc >= 9){
        use_mat = true;
        img_pt_matw_file = argv[8];
    }
    
    omp_set_num_threads(8);
    //Eigen::setNbThreads(8);

    // Read input bundle file
    FILE *f;    
    char buf[256];

    f = fopen(input_bundle_file, "r");
    if (f == NULL) {
        perror ("Error opening bundler file %s!", input_bundle_file);
        return -1;
    }

    int line_num = 1;
    int num_img = 0, num_pts = 0;
    int start_line = 1000;
    
    while(fgets(buf, 256, f)) {
        /* Remove trailing newline */
        if (buf[strlen(buf) - 1] == '\n')
            buf[strlen(buf) - 1] = 0;

        if (line_num == 1){
            line_num++;            
            continue;
        } else if (line_num == 2){
            sscanf(buf, "%d %d", &num_img, &num_pts);
            start_line = 2 + num_img * 5;
        } else if (line_num > start_line) {
            if ((line_num - start_line - 1) % 3 == 0){
                
            }
            
                
        }

        sscanf(buf, "%d %d %g\n", &img_id, &pt_id, &dum);
        img_id -= 1;
        pt_id -= 1;

        points[pt_id].m_imgs.push_back(img_id);
        points[pt_id].m_probs.push_back(dum);
        line_num++;
    }
    fclose(f);



    // Initialize points
    time_t start_time, finish_time;
    std::vector<Point> points;
    points.reserve(num_pts);
    for (int i = 0; i < num_pts; i++) {
        Point pt;
        points.push_back(pt);
    }
    printf("Finished initializing points.\n");
    fflush(stdout);    

    // Reading weighted image point matrix
    if (use_mat) {
        std::vector<T> tripletList;
        SpMat ip_mat(num_pts, num_img);

        printf("Start reading image point matrix...\n");
        fflush(stdout);    

        f = fopen(img_pt_matw_file, "r");

        if (f == NULL) {
            perror ("Error opening image point file %s!", img_pt_matw_file);
            return -1;
        }

        time(&start_time);
        int img_id = 0, pt_id = 0;
        float dum = 0.0f;
        while(fgets(buf, 256, f)) {
            /* Remove trailing newline */
            if (buf[strlen(buf) - 1] == '\n')
                buf[strlen(buf) - 1] = 0;

            sscanf(buf, "%d %d %g\n", &img_id, &pt_id, &dum);
            img_id -= 1;
            pt_id -= 1;

            points[pt_id].m_imgs.push_back(img_id);
            points[pt_id].m_probs.push_back(dum);
        }
        fclose(f);
        time(&finish_time);
        printf("Finished reading image point matrix, time used: %.3f seconds\n",
               difftime(finish_time, start_time));
        fflush(stdout);
    }
    
    // k cover algorithm to initiate
    printf("\nStarting K-cover algorithm...\n");
    fflush(stdout);    

    // image statistics
    std::vector<bool> covered_imgs; covered_imgs.reserve(num_img);
    std::vector<int> covered_times; covered_times.reserve(num_img);
    std::vector<float> covered_mass; covered_mass.reserve(num_img);
    std::vector<double> covered_prob; covered_prob.reserve(num_img);
    for (int i = 0; i < num_img; i++){
        covered_imgs.push_back(false);
        covered_times.push_back(0);
        covered_mass.push_back(0.0f);
        covered_prob.push_back(0.0);
    }

    // point statistics
    std::vector<bool> chosen_pts; chosen_pts.reserve(num_pts);
    std::vector<double> max_gains; max_gains.reserve(num_pts);
    for (int i = 0; i < num_pts; i++){
        chosen_pts.push_back(false);
        max_gains.push_back(0.0f);
    }

    for (int i = 0; i < num_pts; i++){
        for (unsigned int j = 0; j < points[i].m_probs.size();j++){
            max_gains[i] += points[i].m_probs[j];
            if (points[i].m_probs[j] > 1.0f){
                printf("pt %d, img %d, prob %f", i, j, points[i].m_probs[j]);
            }
        }
    }
    printf("First 3 max_gains: %f %f %f\n", max_gains[0], max_gains[1], max_gains[2]);
    fflush(stdout);

    // image-point relations
    SpMat covered_mat(num_img, num_pts); // current ip_mat
    SpMat prob_mat(num_img, num_pts); // current P(m) distribution

    // round statistics
    int num_img_covered = 0, num_pts_chosen = 0;
    int max_pt_id = -1;
    float max_cover = 0.0f, ip_prob = 0.0f;
    
    printf("Finished initialization, starting first covering...\n");
    fflush(stdout);

    time(&start_time);
    while (num_img_covered < num_img){
        max_cover = 0.0f;
        max_pt_id = -1;

        for (int i = 0; i < num_pts; i++){
            if (chosen_pts[i] || max_gains[i] <= max_cover)
                continue;

            double new_cover = 0.0;
            for (unsigned int j = 0; j < points[i].m_imgs.size(); j++) {
                if (!covered_imgs[points[i].m_imgs[j]])
                    new_cover += points[i].m_probs[j];
            }
            max_gains[i] = new_cover;

            if (new_cover > max_cover){
                max_cover = new_cover;
                max_pt_id = i;
            }
        }
        
        if(max_pt_id < 0) {
            printf("Adding more points won't cover new images, exiting.\n");
            break;
        }

        // adding point max_pt_id
        chosen_pts[max_pt_id] = true;
        num_pts_chosen++;

        for (unsigned int i = 0; i < points[max_pt_id].m_imgs.size(); i++) {
            img_id = points[max_pt_id].m_imgs[i];
            ip_prob = points[max_pt_id].m_probs[i];

            covered_times[img_id]++;
            covered_mass[img_id] += ip_prob;
            covered_mat.coeffRef(img_id, max_pt_id) = ip_prob * ip_weight;

            if (!covered_imgs[img_id] && covered_mass[img_id] >= kc){
                covered_imgs[img_id] = true;
                num_img_covered++;
                if (num_img_covered % 500 == 0){
                    printf("Covered %d images (%f%%)...\n", 
                           num_img_covered, (float)num_img_covered*100.0/num_img);
                    fflush(stdout);
                }
            }
        }

        if (num_pts_chosen % 1000 == 0){
            printf("%d out of %d images covered (kc=%d, %f%%), "
                   "number of chosen points: %d (%f%%)\n",
                   num_img_covered, num_img, kc, (float)num_img_covered*100.0/num_img,
                   num_pts_chosen, (float)num_pts_chosen*100.0/num_pts);
        }
    }
    time(&finish_time);
    printf("Finished 1st covering, time used: %.3f seconds, %d out of %d images covered (kc=%d), "
           "number of chosen points: %d (%f%%)\n",
           difftime(finish_time, start_time), num_img_covered, num_img, kc,
           num_pts_chosen, (float)num_pts_chosen*100.0/num_pts);
    fflush(stdout);

    //////////////////////////////////////////////////////////////////////
    // Stage 2: adding points according to survival function (prob. of
    // an image seeing more than kc points)
    //////////////////////////////////////////////////////////////////////
    if (!only_add) {
        num_img_covered = 0;
        for (int i = 0; i < num_img; i++) {
            double surv = getSurvFunc(kc, covered_times[i], 0, i,
                                      0.0f, covered_mat, prob_mat, 
                                      false, true);
            covered_prob[i] = surv;
            covered_imgs[i] = (surv >= min_prob);

            if (covered_imgs[i]){
                num_img_covered++;
                //printf("Image %d covered, covered times: %d, covered prob: %f\n", i, covered_times[i], surv);
            } //else {
            //printf("Image %d not covered, covered times: %d, covered prob: %f\n", i, covered_times[i], surv);
            //}
        
            if ((i+1) % 500 == 0){
                printf("Evaluated survival functions of %d images (%f%%)...\n", 
                       i+1, (float)(i+1)*100.0/num_img);
                fflush(stdout);
            }
        }
        printf("\nAfter initial covering, %d images are covered with min prob %f (%f%%).\n", 
               num_img_covered, min_prob, (float)num_img_covered*100.0/num_img);
        printf("Starting 2nd stage...\n");
        fflush(stdout);

        for (int i = 0; i < num_pts; i++){
            max_gains[i] = 0.0;
        }
        printf("Finished resetting score array max_gains...\n");
        fflush(stdout);
    
        time(&start_time);
        double tar_num_img = pct * num_img;
        while (num_img_covered < tar_num_img){
#pragma omp parallel for
            for (int i = 0; i < num_pts; i++){
                if (chosen_pts[i])
                    continue;
            
                double new_cover = 0.0;
                for (unsigned int j = 0; j < points[i].m_imgs.size(); j++) {
                    int loc_img_id = points[i].m_imgs[j];

                    if (!covered_imgs[loc_img_id]){
                        new_cover += prob_mat.coeffRef(loc_img_id, kc-1) * points[i].m_probs[j];
                        /*
                          double surv = getSurvFunc(kc, covered_times[img_id], i, img_id,
                          points[i].m_probs[j], covered_mat, prob_mat, 
                          true, false);
                          double delta = surv - covered_prob[img_id];

                          double surv = getSurvFuncDirect(kc, covered_times[img_id], i, img_id,
                          points[i].m_probs[j], covered_mat, true);
                          double delta2 = surv - covered_prob[img_id];
                          printf("delta new: %.5e\ndelta original: %.5e\n", delta, delta2);
                        */
                    }
                }
                max_gains[i] = new_cover;
            }
            /*
              for (int i = 0; i < num_pts; i++){
              if (chosen_pts[i])
              continue;
              printf("%d: %.5e\n", i, max_gains[i]);
              }
              exit(0);
            */

            max_cover = 0.0f;
            max_pt_id = -1;
            for (int i = 0; i < num_pts; i++){
                if (chosen_pts[i])
                    continue;

                if (max_gains[i] > max_cover){
                    max_cover = max_gains[i];
                    max_pt_id = i;
                }
            }
        
            //printf("max_cover: %f, max_pt_id: %d\n", max_cover, max_pt_id);
            //fflush(stdout);
            if(max_pt_id < 0) {
                printf("Adding more points won't increase probabilities, exiting.\n");
                break;
            }

            // adding point max_pt_id
            chosen_pts[max_pt_id] = true;
            num_pts_chosen++;

            for (unsigned int i = 0; i < points[max_pt_id].m_imgs.size(); i++) {
                img_id = points[max_pt_id].m_imgs[i];
                ip_prob = points[max_pt_id].m_probs[i] * ip_weight;

                //covered_mat.coeffRef(img_id, max_pt_id) = ip_prob;
                //covered_times[img_id]++;
                //covered_mass[img_id] += points[max_pt_id].m_probs[i];
                //covered_prob[img_id] = getSurvFunc(kc, covered_times[img_id], 0, img_id,
                //                                   0.0f, covered_mat, prob_mat, 
                //                                   false, true);
                //double surv = getSurvFuncDirect(kc, covered_times[img_id], max_pt_id, img_id, 
                //ip_prob, covered_mat, false);
                //printf("surv new: %.5e, surv original: %.5e, difference: %.5e\n", 
                //       covered_prob[img_id], surv, covered_prob[img_id]-surv);
                covered_prob[img_id] += prob_mat.coeffRef(img_id, kc - 1) * ip_prob;
                if (!covered_imgs[img_id]){
                    getSurvFunc(kc, covered_times[img_id], max_pt_id, img_id,
                                ip_prob, covered_mat, prob_mat,
                                true, true);
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

            if (num_pts_chosen % 1000 == 0){
                printf("%d out of %d images covered (kc=%d, %f%%), "
                       "number of chosen points: %d (%f%%)\n",
                       num_img_covered, num_img, kc, (float)num_img_covered * 100.0 / num_img,
                       num_pts_chosen, (float)num_pts_chosen*100.0/num_pts);
            }
        }
        time(&finish_time);

        printf("Finished 2nd covering, time used: %.3f seconds, %d out of %d images covered (min_prob=%f, kc=%d), "
               "number of chosen points: %d (%f%%)\n",
               difftime(finish_time, start_time), num_img_covered, num_img, min_prob, kc,
               num_pts_chosen, (float)num_pts_chosen*100.0/num_pts);
    } else {
        printf("Skipping 2nd covering...");
    }
    
    
    FILE *pt_idx_f = fopen(pt_idx_path, "w");
    if (pt_idx_f == NULL){
        printf("[reduce.cpp] Error opening file %s for writing.\n", pt_idx_path);
        return -1;
    }

    for (int i = 0; i < num_pts; i++){
        if (chosen_pts[i])
            fprintf(pt_idx_f, "%d\n", i);
    }
    fclose(pt_idx_f);

    printf("Finished writing %s\n", pt_idx_path);
    fflush(stdout);

    printf("All done.\n");
    fflush(stdout);
}

// adding update_prob (false true): initial round computing all P(m)
// adding update_prob (true false): incremental rounds evaluation
// adding update_prob (true true): incremental rounds addition
double getSurvFunc(int kc, int cov_count, int pt_idx, int img_idx,
                   float p, SpMat& covered_mat, SpMat& prob_mat,
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

                for (SpMat::InnerIterator it(covered_mat, img_idx); it; ++it){
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
        for (SpMat::InnerIterator it(prob_mat, img_idx); it; ++it){
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
                         float p, SpMat& covered_mat, bool adding)
{
    int N = cov_count;
    if (adding)
        N++;

    if (N < kc)
        return 0.0f;
    else if (N == kc){
        double prod = 1.0;
        for (SpMat::InnerIterator it(covered_mat, img_idx); it; ++it) {
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

        for (SpMat::InnerIterator it(covered_mat, img_idx); it; ++it){
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


