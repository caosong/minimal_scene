#include <assert.h>

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

                       
