/* reduce.h */
#define NUM_THREADS 16

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <vector>
#include <list>
#include <omp.h>

#include <flann/flann.hpp>

int knn_search_rnd_kd_forest(flann::Matrix<float>& qmat, int offset, int num_pts, int dim, int nn,
                             flann::Index<flann::L2<float> >& index, int check_num,
                             std::vector<float>& error_vec, int& changed, unsigned int* assignments)
{
    int j, k;
    flann::Matrix<float> query(new float[dim * num_pts], num_pts, dim);
    flann::Matrix<int> indices(new int[ nn * num_pts], num_pts, nn);
    flann::Matrix<float> dists(new float[nn * num_pts], num_pts, nn);

    for (k = 0; k < num_pts; k++){
        for (j = 0; j < dim; j++) {
            query[k][j] = qmat[offset + k][j];
        }
    }

    // knn search
    flann::SearchParams para(check_num);
    para.cores = NUM_THREADS;
    index.knnSearch(query, indices, dists, nn, para);

    for (k = 0; k < num_pts; k++){
        if (offset + k != indices[k][0])
            changed++;
        
        error_vec.push_back(dists[k][1]);
        assignments[offset + k] = indices[k][1];
    }
        
    delete[] query.ptr();
    delete[] indices.ptr();
    delete[] dists.ptr();
    
    return 0;
}

class Point {
public:
    std::vector<int> m_imgs;
    std::vector<float> m_probs;
    std::vector<int> m_pos;
    std::vector<unsigned char*> m_desc;
    float* m_mean_desc_db;
    unsigned char* m_mean_desc;
    int update_times;
    int m_pt_id;

    float min_sq_dist_val;
    int min_dist_id;
    int last_check_num;
    
    Point(int pt_id)
    {
        m_pt_id = pt_id;
        update_times = 0;

        m_mean_desc = NULL;
        //m_mean_desc = new unsigned char[128];
        //for (int j = 0; j < 128; j++)
        //    m_mean_desc[j] = 0;

        m_mean_desc_db = new float[128];
        for (int j = 0; j < 128; j++)
            m_mean_desc_db[j] = 0.0f;

        min_sq_dist_val = INFINITY;
        min_dist_id = 0;
        last_check_num = 0;
    }

    bool init_mean()
    {
        delete[] m_mean_desc_db;
        m_mean_desc_db = NULL;
        
        m_mean_desc = new unsigned char[128];
        for (int j = 0; j < 128; j++)
            m_mean_desc[j] = 0;
        return true;
    }
    

    /* either use update_mean + convert_mean */
    bool update_mean(unsigned char* desc)
    {
        float factor = 1.0f / (update_times + 1.0f);
        //printf("in update_mean, factor: %.2f\n", factor);
        //fflush(stdout);
        
        for (int j = 0; j < 128; j++){
            m_mean_desc_db[j] = (1-factor) * m_mean_desc_db[j] + factor * desc[j];
            //double intpart = 0.0;
            //m_mean_desc[j] = (unsigned char)(modf(m_mean_desc_db[j],&intpart)>=.5)?ceil(m_mean_desc_db[j]):floor(m_mean_desc_db[j]);            
        }
        update_times++;

        return true;
    }
    
    bool convert_mean()
    {
        if (m_mean_desc != NULL){
            delete[] m_mean_desc;
        }
        m_mean_desc = new unsigned char[128];
        for (int j = 0; j < 128; j++){
            double intpart = 0.0;
            m_mean_desc[j] = (unsigned char)(modf(m_mean_desc_db[j],&intpart)>=.5)?ceil(m_mean_desc_db[j]):floor(m_mean_desc_db[j]);
        }
        
        delete[] m_mean_desc_db;
        m_mean_desc_db = NULL;
        return true;
    }
    
    /* or use compute_mean */
    bool compute_mean(std::vector<float>& diff_vec)
    {
        bool verbose = true;
        //bool verbose = false;
        int min_record = 2;

        if (m_desc.empty())
            return false;

        if (m_mean_desc == NULL){
            m_mean_desc = new unsigned char[128];
        }
        
        int desc_count = (int)m_desc.size();
        if (verbose)
            printf("desc_count: %d\n", desc_count);

        for (int i = 0; i < desc_count; i++){
            if (verbose) {
                printf("img: %d; pos: %d\n", m_imgs[i], m_pos[i]);
                printf("desc %d:", i);            
            }
            for (int j = 0; j < 128; j++){
                m_mean_desc_db[j] += (float)m_desc[i][j];
                if (verbose)
                    printf(" %d", (int)m_desc[i][j]);
            }
            if (verbose)
                printf("\n");
        }

        // compute mean
        for (int j = 0; j < 128; j++){
            double intpart = 0.0;
            m_mean_desc_db[j] /= desc_count;
            m_mean_desc[j] = (unsigned char)(modf(m_mean_desc_db[j],&intpart)>=.5)?ceil(m_mean_desc_db[j]):floor(m_mean_desc_db[j]);
        }
        
        if (verbose) {
            printf("mean:");
            for (int j = 0; j < 128; j++){
                printf(" %.2f", m_mean_desc_db[j]);
            }
            printf("\n");        

            printf("mean_desc:");
            print_desc(m_mean_desc);
        }
        
        if (desc_count >= min_record) {
            for (int i = 0; i < desc_count; i++){
                float dist_mean = get_sift_float_dist(m_desc[i], m_mean_desc_db);
                diff_vec.push_back(dist_mean);
                printf("%d-th descriptor's distance to mean: %f\n", i, dist_mean);
            }
        }

        delete[] m_mean_desc_db;
        m_mean_desc_db = NULL;
        return true;
    }

    float get_pct_closer_to(std::vector<Point>& points, std::vector<int>& idx_list)
    {
        int total_desc_count = (int)m_desc.size();
        int total_pt_count = (int)idx_list.size();
        int closer_to_count = 0;
        /*
        std::vector<bool> closer_vec;
        for (int i = 0; i < total_desc_count; i++) {
            closer_vec.push_back(false);
        }
        */
#pragma omp parallel for
        for (int i = 0; i < total_desc_count; i++) {
            float self_dist = get_sift_sq_dist(m_desc[i], m_mean_desc);

            float min_partner_dist = INFINITY;
            for (int j = 0; j < total_pt_count; j++){
                float partner_dist = get_sift_sq_dist(m_desc[i], points[j].m_mean_desc);
                if (min_partner_dist > partner_dist)
                    min_partner_dist = partner_dist;
            }

            if(min_partner_dist < self_dist){
                #pragma omp atomic
                closer_to_count++;
            }
            
        }

        return (float)closer_to_count/total_desc_count;
    }

    float get_min_dist_to_points(std::vector<Point>& points, std::vector<int>& idx_list)
    {
        int list_len = (int) idx_list.size();
        std::vector<float> dist_vec;
        dist_vec.reserve(list_len - last_check_num);

        //printf("last_check_num: %d, list_len: %d ,min_sq_dist_val: %f, min_dist_id: %d, ",
        //       last_check_num, list_len, min_sq_dist_val, min_dist_id);
        for (int i = last_check_num; i < list_len; i++){
            dist_vec.push_back(INFINITY);
        }

#pragma omp parallel for        
        for (int i = last_check_num; i < list_len; i++){
            int pt_id = idx_list[i];
            if (pt_id == m_pt_id)
                continue;
            
            dist_vec[i-last_check_num] = get_sift_sq_dist(m_mean_desc, points[pt_id].m_mean_desc);
        }

        for (int i = last_check_num; i < list_len; i++){
            int offset = i - last_check_num;
            
            if (dist_vec[offset] < min_sq_dist_val){
                min_sq_dist_val = dist_vec[offset];
                min_dist_id = idx_list[i];
            }
        }
        dist_vec.clear();

        last_check_num = list_len;
        return sqrt(min_sq_dist_val);
    }
    
    float get_min_desc_dist(Point& partner)
    {
        float min_dist = INFINITY;

        for (unsigned int i = 0; i < partner.m_desc.size(); i++) {
            float new_dist = get_sift_dist(partner.m_desc[i], m_mean_desc);
            if (new_dist < min_dist)
                min_dist = new_dist;
        }
        return min_dist;
    }

    float get_max_desc_dist(Point& partner)
    {
        float max_dist = 0.0f;

        for (unsigned int i = 0; i < partner.m_desc.size(); i++) {
            float new_dist = get_sift_dist(partner.m_desc[i], m_mean_desc);
            if (new_dist > max_dist)
                max_dist = new_dist;
        }

        return max_dist;
    }
    
    float get_mean_dist(Point& partner)
    {
        return get_sift_dist(m_mean_desc, partner.m_mean_desc);
    }
    
    float get_desc_dist(Point& partner, int self_id, int partner_id)
    {
        if (self_id < 0 || self_id >= (int)m_desc.size())
            return -1.0;
        if (partner_id < 0 || partner_id >= (int) partner.m_desc.size())
            return -1.0;

        unsigned char* self_desc = m_desc[self_id];
        unsigned char* partner_desc = partner.m_desc[partner_id];        

        return get_sift_dist(self_desc, partner_desc);
    }

    float get_sift_float_dist(unsigned char* desc, float* tar)
    {
        float dist = 0.0;
        for (int i = 0; i < 128; i++){
            dist += (desc[i] - tar[i]) * (desc[i] - tar[i]);
        }
        dist = sqrt(dist);
        return dist;
    }

    float get_sift_dist(unsigned char* desc1, unsigned char* desc2)
    {
        float dist = get_sift_sq_dist(desc1, desc2);
        return sqrt(dist);
    }

    float get_sift_sq_dist(unsigned char* desc1, unsigned char* desc2)
    {
        float dist = 0.0;
        for (int i = 0; i < 128; i++){
            dist += (desc1[i] - desc2[i]) * (desc1[i] - desc2[i]);
        }

        return dist;
    }

    void print_desc(unsigned char* desc)
    {
        for (int j = 0; j < 128; j++){
            printf(" %d", (int)desc[j]);
        }
        printf("\n");
    }

    void fprint_desc(unsigned char* desc, FILE* f)
    {
        for (int j = 0; j < 128; j++){
            fprintf(f, " %d", (int)desc[j]);
        }
        fprintf(f, "\n");
    }

    void fprint_bin_desc(unsigned char* desc, FILE* f)
    {
        fwrite(desc, sizeof(unsigned char), 128, f);
    }

    void print_desc_db(float* desc)
    {
        for (int j = 0; j < 128; j++){
            printf(" %.3f", desc[j]);
        }
        printf("\n");
    }

    void print_mean()
    {
        if (m_mean_desc_db != NULL){
            printf("mean (float):");
            print_desc_db(m_mean_desc_db);
        }
        
        if (m_mean_desc != NULL) {
            printf("mean (unsigned char):");
            print_desc(m_mean_desc);
        }
    }

    void print_all_desc()
    {
        if (m_desc.empty())
            return;
            
        int total_desc_count = (int)m_desc.size();
        
        for (int i = 0; i < total_desc_count; i++) {
            printf("desc %d:",i);
            print_desc(m_desc[i]);
        }
    }

    void fprint_mean(FILE* f)
    {
        if (m_mean_desc != NULL) {
            fprint_desc(m_mean_desc, f);
        }
    }

    void fprint_bin_mean(FILE* f)
    {
        if (m_mean_desc != NULL) {
            fprint_bin_desc(m_mean_desc, f);
        }
    }
};

class PointSet 
{
public:
    int m_card_sum;
    std::vector<int> m_pt_id_vec;
    std::vector<int> m_img_cover;
    std::list<int> m_img_union;

    PointSet()
    {
        m_card_sum = 0;
    }

    bool contains(int pt_idx)
    {
        for (int i = 0; i < (int) m_pt_id_vec.size(); i++){
            if (m_pt_id_vec[i] == pt_idx)
                return true;
        }
        return false;
    }
    
    void recordHist(std::vector<int>& hist)
    {
        m_img_cover.clear();
        std::list<int>::iterator it;

        for (it = m_img_union.begin(); it != m_img_union.end(); it++){
            m_img_cover.push_back(hist[*it]);
        }
    }

    void printUnion()
    {
        std::list<int>::iterator it;
        printf("Current union:");

        for (it = m_img_union.begin(); it != m_img_union.end(); it++){
            printf(" %d", *it);
        }
        printf("\n");
    }
    
    void printVector(std::vector<int>& new_list)
    {
        std::vector<int>::iterator it = new_list.begin();
        printf("New list:");
        
        for (; it != new_list.end(); it++){
            printf(" %d", *it);
        }
        printf("\n");
    }

    bool addPtIdx(int pt_idx, std::vector<int>& new_list)
    {
        m_card_sum += new_list.size();
        m_pt_id_vec.push_back(pt_idx);
        //printf("Adding point index %d, old list size %d, new list length %d.\n", 
        //       pt_idx, (int)m_img_union.size(), (int)new_list.size());
        //fflush(stdout);
        updateUnion(new_list);

        return true;
    }

    bool addPtIdx(int pt_idx, std::vector<int>& new_list, std::vector<int>& hist)
    {
        addPtIdx(pt_idx, new_list);
        for (int i = 0; i < (int) new_list.size(); i++){
            //printf("%d %d\n", new_list[i], (int)m_img_cover.size());
            hist[new_list[i]]++;
        }

        return true;
    }

    void updateUnion(std::vector<int>& new_list)
    {
        int c1 = 0;
        int n1 = (int)new_list.size();
        std::list<int>::iterator it = m_img_union.begin();

        while (c1 < n1 && it != m_img_union.end()) {
            while (c1 < n1 && new_list[c1] < *it){
                m_img_union.insert(it, new_list[c1]);
                c1++;
            }
            
            if (c1 >= n1)
                break;

            if (new_list[c1] == *it) {
                c1++;
            }
            it++;
        }
        
        while (c1 < n1) {
            m_img_union.push_back(new_list[c1]);
            c1++;
        }
    }

    /*
    float getPointGain(int pt_idx, std::vector<int>& new_list)
    {
        if (contains(pt_idx))
            return 1.0f;
        if (m_img_union.empty())
            return 0.0f;
        float union_size = (float) m_img_union.size();
        return (float)getIntersectionSize(new_list)/union_size;
    }
    */

    int getUnionSize(std::vector<int>& new_list)
    {
        int inter_size = getIntersectionSize(new_list);
        return (int)new_list.size()+(int)m_img_union.size()-inter_size;
    }

    int getIntersectionSize(std::vector<int>& new_list)
    {
        int c1 = 0;
        int n1 = (int)new_list.size();

        std::list<int>::iterator it = m_img_union.begin();
        int inter_size = 0;

        while (c1 < n1 && it != m_img_union.end()) {
            while (c1 < n1 && new_list[c1] < *it)
                c1++;

            
            if (c1 >= n1)
                break;

            if (new_list[c1] == *it) {
                inter_size++;
            }
            it++;
        }

        return inter_size;
    }
    
};

class SuperPoint : public PointSet {
public:
    int m_sp_id;
    int m_pt_id;
    float m_score;
    float m_wscore;
    float m_min_inc;

    bool m_chosen;
    int m_new_cover;
    float m_new_gain;
    float m_gain_per_pt;
    
    void init(int sp_id, int num_img)
    {
        m_sp_id = sp_id;
        m_score = 0.0f;
        m_wscore = 0.0f;
        m_min_inc = 0.0f;

        m_chosen = false;
    }
    
    // Initilize using a single point
    SuperPoint(int pt_idx, int sp_id, int num_img, std::vector<int>& new_list)
    {
        m_pt_id = pt_idx;
        init(sp_id, num_img);
        m_new_cover = (int)new_list.size();

        std::vector<int> hist;
        hist.reserve(num_img);
        for (int i = 0; i < num_img; i++) {
            hist.push_back(0);
        }

        addPtIdx(pt_idx, new_list, hist);
        recordHist(hist);
        hist.clear();
    }

    // Initialize using a record
    SuperPoint(const char* buf, int num_img, std::vector<Point>& points)
    {
        int sp_size = 0, offset = 0;
        int pt_idx = 0, idx_offset = 0;
        sscanf(buf, "<super point %d> score: %f, wscore: %f, card_sum: %d, min_inc: %f, pt_idx(%d):%n",
               &m_sp_id, &m_score, &m_wscore, &m_card_sum, &m_min_inc, &sp_size, &offset);

        m_wscore = m_score * log(m_card_sum);
        //m_wscore = m_score * m_card_sum;
        m_chosen = false;
        m_new_cover = m_card_sum;
        m_new_gain = m_wscore;
        
        std::vector<int> hist;
        hist.reserve(num_img);
        for (int i = 0; i < num_img; i++) {
            hist.push_back(0);
        }

        m_card_sum = 0;
        m_pt_id_vec.reserve(sp_size);
        bool first = true;
        while (1 == sscanf(buf + offset, " %d%n", &pt_idx, &idx_offset)) {
            addPtIdx(pt_idx, points[pt_idx].m_imgs, hist);
            offset += idx_offset;
            if (first){
                m_pt_id = pt_idx;
                first = false;
            }
        }
        recordHist(hist);
        hist.clear();
    }

    bool normalizeWScore(float max_wscore)
    {
        m_wscore /= max_wscore;
        m_new_gain = m_wscore;
        return true;
    }
    
    bool updateWMinInc(float cur_score)
    {
        float k = (float)m_pt_id_vec.size() + 1.0f;
        m_min_inc = k * sqrt((float)m_img_union.size() * cur_score) - m_card_sum;

        //printf("current wmin_inc: %.2f\n", m_min_inc);
        return true;
    }
    
    bool updateMinInc(float cur_score)
    {
        float k = (float)m_pt_id_vec.size() + 1.0f;
        m_min_inc = k * (float)m_img_union.size() * cur_score - m_card_sum;

        //printf("current min_inc: %.2f\n", m_min_inc);
        return true;
    }

    bool resetMinInc()
    {
        m_min_inc = 0.0f;
        return true;
    }

    bool addPtIdxUpScore(int pt_idx, std::vector<int>& new_list, float score, float wscore)
    {
        m_score = score;
        m_wscore = wscore;
        return addPtIdx(pt_idx, new_list);
    }
    
    float getNewPointCoverGain(std::vector<bool>& covered_pts)
    {
        int num_uncovered = 0;
        for (int i = 0; i < (int) m_pt_id_vec.size(); i++){
            if (!covered_pts[m_pt_id_vec[i]])
                num_uncovered++;
        }
        m_new_gain = m_wscore * num_uncovered / m_pt_id_vec.size();
        return m_new_gain;
    }
    
    int getNewCoverNum(std::vector<bool>& covered_imgs, std::vector<Point>& points)
    {
        int new_cover = 0;
        int pos = 0;
        
        std::list<int>::iterator it;        
        for (it = m_img_union.begin(); it != m_img_union.end(); it++,pos++){
            if (!covered_imgs[*it])
                new_cover += m_img_cover[pos];
        }

        m_new_cover = new_cover;
        return new_cover;
    }
    
    void fprintStatus(FILE* f)
    {
        fprintf(f, "<super point %d> score: %.2f, wscore: %.2f, card_sum: %d, min_inc: %.2f, pt_idx(%d):",
                m_sp_id, m_score, m_wscore, m_card_sum, m_min_inc, (int)m_pt_id_vec.size());

        for (int i = 0; i < (int) m_pt_id_vec.size(); i++)
            fprintf(f, " %d", m_pt_id_vec[i]);
        fprintf(f, "\n");
    }
    
    void printStatus()
    {
        printf("<super point %d> score: %.2f, wscore: %.2f, card_sum: %d, min_inc: %.2f, pt_idx(%d):",
               m_sp_id, m_score, m_wscore, m_card_sum, m_min_inc, (int)m_pt_id_vec.size());

        for (int i = 0; i < (int) m_pt_id_vec.size(); i++)
            printf(" %d", m_pt_id_vec[i]);
        printf("\n");
    }
    
    bool getWScore(std::vector<int>& new_list, float& score, float& wscore, int& card_sum)
    {
        if (new_list.size() < m_min_inc){
            score = 0.0f;
            wscore = 0.0f;
            return false;
        }
            
        int k = m_pt_id_vec.size()+1;
        card_sum = m_card_sum + new_list.size();

        int union_size = getUnionSize(new_list);

        /* printUnion();
        printVector(new_list);
        printf("Union size: %d\n", union_size);
        updateUnion(new_list);
        printUnion();
        exit(0); */

        if (union_size == 0){
            printf("union set is empty!");
            score = 0.0f;
            wscore = 0.0f;
            return false;
        }
        
        score = (float)card_sum / (float)(union_size * k);
        wscore = score * card_sum / k;
        return true;
    }
};
