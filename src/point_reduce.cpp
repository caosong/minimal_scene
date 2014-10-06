/* point_reduce.cpp */
#include "reduce.h"
#include "kccs.h"
#include "keys.h"
#include "point_reduce.h"

int main(int argc, char **argv)
{
    if (argc < 8 || argc > 20) {
        printf("Usage: %s <img_pt_mat file> <list keys file> "
               "<num images> <num points> "
               "<k cover> <percentage> "
               "<point idx file> [record dist] [use k-cover] "
               "[reduce memory] [ip weight] [min prob] [threshold as prob] "
               "[none 0; read mean 1; write mean 2] [mean file] [use binary] "
               "[load pt idx] [cdf file] [cdf as prob]\n", argv[0]);
        return 1;
    }

    const char *img_pt_mat_file = argv[1];
    const char *list_keys_file = argv[2];
    int num_img = atoi(argv[3]);
    int num_pts = atoi(argv[4]);
    int kc = atoi(argv[5]);
    float pct = atof(argv[6]);
    const char *pt_idx_path = argv[7];

    float min_dist_cur = 180.0f;

    bool record_dist = true;
    if (argc >= 9){
        record_dist = (atoi(argv[8]) > 0);
    }

    bool use_approx = false;
    bool use_direct = false;
    bool use_kc = false;
    if (argc >= 10){
        use_kc = (atoi(argv[9]) > 0);
    }

    bool use_rd = false;
    bool use_rm = false;
    if (argc >= 11){
        use_rm = (atoi(argv[10]) > 0);
    }
    
    float ip_weight = 1.0f;
    float min_prob = 0.99;
    if (argc >= 12){
        ip_weight = atof(argv[11]);
    }
    if (argc >= 13){
        min_prob = atof(argv[12]);
    }
    bool use_pb = (ip_weight > 0.0) && (min_prob > 0.0);

    bool thresh_as_prob = false;
    if (argc >= 14){
        thresh_as_prob = (atoi(argv[13]) > 0);
        if (thresh_as_prob)
            printf("use threshold discounted value as prob\n");
        else
            printf("use constant %.2f as prob\n", ip_weight);
    }
    
    bool write_mean = false, read_mean = false;
    if (argc >= 15){
        int mean_mode = atoi(argv[14]);
        if (mean_mode == 1){
            read_mean = true;
            printf("read existing mean file ");
        } else if (mean_mode == 2) {
            write_mean = true;
            printf("will write new mean file ");
        }
    }
    
    if ((read_mean || write_mean) && argc < 16){
        printf("mean file path required!\n");
        exit(1);
    }

    const char *mean_path = NULL;
    if (argc >= 16) {
        mean_path = argv[15];
        printf("%s\n", mean_path);
    }

    bool use_bin = false;
    if (argc >= 17) {    
        use_bin = (atoi(argv[16]) > 0);
    }

    if (read_mean || write_mean){
        if (use_bin)
            printf("using binary format\n");
        else
            printf("using ascii format\n");
    }

    bool use_load = false;
    if (argc >= 18){
        use_load = (atoi(argv[17]) > 0);
    }
    
    const char *cdf_path = NULL;
    bool use_cdf = false;
    if (argc >= 19){
        use_cdf = true;
        cdf_path = argv[18];
    }

    bool cdf_as_prob = false;
    if (argc >= 20){
        cdf_as_prob = (atoi(argv[19]) > 0);
    }

    FILE *f = NULL;
    char buf[512];
    
    /* Read the cdf file */
    std::vector<double> diff_cdf;
    int max_cdf_key = 0;

    if (use_cdf && cdf_path != NULL) {
        printf("Start reading difference cdf...\n");
        f = fopen(cdf_path, "r");
        if (f == NULL) {
            printf("Error opening file %s for reading\n", cdf_path);
            return 1;
        }
        diff_cdf.push_back(0.0);
        
        while (fgets(buf, 512, f)) {
            /* Remove trailing newline */
            if (buf[strlen(buf) - 1] == '\n')
                buf[strlen(buf) - 1] = 0;
        
            diff_cdf.push_back(atof(buf));
        }
        fclose(f);

        max_cdf_key = (int)diff_cdf.size()-1;
        printf("Finished cdf reading, first 10 elements of cdf:\n");
        for (int i = 0; i < 10; i++)
            printf("%.5e\n", diff_cdf[i]);
    }
    
    /* Read the list of key file paths */
    std::vector<std::string> key_files;
    
    f = fopen(list_keys_file, "r");
    if (f == NULL) {
        printf("Error opening file %s for reading\n", list_keys_file);
        return 1;
    }

    while (fgets(buf, 512, f)) {
        /* Remove trailing newline */
        if (buf[strlen(buf) - 1] == '\n')
            buf[strlen(buf) - 1] = 0;
        
        key_files.push_back(std::string(buf));
    }
    fclose(f);

    omp_set_num_threads(10);
    time_t start_time, finish_time;

    /* Initialize points */
    std::vector<Point> points;
    points.reserve(num_pts);
    for (int i = 0; i < num_pts; i++) {
        Point pt(i);
        points.push_back(pt);
    }
    printf("Finished initializing points.\n");
    fflush(stdout);

    /* Read image point matrix */
    printf("Start reading image point matrix...\n");

    f = fopen(img_pt_mat_file, "r");
    if (f == NULL) {
        printf("Error opening image point file %s\n", img_pt_mat_file);
        return -1;
    }

    unsigned char* desc;
    int img_id = 0, pt_id = 0, pos = 0;

    std::vector<TI> tripletIntList;
    SpColIntMat ip_mat(num_pts, num_img);

    time(&start_time);    
    while(fgets(buf, 512, f)) {
        /* Remove trailing newline */
        if (buf[strlen(buf) - 1] == '\n')
            buf[strlen(buf) - 1] = 0;

        sscanf(buf, "%d %d %d", &img_id, &pt_id, &pos);
        //printf("img id: %d, pt id: %d, pos: %d\n", img_id, pt_id, pos);
        if (!read_mean) {
            tripletIntList.push_back(TI(pt_id, img_id, (int) points[pt_id].m_imgs.size() + 1));
        }
        
        points[pt_id].m_imgs.push_back(img_id);
        points[pt_id].m_pos.push_back(pos);
        points[pt_id].m_probs.push_back(1.0f);

        if (!use_rm && !read_mean) {
            /* Have to read keys from scratch */
            desc = new unsigned char[128];
            for (int j = 0; j < 128; j++)
                desc[j] = 0;

            points[pt_id].m_desc.push_back(desc);
        }
    }

    if (!read_mean){
        ip_mat.setFromTriplets(tripletIntList.begin(),tripletIntList.end());
        tripletIntList.clear();
    }
    
    fclose(f);
    time(&finish_time);
    printf("Finished reading image point matrix, time used: %.3f seconds\n",
           difftime(finish_time, start_time));
    fflush(stdout);
    std::vector<float> diff_vec;
            
    if (!read_mean){
        /* Read the list of keys files directly */
        std::vector<bool> touched_pts; touched_pts.reserve(num_pts);
        for (int i = 0; i < num_pts; i++)
            touched_pts.push_back(false);

        clock_t start = clock();    
        for (int i = 0; i < (int) key_files.size(); i++) {
            unsigned char *keys_i;
            int num_i = ReadKeyFile(key_files[i].c_str(), &keys_i);
            if (num_i == 0){
                printf("Image %d has an empty key file!\n", i);
                continue;
            }

            assert(i < ip_mat.outerSize());
            for (SpColIntMat::InnerIterator it(ip_mat, i); it; ++it){
                int pt_id = it.row();
                int in_pt_pos = it.value() - 1;
                int in_img_pos = points[pt_id].m_pos[in_pt_pos];
                touched_pts[pt_id] = true;
                
                if (use_rm) {
                    points[pt_id].update_mean(keys_i + in_img_pos * 128 * sizeof(unsigned char));
                } else {
                    unsigned char* pt_desc = points[pt_id].m_desc[in_pt_pos];
                    for (int j = 0; j < 128; j++) 
                        pt_desc[j] = keys_i[in_img_pos*128 + j];
                }
            }
            
            clock_t end = clock();    
            if ((i+1) % 50 == 0){
                printf("Reading %d (%.2f %%) keys took %0.3fs\n", 
                       i+1, (float) i / key_files.size()*100, (end - start) / ((double) CLOCKS_PER_SEC));
                fflush(stdout);
            }
            delete[] keys_i;
        }
        printf("Finished reading key files.\n");
        fflush(stdout);

        if (use_rm){
            for (int i = 0; i < num_pts; i++){
                points[i].convert_mean();
            }
        } else {
            for (int i = 0; i < num_pts; i++){
                if (touched_pts[i])
                    points[i].compute_mean(diff_vec);
            }
            printf("Finished computing key means for points. Diff_vec length: %d.\n",
                   (int) diff_vec.size());
            fflush(stdout);    
        }

        printf("Finished converting mean to unsigned char.\n");
        fflush(stdout);
    } else {
        f = fopen(mean_path, "r");
        if (f == NULL) {
            printf("Error opening file %s for reading\n", mean_path);
            return 1;
        }

        for (int i = 0; i < num_pts; i++)
            points[i].init_mean();
        printf("Finished initializaing mean descriptors of points.\n");
        fflush(stdout);

        if (use_bin){
            int read_num_pts = 0;
            fread(&read_num_pts, sizeof(int), 1, f);
            assert(read_num_pts = num_pts);

            for (int i = 0; i < num_pts; i++){
                fread(points[i].m_mean_desc, sizeof(unsigned char), 128, f);
                if ((i+1) % 10000 == 0){
                    printf("Read mean descriptors for %d/%d points ... %.3f%%\n", 
                           (i+1), num_pts, (float)(i+1) * 100.0f/num_pts);
                }
            }
        } else {
            int read_pt_id = 0;
            while(fgets(buf, 2048, f)) {
                /* Remove trailing newline */
                if (buf[strlen(buf) - 1] == '\n')
                    buf[strlen(buf) - 1] = 0;

                int values[128];
                sscanf(buf, "%d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d",
                       &(values[0]), &(values[1]), &(values[2]), &(values[3]), &(values[4]), &(values[5]), &(values[6]), &(values[7]), &(values[8]), &(values[9]), &(values[10]), &(values[11]), &(values[12]), &(values[13]), &(values[14]), &(values[15]), &(values[16]), &(values[17]), &(values[18]), &(values[19]), &(values[20]), &(values[21]), &(values[22]), &(values[23]), &(values[24]), &(values[25]), &(values[26]), &(values[27]), &(values[28]), &(values[29]), &(values[30]), &(values[31]), &(values[32]), &(values[33]), &(values[34]), &(values[35]), &(values[36]), &(values[37]), &(values[38]), &(values[39]), &(values[40]), &(values[41]), &(values[42]), &(values[43]), &(values[44]), &(values[45]), &(values[46]), &(values[47]), &(values[48]), &(values[49]), &(values[50]), &(values[51]), &(values[52]), &(values[53]), &(values[54]), &(values[55]), &(values[56]), &(values[57]), &(values[58]), &(values[59]), &(values[60]), &(values[61]), &(values[62]), &(values[63]), &(values[64]), &(values[65]), &(values[66]), &(values[67]), &(values[68]), &(values[69]), &(values[70]), &(values[71]), &(values[72]), &(values[73]), &(values[74]), &(values[75]), &(values[76]), &(values[77]), &(values[78]), &(values[79]), &(values[80]), &(values[81]), &(values[82]), &(values[83]), &(values[84]), &(values[85]), &(values[86]), &(values[87]), &(values[88]), &(values[89]), &(values[90]), &(values[91]), &(values[92]), &(values[93]), &(values[94]), &(values[95]), &(values[96]), &(values[97]), &(values[98]), &(values[99]), &(values[100]), &(values[101]), &(values[102]), &(values[103]), &(values[104]), &(values[105]), &(values[106]), &(values[107]), &(values[108]), &(values[109]), &(values[110]), &(values[111]), &(values[112]), &(values[113]), &(values[114]), &(values[115]), &(values[116]), &(values[117]), &(values[118]), &(values[119]), &(values[120]), &(values[121]), &(values[122]), &(values[123]), &(values[124]), &(values[125]), &(values[126]), &(values[127]));
            
                for (int j = 0; j < 128; j++){
                    points[read_pt_id].m_mean_desc[j] = (unsigned char) values[j];
                }
            
                read_pt_id++;
                if (read_pt_id % 10000 == 0){
                    printf("Read mean descriptors for %d/%d points ... %.3f%%\n", 
                           read_pt_id, num_pts, (float)read_pt_id * 100.0f/num_pts);
                }
            }
        }
        
        fclose(f);
        //points[0].print_mean();
        //exit(1);
    }
    
    if (write_mean){
        FILE *mean_f = fopen(mean_path, "w");
        if (mean_f == NULL){
            printf("[reduce.cpp] Error opening file %s for writing.\n", mean_path);
            return -1;
        }
        if (use_bin)
            fwrite(&num_pts, sizeof(int), 1, mean_f);

        for (int i = 0; i < num_pts; i++){
            if (use_bin)
                points[i].fprint_bin_mean(mean_f);
            else
                points[i].fprint_mean(mean_f);
        }
        fclose(mean_f);

        printf("Finished writing %s\n", mean_path);
        fflush(stdout);
        exit(0);
    }

    if (!use_rm && record_dist) {
        // Compute difference value histograms
        printf("Difference values: \n");
        for (unsigned int i = 0; i < diff_vec.size(); i++) {
            if (diff_vec[i] == 0.0)
                continue;
            printf("%.5e ", diff_vec[i]);
        }
        printf("\n");
        fflush(stdout);
        return 0;

        printf("start allocating means_mat\n");
        fflush(stdout);

        flann::Matrix<float> means_mat(new float[128*num_pts], num_pts, 128);

        printf("finished allocating means_mat\n");
        fflush(stdout);

        // Compute min_dist for each point
        int tree_num = 8, nn = 2, check_num = 1024;
        int chunk_size = 100000; // 100K
        int chunk_num = num_pts / chunk_size;
        int chunk_rest = num_pts % chunk_size;
        int changed = 0;
        unsigned int *assignments = new unsigned int[num_pts];

        if (assignments == NULL){
            printf("Error allocating assignments!");
            exit(-1);
        }

        for (int i = 0; i < num_pts; i++){
            for (int j = 0; j < 128; j++) {
                means_mat[i][j] = points[i].m_mean_desc_db[j];
            }
        }

        flann::Index<flann::L2<float> > index(means_mat, flann::KDTreeIndexParams(tree_num));
        index.buildIndex();
        printf("Finished building random KD index with %d trees (%d checks)...\n", 
               tree_num, check_num); 
        fflush(stdout);

        std::vector<float> min_dist_vec;
        min_dist_vec.reserve(num_pts);
        for (int i = 0; i < chunk_num; i++){    
            knn_search_rnd_kd_forest(means_mat, i * chunk_size, chunk_size, 128, nn,
                                     index, check_num, min_dist_vec, changed, assignments);
        }

        if (chunk_rest > 0) {
            knn_search_rnd_kd_forest(means_mat, chunk_num * chunk_size, chunk_rest, 128, nn,
                                     index, check_num, min_dist_vec, changed, assignments);
        }

        printf("Min distance values:");
        for (int i = 0; i < (int)min_dist_vec.size(); i++) {
            printf(" %.5e", min_dist_vec[i]);
        }
        printf("\n");

        printf("Min distance indices:");
        for (int i = 0; i < num_pts; i++) {
            printf(" %d", assignments[i]);
        }
        printf("\n");

        printf("Changed: %d\n", changed); 
        fflush(stdout);
        exit(0);
    }

    //////////////////////////////////////////////////////////////////////
    // Stage 1: K-Cover Algorithm
    //////////////////////////////////////////////////////////////////////
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
        max_gains.push_back(0.0);
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
    SpRowMat covered_mat(num_img, num_pts); // current ip_mat
    SpRowMat prob_mat(num_img, num_pts); // current P(m) distribution
    std::vector<T> tripletList;

    // round statistics
    int num_img_covered = 0, num_pts_chosen = 0;
    int max_pt_id = -1;
    double max_score = 0.0;
    std::vector<int> selected_ids;

    // kd tree parameters
    std::vector<int> unindexed_ids;
    int chunk_size = 2000, chunk_num = 0;
    int tree_num = 8, nn = 1, check_num = 2048;
    flann::Matrix<float>* means_mat = NULL;
    flann::Index<flann::L2<float> >* index = NULL;

    printf("Finished initialization, starting first covering...\n");
    fflush(stdout);

    if (record_dist && use_load){
        FILE *pt_idx_f = fopen(pt_idx_path, "r");
        if (pt_idx_f == NULL) {
            printf("Error opening image point file %s\n", pt_idx_path);
            return -1;
        }

        int pt_idx = 0;
        while(fgets(buf, 512, pt_idx_f)) {
            if (buf[strlen(buf) - 1] == '\n')
                buf[strlen(buf) - 1] = 0;
            sscanf(buf, "%d", &pt_idx);
            add_point(pt_idx, kc, chosen_pts, selected_ids, num_pts_chosen,
                      covered_times, covered_mass, points, covered_imgs,
                      num_img_covered, num_img);
        }
        fclose(pt_idx_f);

        printf("Finished loading point index file %s.\n", pt_idx_path);

        char pt_idx_tmp[128];
        strcpy(pt_idx_tmp,pt_idx_path);
        
        FILE *pt_idx_nn_f = fopen(strcat(pt_idx_tmp,".nn"), "w");
        if (pt_idx_nn_f == NULL){
            printf("Error opening file %s for writing.\n", pt_idx_tmp);
            return -1;
        }

        std::vector<float> dist_vec;
        for (unsigned int i = 0; i < selected_ids.size(); i++){
            int pt_id = selected_ids[i];
            float min_pt_dist = points[pt_id].get_min_dist_to_points(points, selected_ids);
            dist_vec.push_back(min_pt_dist);
            fprintf(pt_idx_nn_f, "%d\n", points[pt_id].min_dist_id);
        }
        fclose(pt_idx_nn_f);

        printf("Distances to nearest neighbor: \n");
        for (unsigned int i = 0; i < dist_vec.size(); i++){
            printf("%.5e ", dist_vec[i]);
        }
        printf("\n");
        return 0;
    }
    
    // Start K-cover
    double time_dist = 0.0;
    time(&start_time);
    while (num_img_covered < num_img){
        max_score = 0.0;
        max_pt_id = -1;
        double time_dist_round = 0.0;
        clock_t start_round = clock();        
        for (int i = 0; i < num_pts; i++){
            /* pruning from submodularity */
            if (chosen_pts[i] || max_gains[i] <= max_score)
                continue;

            /* compute new cover */
            double new_cover = 0.0;
            for (unsigned int j = 0; j < points[i].m_imgs.size(); j++) {
                if (!covered_imgs[points[i].m_imgs[j]])
                    new_cover += points[i].m_probs[j];
            }

            if (new_cover <= max_score){
                /* update max_gain for pruning */
                max_gains[i] = new_cover;
                continue;
            }

            /* compute nearest distance to currently chosen set of points */
            float discount = 1.0f;
            bool first_round = selected_ids.empty();
            if (first_round) {
                discount = 1.0f;
                printf("First round, checking pt %d, skipping computing discount factor...\n", i);
            } else if (use_kc) {
                discount = 1.0f;
            } else if (use_direct) {
                discount = 1 - points[i].get_pct_closer_to(points, selected_ids);
                //printf("percentage closer to chosen points for point %d is %f\n", i, 1-discount);
                //fflush(stdout);
            } else if (!use_approx) {
                /*
                float min_dist = INFINITY;
                //#pragma omp parallel for
                for (unsigned int j = 0; j < selected_ids.size(); j++){
                    float new_dist = 0.0f;
                    if (use_cdf)
                        new_dist = points[selected_ids[j]].get_max_desc_dist(points[i]);
                    else
                        new_dist = points[selected_ids[j]].get_mean_dist(points[i]);
                    //#pragma omp critical
                    //{
                        if (new_dist < min_dist)
                            min_dist = new_dist;
                        //}
                }
                */
                clock_t start = clock();
                //printf("Evaluating point %d, max_score: %f, max_gain: %f --- ", 
                //       i, max_score, max_gains[i]);
                float min_dist = points[i].get_min_dist_to_points(points, selected_ids);
                clock_t end = clock();    
                time_dist_round += (end - start) / ((double) CLOCKS_PER_SEC);
                
                if (use_cdf) {
                    discount = get_cdf_weight(diff_cdf, max_cdf_key, min_dist);
                } else {
                    discount = get_thresh_weight(min_dist_cur, min_dist);
                }
                //printf("point %d, min_dist %f,  min_id %d, discount %f\n", i, min_dist, points[i].min_dist_id, discount);
                //printf("min_dist %f, min_id %d, discount %f, ", min_dist, points[i].min_dist_id, discount);
            } else {
                float min_dist = INFINITY;
                for (int j = 0; j < (int)unindexed_ids.size(); j++){
                    float new_dist = points[unindexed_ids[j]].get_mean_dist(points[i]);
                    if (new_dist < min_dist){
                        min_dist = new_dist;
                    }
                }

                if (index != NULL){
                    flann::Matrix<float> query(new float[128], 1, 128);
                    flann::Matrix<int> indices(new int[nn], 1, nn);
                    flann::Matrix<float> dists(new float[nn], 1, nn);
                    flann::SearchParams para(check_num);
                    para.cores = NUM_THREADS;

                    for (int j = 0; j < 128; j++) {
                        query[0][j] = points[i].m_mean_desc_db[j];
                    }

                    (*index).knnSearch(query, indices, dists, nn, para);
                    float new_dist = dists[0][0];

                    if (new_dist < min_dist){
                        min_dist = new_dist;
                    }

                    delete[] query.ptr();
                    delete[] dists.ptr();
                    delete[] indices.ptr();
                }

                if ((int) selected_ids.size() > chunk_size) {
                    printf("min_dist approx: %.5e, ", min_dist);

                    float min_dist_ver = INFINITY;
                    for (int j = 0; j < (int)selected_ids.size(); j++){
                        float new_dist = points[selected_ids[j]].get_mean_dist(points[i]);
                        if (new_dist < min_dist_ver){
                            min_dist_ver = new_dist;
                        }
                    }

                    printf("min_dist real: %.5e\n", min_dist_ver);
                }
            }
            
            /* compute score for the current point */
            double score = new_cover * discount;
            //printf("discount: %f, score: %f\n", discount, score);
            
            /* update max_gain for pruning */
            max_gains[i] = score;

            if (score > max_score) {
                max_score = score;
                max_pt_id = i;
            }
        }

        clock_t end_round = clock();    
        if(max_pt_id < 0) {
            printf("Adding more points won't cover new images, exiting.\n");
            break;
        }

        /* adding point max_pt_id */
        printf("Adding point %d in stage 1 with score %f, "
               "time used: %.3f, time used for computing distance: %.3f\n",
               max_pt_id, max_score, (end_round - start_round) / ((double) CLOCKS_PER_SEC), time_dist_round);
        time_dist += time_dist_round;

        /* adding point max_pt_id */
        add_point(max_pt_id, kc, chosen_pts, selected_ids, num_pts_chosen,
                  covered_times, covered_mass, points, covered_imgs,
                  num_img_covered, num_img);
        
        /* construct KD tree if necessary */
        if (use_approx){
            unindexed_ids.push_back(max_pt_id);
            if ((int)unindexed_ids.size() == chunk_size){
                if (means_mat != NULL){
                    delete[] (*means_mat).ptr();
                }

                if (index != NULL) {
                    delete index;
                }
                chunk_num++;
                int num_idx_points = (int) selected_ids.size();

                printf("Start building index, %d points selected, %d unindexed ...\n", 
                       num_idx_points, (int)unindexed_ids.size());
                fflush(stdout);

                unindexed_ids.clear();

                means_mat = new flann::Matrix<float> (new float[128 * num_idx_points], num_idx_points, 128);
                printf("Finished initializing means_mat.\n");
                fflush(stdout);            

                for (int i = 0; i < num_idx_points; i++){
                    for (int j = 0; j < 128; j++) {
                        (*means_mat)[i][j] = points[selected_ids[i]].m_mean_desc_db[j];
                    }
                }
                printf("Finished writing means_mat.\n");
                fflush(stdout);            

                index = new flann::Index<flann::L2<float> >((*means_mat), flann::KDTreeIndexParams(tree_num));
                printf("Finished initializing index.\n");
                fflush(stdout);            

                (*index).buildIndex();
                printf("Finished building random KD index for %d points with %d trees (%d checks)...\n", 
                       chunk_size*chunk_num, tree_num, check_num); 
                fflush(stdout);
            }
        }

        if (num_pts_chosen % 1000 == 0){
            printf("%d out of %d images covered (kc=%d, %f%%), "
                   "number of chosen points: %d (%f%%)\n",
                   num_img_covered, num_img, kc, (float)num_img_covered*100.0/num_img,
                   num_pts_chosen, (float)num_pts_chosen*100.0/num_pts);
        }

        if (use_rd) {
            reduce_redundancy(kc, chosen_pts, selected_ids, num_pts_chosen, covered_times,
                              covered_mass, points);
        }
    }
    time(&finish_time);
    printf("Finished 1st covering, time used: %.3f seconds, "
           "time used for computing distances: %.3f seconds, "
           "%d out of %d images covered (kc=%d), "
           "number of chosen points: %d (%f%%)\n",
           difftime(finish_time, start_time), time_dist, num_img_covered, num_img, kc,
           num_pts_chosen, (float)num_pts_chosen*100.0/num_pts);
    fflush(stdout);

    if (record_dist){
        std::vector<float> dist_vec;
        for (unsigned int i = 0; i < selected_ids.size(); i++){
            int pt_id = selected_ids[i];
            float min_pt_dist = points[pt_id].get_min_dist_to_points(points, selected_ids);
            dist_vec.push_back(min_pt_dist);
        }
        printf("Distances to nearest neighbor: \n");
        for (unsigned int i = 0; i < dist_vec.size(); i++){
            printf("%.5e ", dist_vec[i]);
        }
        printf("\n");
    }

    //////////////////////////////////////////////////////////////////////
    // Stage 2: adding points according to survival function (prob. of
    // an image seeing more than kc points)
    //////////////////////////////////////////////////////////////////////
    if (use_pb) {
        /* set covered_mat : probability weighted ip_mat */
        if (cdf_as_prob || thresh_as_prob){
            for (unsigned int i = 0; i < selected_ids.size(); i++){
                int pt_id = selected_ids[i];
                float pt_prob = 0.0;
                float min_pt_dist = points[pt_id].get_min_dist_to_points(points, selected_ids);
                
                if (cdf_as_prob) {
                    pt_prob = get_cdf_weight(diff_cdf, max_cdf_key, min_pt_dist);
                } else if (thresh_as_prob) {
                    pt_prob = get_thresh_weight(min_dist_cur, min_pt_dist);
                }
                
                for (unsigned int j = 0; j < points[pt_id].m_imgs.size(); j++){
                    int img_id_cur = points[pt_id].m_imgs[j];
                    points[pt_id].m_probs[j] = points[pt_id].m_probs[j] * pt_prob;
                    
                    tripletList.push_back(T(img_id_cur, pt_id, points[pt_id].m_probs[j] * ip_weight));
                }

                if ((i+1) % 500 == 0) {
                    printf("Evaluated prob matrix for %d points (%.3f%%)\n", i+1, (i+1) / (float)selected_ids.size()*100.0);
                }
            }
        } else {
            for (unsigned int i = 0; i < selected_ids.size(); i++){
                int pt_id = selected_ids[i];

                for (unsigned int j = 0; j < points[pt_id].m_imgs.size(); j++){
                    int img_id_cur = points[pt_id].m_imgs[j];

                    tripletList.push_back(T(img_id_cur, pt_id, points[pt_id].m_probs[j] * ip_weight));
                }
            }
        }
        covered_mat.setFromTriplets(tripletList.begin(),tripletList.end());
        tripletList.clear();

        //std::vector<double> max_pb_vals; max_pb_vals.reserve(num_img);
        num_img_covered = 0;
        for (int i = 0; i < num_img; i++) {
            double surv = getSurvFunc(kc, covered_times[i], 0, i,
                                      0.0f, covered_mat, prob_mat, 
                                      false, true);

            //double max_pb_val = 0.0;
            printf("Image %d surv %f distribution:", i, surv);
            
            for (SpRowMat::InnerIterator it (prob_mat, i); it; ++it){
                printf(" %f", it.value());
                //if (it.value() > max_pb_val)
                //max_pb_val = it.value();
            }
            printf("\n");
            //max_pb_vals.push_back(max_pb_val);
            //printf("image %d has max_pb_val %f\n", i, max_pb_val);

            covered_prob[i] = surv;
            covered_imgs[i] = (surv >= min_prob);

            if (covered_imgs[i]){
                num_img_covered++;
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

        /* use pb_gains to record utility value for each point */
        std::vector<double> pb_gains; pb_gains.reserve(num_pts);
        //std::vector<bool> dirty_vec; dirty_vec.reserve(num_pts);
        for (int i = 0; i < num_pts; i++) {
            pb_gains.push_back(0.0);
            /*
            dirty_vec.push_back(true);
            max_gains[i] = 0.0;
            for (unsigned int j = 0; j < points[i].m_imgs.size(); j++) {
                if (!covered_imgs[points[i].m_imgs[j]])
                    max_gains[i] += points[i].m_probs[j];
                    //max_gains[i] += points[i].m_probs[j] * max_pb_vals[points[i].m_imgs[j]];
            }
            */
        }
        
        printf("Finished resetting score array pb_gains...\n");
        fflush(stdout);
    
        time(&start_time);
        double tar_num_img = pct * num_img;
        float max_pt_prob = 0.0f;
        
        while (num_img_covered < tar_num_img){
            max_score = 0.0f;
            max_pt_id = -1;
            max_pt_prob = 0.0f;
            
#pragma omp parallel for
            for (int i = 0; i < num_pts; i++){
                //if (chosen_pts[i] || max_gains[i] <= max_score)
                //if (chosen_pts[i] || !dirty_vec[i] || max_gains[i] <= 0)
                if (chosen_pts[i])
                    continue;

                //printf("Evaluating point %d, max_score: %.3f, max_gains: %.3f, ", 
                //       i, max_score, max_gains[i]);
                //printf("Evaluating point %d, previous pb_gain: %.3f, ", 
                //       i, pb_gains[i]);

                pb_gains[i] = 0.0;
                //max_gains[i] = 0.0;

                for (unsigned int j = 0; j < points[i].m_imgs.size(); j++) {
                    int loc_img_id = points[i].m_imgs[j];

                    if (!covered_imgs[loc_img_id]){
                        pb_gains[i] += prob_mat.coeffRef(loc_img_id, kc-1) * points[i].m_probs[j];
                        //max_gains[i] += max_pb_vals[loc_img_id] * points[i].m_probs[j];
                    }
                }
                //printf("after pb_gain: %.3f\n", pb_gains[i]);
                //dirty_vec[i] = false;
                //printf("Skipping point %d, pb_gain: %.3f\n", 
                //       i, pb_gains[i]);
            }

            for (int i = 0; i < num_pts; i++){
                if (chosen_pts[i] || pb_gains[i] <= max_score)
                    continue;

                float pt_prob = 1.0f;                
                float min_dist = INFINITY;

                if (cdf_as_prob || thresh_as_prob) {
                    min_dist = points[i].get_min_dist_to_points(points, selected_ids);

                    if (use_cdf) {
                        pt_prob = get_cdf_weight(diff_cdf, max_cdf_key, min_dist);
                    } else {
                        pt_prob = get_thresh_weight(min_dist_cur, min_dist);
                    }
                }

                // compute score for the current point
                double score = pb_gains[i] * pt_prob;

                //printf("max_gains after %.3f, score %.3f\n", 
                //       max_score, max_gains[i]);

                if (score > max_score){
                    max_score = score;
                    max_pt_id = i;
                    max_pt_prob = pt_prob;
                }
            }

            if(max_pt_id < 0) {
                printf("Adding more points won't increase probabilities, exiting.\n");
                break;
            }

            /* adding point max_pt_id */
            printf("Adding point %d in stage 2 with score %f\n", max_pt_id, max_score);
            add_point_pb(max_pt_id, kc, min_prob, max_pt_prob, covered_prob,
                         prob_mat, covered_mat, chosen_pts, selected_ids, num_pts_chosen, 
                         covered_times, covered_mass, points,
                         covered_imgs, num_img_covered, num_img);
            
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
        printf("Skipping 2nd covering...\n");
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
    return 0;
}
