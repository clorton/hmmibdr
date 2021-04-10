
FILE* open_frequency_file(bool flag, const std::string& filename);
FILE* open_input_file(bool flag, const std::string& filename);
FILE* open_output_file(const std::string& filename);
bool load_bad_samples(const std::string& filename, std::vector<std::string>& samples);
std::vector<std::string> split(const std::string& string, const char delimiter, bool collapse=true);
bool load_good_samples(const std::string& filename, std::vector<std::pair<std::string, std::string>>& samples);
void read_text_file(const std::string& filename, std::vector<std::string>& contents);

FILE* open_frequency_file(bool flag, const std::string& filename)
{
    FILE* file = nullptr;

    if (flag) {
        file = fopen(filename.c_str(), "r");
        if (!file) {
            REprintf("Could not open frequency file %s\n", filename.c_str());
            Rcpp::stop("");
        }
    }
  
  return file;
}

FILE* open_input_file(bool flag, const std::string& filename)
{
  FILE* file = nullptr;

  if (flag) {
      file = fopen(filename.c_str(), "r");
      if (!file) {
          REprintf("Could not open input file %s\n", filename.c_str());
          Rcpp::stop("");
      }
  }

  return file;
}

FILE* open_output_file(const std::string& filename)
{
    FILE* file = fopen(filename.c_str(), "w");

    if (!file) {
        REprintf("Could not open output file %s\n", filename.c_str());
        Rcpp::stop("");
    }

  return file;
}

bool load_bad_samples(const std::string& filename, std::vector<std::string>& samples)
{
    bool wasSuccessful = true;

    std::ifstream input;
    input.open(filename);

    if (!input.is_open()) {
        REprintf("Could not open file of bad samples: %s\n", filename.c_str());
        wasSuccessful = false;
    }
    else {
        std::string line;
        while (std::getline(input, line)) {
            samples.push_back(line);
        }
        input.close();
    }

    return wasSuccessful;
}

bool load_good_samples(const std::string& filename, std::vector<std::pair<std::string, std::string>>& samples)
{
    bool wasSuccessful = true;

    std::ifstream input;
    input.open(filename);

    if (!input.is_open()) {
        REprintf("Could not open file of good sample pairs: %s\n", filename.c_str());
        wasSuccessful = false;
    }
    else {
        std::string line;
        while (std::getline(input, line)) {
            std::vector<std::string> tokens = split(line, '\t');
            switch (tokens.size()) {
                case 1:
                  REprintf("Skipped line with only one record in file of good sample pairs: %s.\n", line.c_str());
                  break;

                case 2:
                  samples.push_back(std::pair<std::string, std::string>(tokens[0], tokens[1]));
                  break;

                default:
                  REprintf("Encountered bad line in %s: '%s'\n", filename.c_str(), line.c_str());
                  break;
            }
        }

        input.close();
    }

    return wasSuccessful;
}

std::vector<std::string> split(const std::string& string, const char delimiter, bool collapse)
{
    std::vector<std::string> splits;

    size_t start = 0;
    size_t stop = std::string::npos;
    do {
        stop = string.find(delimiter, start);
        // std::cout << "start = " << start << ", stop = " << stop << std::endl;
        if ( (stop-start > 0) || !collapse) {
            splits.push_back(string.substr(start, stop-start));
        }
        start = stop != std::string::npos ? stop + 1 : std::string::npos;
    } while (start != std::string::npos);

    return splits;
}

void read_text_file(const std::string& filename, std::vector<std::string>& contents)
{
    std::ifstream file;
    file.open(filename);

    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            contents.push_back(line);
        }
        file.close();
    }
}

int store_sample_names(const std::string& header,
                       std::vector<std::string>& samples,
                       const std::vector<std::string>& filter_list,
                       std::vector<bool>& use_sample)
{
    std::vector<std::string> tokens = split(header, '\t');
    int count = 0;

    for ( int isamp = 0, i = 2; i < tokens.size(); ++i ) {

        samples.push_back(tokens[i]);

        for ( const auto& filtered : filter_list ) {
            if ( tokens[i] == filtered ) {
                use_sample[isamp] = false;
                Rprintf("killing sample %s\n", tokens[i].c_str());
                break;
            }
        }

        if ( use_sample[isamp] ) {
            ++count;
        }

        ++isamp;
    }

    return count;
}

void display_parameters(int niter,
                        int min_snp_sep,
                        int min_inform,
                        float min_discord,
                        float max_discord,
                        float eps,
                        bool iflag2,
                        const std::string& data_file1,
                        const std::string& data_file2,
                        bool freq_flag1,
                        const std::string& freq_file1,
                        bool freq_flag2,
                        const std::string& freq_file2,
                        bool nflag,
                        float k_rec_max,
                        int nsample_use1,
                        int nsample_use2,
                        int nsample1,
                        int nsample2,
                        int& npair,
                        int& npair_report)
{
    Rprintf("Maximum fit iterations allowed: %d\n", niter);
    Rprintf("Minimum marker spacing (bp): %d\n", min_snp_sep);
    Rprintf("Minimum informative markers: %d\n", min_inform);
    Rprintf("Pairs accepted with discordance in range (%.2f%%, %.2f%%)\n",
            min_discord*100, max_discord*100);
    Rprintf("Genotyping error rate: %.2f%%\n", eps*100);
    if (iflag2) {
        Rprintf("Input files: %s and %s\n", data_file1.c_str(), data_file2.c_str());
    }
    else {
        Rprintf("Input file: %s\n", data_file1.c_str());
    }
    Rprintf("Frequency files: ");
    if (freq_flag1) {Rprintf("%s and ", freq_file1.c_str());}
    else {Rprintf("none and ");}
    if (freq_flag2) {Rprintf("%s\n", freq_file2.c_str());}
    else {Rprintf("none\n");}
    if (nflag) {Rprintf("Number of generations capped at %.2f\n", k_rec_max);}

    npair = nsample_use1 * nsample_use2;
    npair_report = nsample_use1 * (nsample_use1-1) / 2;
    Rprintf("pop1 nsample: %d used: %d", nsample1, nsample_use1);
    if (iflag2) {
        npair_report = npair;
        Rprintf("  pop2 nsample: %d used: %d", nsample2, nsample_use2);}
    Rprintf(" Expected pairs: %d\n", npair_report);
}

void determine_pairs_to_compare(int nsample1, bool iflag2, int nsample2, bool gflag, bool** use_pair,
                                const std::vector<std::pair<std::string, std::string>>& good_pairs,
                                const std::vector<std::string> sample1,
                                const std::vector<std::string> sample2)
{
    // All pairs are to be used by default, unless we've read a file of good pairs
    // Except single-pop samples should not be compared with themselves
    for (int isamp = 0; isamp < nsample1; isamp++) {

        // If we are comparing two distinct populations, consider all samples in the other population;
        // otherwise just consider the remaining samples of this population.
        int jstart = (iflag2 ? 0 : isamp+1);

        for (int jsamp = jstart; jsamp < nsample2; jsamp++) {

            // If good pairs were explicitly provided...
            if (gflag) {
                use_pair[isamp][jsamp] = false;
                for ( int ipair = 0; ipair < good_pairs.size(); ++ipair ) {
                    if ( ((good_pairs[ipair].first == sample1[isamp]) && (good_pairs[ipair].second == sample2[jsamp])) ||
                         ((good_pairs[ipair].first == sample2[jsamp]) && (good_pairs[ipair].second == sample1[isamp])) ) {
                        use_pair[isamp][jsamp] = true;
                    }
                }
            }
            // No explicit list of good pairs, don't compare with self unless using two distinct populations.
            else {
                if (isamp != jsamp || iflag2) {
                    use_pair[isamp][jsamp] = true;
                }
            }
        }
    }
}
