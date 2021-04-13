
#include <functional>

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

void increase_storage(
    int*& nall,
    int& max_snp,
    int*& pos,
    int nsample1,
    int** geno1,
    double**& freq1,
    int max_all,
    bool iflag2,
    int nsample2,
    int** geno2,
    double**& freq2
)
{
    nall = (int *)realloc(nall, 2*max_snp*sizeof(int));
    pos = (int *)realloc(pos, 2*max_snp*sizeof(int));

    for (int isamp = 0; isamp < nsample1; ++isamp) {
        geno1[isamp] = (int *)realloc(geno1[isamp], 2*max_snp*sizeof(int));
    }

    freq1 = (double **)realloc(freq1, 2*max_snp*sizeof(double*));

    for (int isnp = max_snp; isnp < 2*max_snp; ++isnp) {
        freq1[isnp] = new double[max_all+1];
    }

    if (iflag2) {

        for (int isamp = 0; isamp < nsample2; ++isamp) {
            geno2[isamp] = (int *)realloc(geno2[isamp], 2*max_snp*sizeof(int));
        }

        freq2 = (double **)realloc(freq2, 2*max_snp*sizeof(double*));

        for (int isnp = max_snp; isnp < 2*max_snp; ++isnp) {
            freq2[isnp] = new double[max_all+1];
        }
    }   // end if iflag2
    else {
        geno2 = geno1;
        freq2 = freq1;
    }

    max_snp *= 2;
}

void parse_input1_line(
    const std::string& input_line,
    int&  chromosome,
    int*  position,
    int   nsnp,
    int   previous_chromosome,
    int   min_snp_sep,
    int&  nskipped,
    bool& killit,
    int&  allele,
    int   maximum_alleles,
    int&  excessive_alleles,
    int** genome,
    const std::vector<bool>& use_sample,
    int*  allele_count,
    int&  total_alleles
)
{
    std::vector<std::string> tokens = split(input_line, '\t');
    for ( int itoken = 0; itoken < tokens.size(); ++itoken ) {

        const std::string& token = tokens[itoken];

        switch (itoken) {
            case 0:
                {
                    try {
                        chromosome = std::stol(token);
                    } catch (...) {
                        REprintf("Invalid chromosome '%s' (must be integer)\n", token.c_str());
                        Rcpp::stop("");
                    }
                }
                break;

            case 1:
                {
                    try {
                        position[nsnp] = std::stol(token);
                    } catch (...) {
                        REprintf("Invalid position '%s' (must be integer)\n", token.c_str());
                        Rcpp::stop("");
                    }

                    if ((chromosome == previous_chromosome && position[nsnp] < position[nsnp - 1]) || chromosome < previous_chromosome) {
                        REprintf("Variants are out of order\n");
                        Rcpp::stop("");
                    }

                    if (nsnp > 0 && chromosome == previous_chromosome && position[nsnp] - position[nsnp - 1] < min_snp_sep) {
                        ++nskipped;
                        killit = true;
                    }
                }
                break;

            default:
                {
                    try {
                        allele = std::stol(token);
                    } catch (...) {
                        REprintf("Invalid allele '%s' (must be integer)\n", token.c_str());
                        Rcpp::stop("");
                    }

                    if (allele > maximum_alleles) {
                        ++excessive_alleles;
                        killit = true;
                    }
                    else {
                        genome[itoken - 2][nsnp] = allele;

                        if (use_sample[itoken - 2]) {
                            if (allele >= 0) {
                                allele_count[allele]++;
                                ++total_alleles;
                            }
                        }
                    }
                }
                break;
        }

        if (killit) break;
    }
}

void parse_input2_line(
    const std::string& input_line,
    int& chromosome,
    int& position,
    int* positions,
    int nsnp,
    int pop1_chromosome,
    int& allele,
    int maximum_alleles,
    bool& killit,
    int& excessive_alleles,
    int** genome,
    const std::vector<bool>& use_sample,
    int* allele_count,
    int& total_alleles
)
{
    std::vector<std::string> tokens = split(input_line, '\t');
    for ( int itoken = 0; itoken < tokens.size(); ++itoken ) {

        const std::string& token = tokens[itoken];

        switch (itoken) {
            case 0:
                try {
                    chromosome = std::stol(token);
                } catch (...) {
                    REprintf("Invalid chromosome '%s' (must be integer)\n", token.c_str());
                    Rcpp::stop("");
                }
                break;

            case 1:
                try {
                    position = std::stol(token);
                } catch (...) {
                    REprintf("Invalid position '%s' (must be integer)\n", token.c_str());
                    Rcpp::stop("");
                }

                if (position != positions[nsnp] || chromosome != pop1_chromosome) {
                    REprintf("Data files do not agree on SNPs\n");
                    Rcpp::stop("");
                }
                break;

            default:
                try {
                    allele = std::stol(token);
                } catch (...) {
                    REprintf("Invalid allele '%s' (must be integer)\n", token.c_str());
                    Rcpp::stop("");
                }

                if (allele > maximum_alleles) {
                    ++excessive_alleles;
                    killit = true;
                }
                else {
                    genome[itoken - 2][nsnp] = allele;

                    if (use_sample[itoken - 2]) {
                        if (allele >= 0) {
                            allele_count[allele]++;
                            ++total_alleles;
                        }
                    }
                }
                break;
        }

        if (killit) break;
    }
}

void parse_frequency_line(
    const std::string& input_line,
    double* frequencies,
    int max_all,
    int& fposition,
    int& fchromosome,
    int chromosome,
    int* positions,
    int nsnp
)
{
    // Clear previous frequencies (since might have skipped previous snp via 'continue')
    memset(frequencies, 0, (max_all + 1) * sizeof(double));

    fposition = fchromosome = 0;

    std::vector<std::string> tokens = split(input_line, '\t');
    for ( int itoken = 0; itoken < tokens.size(); ++itoken ) {

        const std::string& token = tokens[itoken];

        switch (itoken) {
            case 0:
                try {
                    fchromosome = std::stol(token);
                } catch (...) {
                    REprintf("Invalid chromosome '%s' (must be integer)\n", token.c_str());
                    Rcpp::stop("");
                }
                break;

            case 1:
                try {
                    fposition = std::stol(token);
                } catch (...) {
                    REprintf("Invalid position '%s' (must be integer)\n", token.c_str());
                    Rcpp::stop("");
                }
                break;

            default:
                try {
                    frequencies[itoken - 2] = std::stod(token);
                } catch (...) {
                    REprintf("Invalid frequency '%s' (must be numeric)\n", token.c_str());
                    Rcpp::stop("");
                }
                break;
        }
    }

    if (fchromosome != chromosome || fposition != positions[nsnp]) {
        REprintf("Mismatch between data file and frequency file. Data file (chr/pos): %d/%d vs freq file: %d/%d\n",
                 chromosome, positions[nsnp], fchromosome, fposition);
        Rcpp::stop("");
    }
}

void calculate_allele_frequencies(
    int maximum_alleles,
    bool freq_flag1,
    double** freq1,
    int nsnp,
    double* ffreq1,
    int* pop1_allele_counts,
    int pop1_total_alleles,
    bool iflag2,
    bool freq_flag2,
    double** freq2,
    double* ffreq2,
    int* pop2_allele_counts,
    int pop2_total_alleles,
    double& fmean,
    double& maximum_frequency,
    int& majall,
    int* alleles_per_snp
)
{
    // process this variant -- calculate allele frequencies for each pop
    for (int allele_index = 0; allele_index <= maximum_alleles; ++allele_index) {

        if (freq_flag1) { freq1[nsnp][allele_index] = ffreq1[allele_index]; }
        else {            freq1[nsnp][allele_index] = (double) pop1_allele_counts[allele_index] / pop1_total_alleles; }

        if (iflag2) {
            if (freq_flag2) { freq2[nsnp][allele_index] = ffreq2[allele_index]; }
            else {            freq2[nsnp][allele_index] = (double) pop2_allele_counts[allele_index] / pop2_total_alleles; }
        }

        fmean = (freq1[nsnp][allele_index] + freq2[nsnp][allele_index]) / 2;

        if (fmean > maximum_frequency) {
            maximum_frequency = fmean;
            majall = allele_index;
        }

        if (fmean > 0) {
            alleles_per_snp[nsnp]++;
        }
    }
}

void tabulate_differences_by_pair_for_discordance(
    int pop1_nsamples,
    const std::vector<bool>& pop1_use_sample,
    bool pop2_present,
    int pop2_nsamples,
    const std::vector<bool>& pop2_use_sample,
    int** pop1_genome,
    int** pop2_genome,
    int majall,
    int* same_min,
    int* diff,
    int* nmiss_bypair,
    int nsnp,
    int chromosome,
    std::vector<int>& start_chromosome,
    std::vector<int>& end_chromosome
)
// Tabulate differences by pair for calculating discordance
{
    int ipair = 0;
    for (int isamp = 0; isamp < pop1_nsamples; ++isamp) {

        if (!pop1_use_sample[isamp]) {continue;}

        // If 2 pops, need to loop over all combinations, but not if one pop
        int jstart = (pop2_present ? 0 : isamp + 1);
        for (int jsamp = jstart; jsamp < pop2_nsamples; ++jsamp) {

            if (!pop2_use_sample[jsamp]) {continue;}

            if (pop1_genome[isamp][nsnp] != -1 && pop2_genome[jsamp][nsnp] != -1) {
                if (pop1_genome[isamp][nsnp] == pop2_genome[jsamp][nsnp]) {
                    if (pop1_genome[isamp][nsnp] != majall) {
                        same_min[ipair]++;
                    }
                }
                else {diff[ipair]++;}
            }
            else {
                nmiss_bypair[ipair]++;
            }

            ++ipair;
        }
    }

    if (nsnp < start_chromosome[chromosome]) {
        start_chromosome[chromosome] = nsnp;
    }

    if (nsnp > end_chromosome[chromosome]) {
        end_chromosome[chromosome] = nsnp;
    }
}

void filter_pairs_by_discordance_and_markers(
    int nsample1,
    const std::vector<bool>& use_sample1,
    bool iflag2,
    int nsample2,
    const std::vector<bool>& use_sample2,
    int& sum,
    int* diff,
    int* same_min,
    double** discord,
    bool** use_pair,
    double max_discord,
    int min_inform,
    double min_discord,
    int& nuse_pair
)
{
    int ipair = 0;
    for (int isamp = 0; isamp < nsample1; ++isamp) {

        if (!use_sample1[isamp]) {continue;}

        int jstart = (iflag2 ? 0 : isamp+1);
        for (int jsamp = jstart; jsamp < nsample2; ++jsamp) {

            if (!use_sample2[jsamp]) {continue;}

            sum = diff[ipair] + same_min[ipair];

            if (sum == 0) {
                discord[isamp][jsamp]  = 0;
            }
            else {
                discord[isamp][jsamp] = (double) diff[ipair] / sum;
            }

            if (use_pair[isamp][jsamp]) {
                if ( (discord[isamp][jsamp] >  max_discord) ||
                      (sum < min_inform) ||
                      (discord[isamp][jsamp] < min_discord) ) {

                    use_pair[isamp][jsamp] = false;
                }
                else {
                    ++nuse_pair;
                }
            }

            ++ipair;
        }
    }

    Rprintf("sample pairs analyzed (filtered for discordance and informative markers): %d\n", nuse_pair);
}

typedef std::function<void(const std::string&, const std::string&, int, int)> output1_t;
typedef std::function<void(int, int, int)> output2_t;
typedef std::function<void(const std::string&, const std::string&, int, double, double, int, double, int, double, double, double)> output3_t;

void calculate_hmm_ibd(
    int maxlen,
    int nsample1,
    const std::vector<bool>& use_sample1,
    bool iflag2,
    int nsample2,
    const std::vector<bool>& use_sample2,
    int* diff,
    int* same_min,
    bool** use_pair,
    double pi[2],
    double pinit[2],
    double k_rec_init,
    int niter,
    int nchrom,
    const std::vector<int>& start_chr,
    const std::vector<int>& end_chr,
    int* pos,
    int** geno1,
    int** geno2,
    double eps,
    int* nall,
    double** freq1,
    double** freq2,
    double rec_rate,
    const std::vector<std::string>& sample1,
    const std::vector<std::string>& sample2,
    output1_t output_fn1,
    output2_t output_fn2,
    output3_t freq_fn,
    bool nflag,
    double k_rec_max,
    double fit_thresh_dpi,
    double fit_thresh_dk,
    double fit_thresh_drelk,
    double** discord
)
{
    int* psi[2];
    double* phi[2];
    double* b[2];
    double* alpha[2];
    double* beta[2];

    for (int is = 0; is < 2; ++is) {
        psi[is] = new int[maxlen];
        phi[is] = new double[maxlen];
        b[is] = new double[maxlen];
        alpha[is] = new double[maxlen];
        beta[is] = new double[maxlen];
    }

    double* scale = new double[maxlen];
    int* traj = new int[maxlen];

    int ipair = 0;
    for (int isamp = 0; isamp < nsample1; ++isamp) {

        if (!use_sample1[isamp]) {continue;}

        int jstart = (iflag2 ? 0 : isamp+1);
        for (int jsamp = jstart; jsamp < nsample2; ++jsamp) {

            if (!use_sample2[jsamp]) {continue;}

            int sum = diff[ipair] + same_min[ipair];

            double seq_ibd = 0;
            double seq_dbd = 0;
            double count_ibd_fb = 0;
            double count_dbd_fb = 0;
            double seq_ibd_fb = 0;
            double seq_dbd_fb = 0;
            int count_ibd_vit = 0;
            int count_dbd_vit = 0;
            double max_phi = 0;

            if (use_pair[isamp][jsamp]) {

                double last_prob = 0;
                double last_pi = pi[0] = pinit[0];  // initialize with prior
                pi[1] = pinit[1];
                double k_rec = k_rec_init;
                double last_krec = k_rec;
                bool finish_fit = false;

                double fmean, fmeani, fmeanj;

                // Declare these out here because we will output them along with the other data after the loop.
                int ntrans = 0;
                int iter;
                for (iter = 0; iter < niter; ++iter) {

                    double trans_obs = 0;
                    double trans_pred = 0;

                    for (int chr = 1; chr <= nchrom; ++chr) {

                        int nsite = 0;

                        if (end_chr[chr] < 0) {continue;}

                        int snp_ind = 0;
                        for (int isnp = start_chr[chr]; isnp <= end_chr[chr]; ++isnp) {

                            snp_ind = isnp - start_chr[chr];
                            int gi = geno1[isamp][isnp];
                            int gj = geno2[jsamp][isnp];
                            double pright = 1 - eps * (nall[isnp] - 1);

                            if (gi == -1 || gj == -1) {
                                // O = n (missing data)
                                b[0][snp_ind] = b[1][snp_ind] = 1.0;
                            }
                            else if (gi == gj) {
                                // homozygote
                                ++nsite;
                                fmean = (freq1[isnp][gi] + freq2[isnp][gi]) / 2;
                                b[0][snp_ind] = pright * pright * fmean + eps * eps * (1 - fmean); //IBD
                                b[1][snp_ind] = pright * pright * freq1[isnp][gi]*freq2[isnp][gj] +
                                                pright * eps * freq1[isnp][gi] * (1-freq2[isnp][gj]) +
                                                pright * eps * (1 - freq1[isnp][gi]) * freq2[isnp][gj] +
                                                eps * eps * (1-freq1[isnp][gi]) * (1-freq2[isnp][gj]);
                            }
                            else {
                                // O = aA
                                ++nsite;
                                fmeani = (freq1[isnp][gi] + freq2[isnp][gi]) / 2;
                                fmeanj = (freq1[isnp][gj] + freq2[isnp][gj]) / 2;
                                b[0][snp_ind] = pright*eps*(fmeani+fmeanj) +
                                                eps * eps * (1 - fmeani - fmeanj);
                                b[1][snp_ind] = pright * pright * freq1[isnp][gi] * freq2[isnp][gj] +
                                                pright*eps*( freq1[isnp][gi]*(1-freq2[isnp][gj]) + freq2[isnp][gj]*(1-freq1[isnp][gi]) ) +
                                                eps * eps * (1 - freq1[isnp][gi]) * (1 - freq2[isnp][gj]);
                            }

                            if (isnp == start_chr[chr]) {

                                psi[0][snp_ind] = psi[1][snp_ind] = 0;

                                for (int is = 0; is < 2; ++is) {
                                    phi[is][snp_ind] = log(pi[is]) + log(b[is][snp_ind]);
                                    alpha[is][snp_ind] = pi[is] * b[is][snp_ind];
                                    scale[is] = 1.;
                                }
                            }
                            else {
                                double ptrans = k_rec * rec_rate * (pos[isnp] - pos[isnp-1]);
                                double a[2][2];
                                a[0][1] = 1 - pi[0] - (1 - pi[0]) * exp(-ptrans);
                                a[1][0] = 1 - pi[1] - (1 - pi[1]) * exp(-ptrans);
                                a[0][0] = 1 - a[0][1];
                                a[1][1] = 1 - a[1][0];

                                for (int js = 0; js < 2; ++js) {    // index over state of current snp

                                    double maxval = -10000000;
                                    alpha[js][snp_ind] = scale[snp_ind] = 0;

                                    for (int is = 0; is < 2; ++is) {    // index over state of previous snp

                                        if (phi[is][snp_ind-1] + log(a[is][js]) > maxval ) {
                                            maxval = phi[is][snp_ind-1] + log(a[is][js]);
                                            psi[js][snp_ind] = is;
                                        }

                                        phi[js][snp_ind] = maxval + log(b[js][snp_ind]);
                                        alpha[js][snp_ind] += alpha[is][snp_ind-1] * a[is][js] * b[js][snp_ind];
                                    }

                                    scale[snp_ind] += alpha[js][snp_ind];
                                }

                                for (int js = 0; js < 2; ++js) {    // scale alpha to prevent underflow (Rabiner eqns 92)
                                    alpha[js][snp_ind] /= scale[snp_ind];
                                }
                            }   // end if initializing/continuing
                        }  // end snp loop

                        int last_snp = snp_ind;
                        double max_phiL = phi[1][snp_ind];

                        if (phi[0][snp_ind] > phi[1][snp_ind]) {max_phiL = phi[0][snp_ind];}

                        int max = (phi[1][snp_ind] > phi[0][snp_ind]) ? 1 : 0;
                        traj[snp_ind] = max;
                        max_phi += max_phiL;

                        // Loop backward to calculate betas (backward part of forward-backward)
                        snp_ind = end_chr[chr] - start_chr[chr];
                        beta[0][snp_ind] = beta[1][snp_ind] = 1;

                        for (int isnp = end_chr[chr]-1; isnp >= start_chr[chr]; isnp--) {

                            int snp_ind = isnp - start_chr[chr];
                            double ptrans = k_rec * rec_rate * (pos[isnp+1] - pos[isnp]);
                            double a[2][2];
                            a[0][1] = 1 - pi[0] - (1 - pi[0]) * exp(-ptrans);
                            a[1][0] = 1 - pi[1] - (1 - pi[1]) * exp(-ptrans);
                            a[0][0] = 1 - a[0][1];
                            a[1][1] = 1 - a[1][0];

                            for (int is = 0; is < 2; ++is) {    // index over state of current snp

                                beta[is][snp_ind] = 0;

                                for (int js = 0; js < 2; ++js) {
                                    // index over state of previous snp (= snp+1, since looping backward)
                                    beta[is][snp_ind] += beta[js][snp_ind+1] * a[is][js] * b[js][snp_ind+1] / scale[snp_ind];
                                }
                            }
                        }

                        // inverse loop for Viterbi
                        for (int snp_ind = last_snp-1; snp_ind >= 0; snp_ind--) {
                            traj[snp_ind] = psi[max][snp_ind+1];
                            max = traj[snp_ind];
                        }

                        // tabulate FB results
                        for (int snp_ind = 0; snp_ind <= last_snp; ++snp_ind) {

                            double p_ibd = alpha[0][snp_ind] * beta[0][snp_ind] /
                                           (alpha[0][snp_ind] * beta[0][snp_ind] + alpha[1][snp_ind] * beta[1][snp_ind]);

                            count_ibd_fb += p_ibd;
                            count_dbd_fb += (1-p_ibd);

                            if (snp_ind < last_snp) {

                                int isnp = snp_ind + start_chr[chr];
                                int delpos = pos[isnp+1] - pos[isnp];
                                double ptrans = k_rec * rec_rate * delpos;
                                double a[2][2];
                                a[0][1] = 1 - pi[0] - (1 - pi[0]) * exp(-ptrans);
                                a[1][0] = 1 - pi[1] - (1 - pi[1]) * exp(-ptrans);
                                a[0][0] = 1 - a[0][1];
                                a[1][1] = 1 - a[1][0];
                                double xi[2][2];
                                double xisum = 0;

                                for (int is = 0; is < 2; ++is) {
                                    for (int js = 0; js < 2; ++js) {
                                        xi[is][js] = alpha[is][snp_ind] * a[is][js] * b[js][snp_ind+1] * beta[js][snp_ind+1];
                                        xisum += xi[is][js];
                                    }
                                }

                                double gamma[2] = { 0, 0 };
                                for (int is = 0; is < 2; ++is) {
                                    for (int js = 0; js < 2; ++js) {
                                        xi[is][js] /= xisum;
                                        gamma[is] += xi[is][js];
                                    }
                                }

                                // xi(0,0) in Rabiner notation = prob of being in state 0 at t and 0 at t+1
                                seq_ibd_fb += delpos * xi[0][0];
                                seq_dbd_fb += delpos * xi[1][1];
                                trans_obs += (xi[0][1] + xi[1][0]);
                                trans_pred += gamma[0] * a[0][1] + gamma[1] * a[1][0];
                            }
                        }

                        // print final Viterbi trajectory
                        // start
                        if (iter == niter - 1 || finish_fit) {
                            if (nsite > 0) {

                                // fprintf(outf, "%s\t%s\t%d\t%d", sample1[isamp].c_str(), sample2[jsamp].c_str(), chr, pos[0+start_chr[chr]]);
                                output_fn1( sample1[isamp], sample2[jsamp], chr, pos[0+start_chr[chr]]);

                                int start_snp = 0;
                                for (int isnp = 1; isnp < end_chr[chr] - start_chr[chr] + 1; ++isnp) {

                                    if (traj[isnp] != traj[isnp-1]) {

                                        ++ntrans;
                                        int add_seq = pos[isnp - 1 + start_chr[chr]] - pos[start_snp + start_chr[chr]] + 1;

                                        if (traj[isnp-1] == 0) { seq_ibd += add_seq; }
                                        else                   { seq_dbd += add_seq; }

                                        // end one and start another
                                        // fprintf(outf, "\t%d\t%d\t%d\n", pos[isnp - 1 + start_chr[chr]], traj[isnp-1], isnp - start_snp);
                                        // fprintf(outf, "%s\t%s\t%d\t%d", sample1[isamp].c_str(), sample2[jsamp].c_str(), chr, pos[isnp + start_chr[chr]]);

                                        output_fn2(pos[isnp - 1 + start_chr[chr]], traj[isnp-1], isnp - start_snp);
                                        output_fn1(sample1[isamp], sample2[jsamp], chr, pos[isnp + start_chr[chr]]);

                                        start_snp = isnp;
                                    }
                                }

                                int isnp = end_chr[chr] - start_chr[chr];
                                int add_seq = pos[isnp + start_chr[chr]] - pos[start_snp + start_chr[chr]] + 1;

                                if (traj[isnp] == 0) { seq_ibd += add_seq; }
                                else                 { seq_dbd += add_seq; }

                                // fprintf(outf, "\t%d\t%d\t%d\n", pos[end_chr[chr]], traj[isnp], isnp - start_snp + 1);
                                output_fn2(pos[end_chr[chr]], traj[isnp], isnp - start_snp + 1);

                                // Tabulate sites by state for Viterbi trajectory
                                for (int isnp = 0; isnp < end_chr[chr] - start_chr[chr] + 1; ++isnp) {
                                    if (traj[isnp] == 0) { ++count_ibd_vit; }
                                    else                 { ++count_dbd_vit; }
                                }
                            }
                        }
                    }  // end chrom loop

                    // quit if fit converged on previous iteration; otherwise, update parameters
                    if (finish_fit) {break;}

                    if ((count_ibd_fb + count_dbd_fb) == 0) {
                        REprintf("Insufficient information to estimate parameters.\n");
                        REprintf("Do you have only one variant per chromosome?\n");
                        Rcpp::stop("");
                    }

                    pi[0] = count_ibd_fb / (count_ibd_fb + count_dbd_fb);
                    k_rec *= trans_obs / trans_pred;

                    if (nflag && k_rec > k_rec_max) {k_rec = k_rec_max;}

                    // ad hoc attempt to avoid being trapped in extremum
                    if (iter < niter-1 && !finish_fit) {

                        if (pi[0] < 1e-5) {pi[0] = 1e-5;}
                        else if (pi[0] > 1- 1e-5) {pi[0] = 1 - 1e-5;}

                        if (k_rec < 1e-5) {k_rec = 1e-5;}
                    }

                    pi[1] = 1 - pi[0];
                    double delpi = pi[0] - last_pi;
                    double delk = k_rec - last_krec;

                    if (nflag && k_rec > k_rec_max) {delk = k_rec_max - last_krec;}

                    last_pi = pi[0];
                    last_krec = k_rec;
                    last_prob = max_phi;

                    // Evaluate fit
                    if ( (fabs(delpi) < fit_thresh_dpi) &&
                         (fabs(delk) < fit_thresh_dk || fabs(delk/k_rec) < fit_thresh_drelk) ) {

                        finish_fit = true;
                    }
                }  // end parameter fitting loop

/*
                fprintf(
                    pf,
                    "%s\t%s\t%d\t%.4f\t%.5e\t%d\t%.3f\t%d\t%.5f\t%.5f\t%.5f\n",
                    sample1[isamp].c_str(),
                    sample2[jsamp].c_str(),
                    sum,
                    (discord[isamp][jsamp]),
                    max_phi,
                    iter,
                    k_rec,
                    ntrans,
                    seq_ibd / (seq_ibd + seq_dbd),
                    count_ibd_fb / (count_ibd_fb + count_dbd_fb),
                    double(count_ibd_vit) / (count_ibd_vit + count_dbd_vit)
                    );
*/
                freq_fn(
                    sample1[isamp],
                    sample2[jsamp],
                    sum,
                    discord[isamp][jsamp],
                    max_phi,
                    iter,
                    k_rec,
                    ntrans,
                    seq_ibd / (seq_ibd + seq_dbd),
                    count_ibd_fb / (count_ibd_fb + count_dbd_fb),
                    double(count_ibd_vit) / (count_ibd_vit + count_dbd_vit)
                );
            }   // end if use pair

            ++ipair;

            if ((ipair%10000) == 0) {Rprintf("Starting pair %d\n", ipair);}
        }
    }
}
