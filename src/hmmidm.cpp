#include <fstream>
#include <string>
#include <vector>

#include "hmmidm.h"
#include "refactor.h"

void parse_parameters(
    Rcpp::List& param_list,
    bool& freq_flag1, std::string& freq_file1,
    bool& freq_flag2, std::string& freq_file2,
    bool& bflag, std::string& bad_file,
    bool& gflag, std::string& good_file,
    bool& mflag, int& niter,
    bool& nflag, double& k_rec_max,
    bool& iflag1, std::string& data_file1,
    bool& iflag2, std::string& data_file2,
    bool& oflag, std::string& out_filebase);

//' @title
//' hmmidm
//' @description
//' hmmIBD implementation from https://github.com/glipsnort/hmmIBD
//'
//' @param param_list A list of parameters created with \code{hmm_idm}
//'
//' @details
//' \code{hmmidm_cpp} implements hidden Markov model for detecting segments
//'   of shared ancestry (identity by descent) in genetic sequence data.
//'
//' @export
// [[Rcpp::export]]
int hmmidm_cpp(Rcpp::List param_list) {

    /* User-settable parameters */
    double eps         = Rcpp::as<double>(param_list["eps"]);         // error rate in genotype calls
    int min_inform     = Rcpp::as<int>(param_list["min_inform"]);     // minimum number of informative sites in a pairwise
    //  comparison (those w/ minor allele)
    double min_discord = Rcpp::as<double>(param_list["min_discord"]); // minimum discordance in comparison; set > 0 to skip identical pairs
    double max_discord = Rcpp::as<double>(param_list["max_discord"]); // set < 1 to skip unrelated pairs
    int nchrom         = Rcpp::as<int>(param_list["nchrom"]);         // 14 for falciparum
    int min_snp_sep    = Rcpp::as<int>(param_list["min_snp_sep"]);    // skip next snp(s) if too close to last one; in bp
    double rec_rate    = Rcpp::as<double>(param_list["rec_rate"]);    // 7.4e-5 cM/bp or 13.5 kb/cM Miles et al, Genome Res 26:1288-1299 (2016)
    //  const double rec_rate = 5.8e-7;   // 5.8e-5 cM/bp, or 17kb/cM

    /* end user-settable parameters */
    const double fit_thresh_dpi   = .001;
    const double fit_thresh_dk    = .01;
    const double fit_thresh_drelk = .001;

    double k_rec_init = 1.0;          // starting value for N generation parameter
    double k_rec_max = 0;  // working and max value for same
    const int max_all = 8;
    int niter = 5;    // maximum number of iterations of fit; can be overriden by -m
    int max_snp = 30000;
    std::vector<std::string> sample1;
    std::vector<std::string> sample2;
    std::vector<std::string> bad_samp;
    std::vector<std::pair<std::string, std::string>> good_pairs;
    int nsample1=0, nsample2=0, chr, all;
    int **geno1, **geno2, chr2, pos2, majall, npair_report;
    double **discord, fmean;
    double **freq1=nullptr, **freq2=nullptr, *ffreq1=nullptr, *ffreq2=nullptr;
    double pinit[2], pi[2];

    double maxfreq;
    FILE *outf=nullptr, *pf=nullptr;
    int *diff=nullptr, *same_min=nullptr, *allcount1, *allcount2;
    std::vector<bool> use_sample1;
    std::vector<bool> use_sample2;
    int nsample_use2;
    int nsample_use1, nsnp, npair, isnp, *pos;
    int *nmiss_bypair=nullptr, totall1, totall2, maxlen;
    std::vector<int> start_chr;
    std::vector<int> end_chr;
    bool **use_pair = nullptr;
    int *nall=nullptr, nuse_pair=0;
    bool killit;
    int ntri=0, ex_all=0;
    int prev_chrom, nskipped=0;

    std::string freq_file1;
    std::string freq_file2;
    std::string bad_file;
    std::string good_file;
    std::string data_file1;
    std::string data_file2;
    std::string out_filebase;
    std::string file;

    bool freq_flag1 = false;
    bool freq_flag2 = false;
    bool bflag      = false;
    bool gflag      = false;
    bool mflag      = false;
    bool nflag      = false;
    bool iflag1     = false;
    bool iflag2     = false;
    bool oflag      = false;

    pinit[0] = 0.5;  // flat prior
    pinit[1] = 0.5;

    char usage_string[512];
    strcpy(usage_string, "Usage: hmm -i <input file, pop1> -o <output filename> [-I <input file, pop2>] \n");
    strcat(usage_string, "[-m <max fit iter>] [-f <allele freq file, pop1>] [-F <allele freq file, pop2>]\n");
    strcat(usage_string,  "[-b <file with samples to skip>] [-n <max N generation>]");
    strcat(usage_string, "  [-g <file with sample pairs to use>]\n");

    opterr = 0;

    parse_parameters(
        param_list,
        freq_flag1, freq_file1,
        freq_flag2, freq_file2,
        bflag, bad_file,
        gflag, good_file,
        mflag, niter,
        nflag, k_rec_max,
        iflag1, data_file1,
        iflag2, data_file2,
        oflag, out_filebase);

    if (freq_flag2 && !iflag2) {
        REprintf("Inconsistent options: frequency file for 2nd population specified");
        REprintf(" with no data file for that population\n");
        Rcpp::stop("");
    }

    allcount1 = new int[max_all + 1];

    nall = new int[max_snp];
    freq1 = new double*[max_snp];
    pos = new int[max_snp];

    start_chr.resize(nchrom+1, INT_MAX);
    start_chr[0] = 0;
    end_chr.resize(nchrom+1, INT_MIN);
    end_chr[0] = 0;

    for (isnp = 0; isnp < max_snp; ++isnp) {
        freq1[isnp] = new double[max_all+1];
    }
    ffreq1 = new double[max_all+1];
    if (iflag2) {
        allcount2 = new int[max_all+1];
        freq2 = new double*[max_snp];
        for (isnp = 0; isnp < max_snp; ++isnp) {
            freq2[isnp] = new double[max_all+1];
        }
        ffreq2 = new double[max_all+1];
    }
    else {
        allcount2 = allcount1;
        freq2 = freq1;
        ffreq2 = nullptr;
    }

    // clorton - input files, frequency files, good samples, bad samples, output files

    std::vector<std::string> file1_data;
    std::vector<std::string> file2_data;
    read_text_file(data_file1, file1_data);
    if ( iflag2 ) read_text_file(data_file2, file2_data);

    std::vector<std::string> freq1_data;
    std::vector<std::string> freq2_data;
    if ( freq_flag1 ) read_text_file(freq_file1, freq1_data);
    if ( freq_flag2 ) read_text_file(freq_file2, freq2_data);

    if (gflag) gflag = load_good_samples(good_file, good_pairs);
    if (bflag) bflag = load_bad_samples(bad_file, bad_samp);

    outf = open_output_file(out_filebase + ".hmm.txt");
    pf = open_output_file(out_filebase + ".hmm_fract.txt");

    const auto& header1 = file1_data[0];
    nsample1 = split(header1, '\t').size() - 2; // Ignore first two columns, 'chrom' and 'pos'

    geno1 = new int*[nsample1];
    use_sample1.resize(nsample1, true);
    discord = new double*[nsample1];
    use_pair = new bool*[nsample1];  // nsamp1 x nsamp2
    for (int isamp = 0; isamp < nsample1; ++isamp) {
        geno1[isamp] = new int[max_snp];
    }
    prev_chrom = -1;

    nsample_use1 = store_sample_names(header1, sample1, bad_samp, use_sample1);

    if ( iflag2 ) {

        const auto& header2 = file2_data[0];
        nsample2 = split(header2, '\t').size() - 2; // Ignore first two columns, 'chrom' and 'pos'

        geno2 = new int*[nsample2];
        use_sample2.resize(nsample2, true);
        for (int isamp = 0; isamp < nsample2; ++isamp) {
            geno2[isamp] = new int[max_snp];
        }

        nsample_use2 = store_sample_names(header2, sample2, bad_samp, use_sample2);
    }
    else {
        sample2 = sample1;
        nsample2 = nsample1;
        geno2 = geno1;
        use_sample2 = use_sample1;
        nsample_use2 = nsample_use1;
    }

    for (int isamp = 0; isamp < nsample1; ++isamp) {
        use_pair[isamp] = new bool[nsample2];
        discord[isamp] = new double[nsample2];
    }

    display_parameters(niter,
                      min_snp_sep,
                      min_inform,
                      min_discord,
                      max_discord,
                      eps,
                      iflag2,
                      data_file1,
                      data_file2,
                      freq_flag1,
                      freq_file1,
                      freq_flag2,
                      freq_file2,
                      nflag,
                      k_rec_max,
                      nsample_use1,
                      nsample_use2,
                      nsample1,
                      nsample2,
                      npair,
                      npair_report);

    nmiss_bypair = new int[npair];  memset(nmiss_bypair, 0, npair*sizeof(int));
    same_min = new int[npair];      memset(same_min,     0, npair*sizeof(int));
    diff = new int[npair];          memset(diff,         0, npair*sizeof(int));

    determine_pairs_to_compare(nsample1, iflag2, nsample2, gflag, use_pair,
                               good_pairs, sample1, sample2);

    nsnp = 0;

    for ( int iline = 1; iline < file1_data.size(); ++iline ) {

        // See if we need more space...
        if (nsnp == max_snp) {
            increase_storage(nall, max_snp, pos, nsample1, geno1, freq1, max_all, iflag2, nsample2, geno2, freq2);
        }

        totall1 = totall2 = 0;
        killit = false;

        memset(allcount1, 0, (max_all+1)*sizeof(int));
        memset(allcount2, 0, (max_all+1)*sizeof(int));

        // Parse line, pop1
        const std::string& input_line = file1_data[iline];
        parse_input1_line(input_line, chr, pos, nsnp, prev_chrom, min_snp_sep, nskipped, killit,
                         all, max_all, ex_all, geno1, use_sample1, allcount1, totall1);

        // If allele frequencies are not supplied separately and we didn't observe any alleles, cannot determine frequency.
        if (!freq_flag1 && totall1 == 0) {
            killit = true;
        }

        // Parse line, pop2
        if (iflag2) {

            const std::string& input_line = file2_data[iline];
            parse_input2_line(input_line, chr2, pos2, pos, nsnp, chr, all, max_all, killit, ex_all, geno2, use_sample2, allcount2, totall2);

            // If allele frequencies are not supplied separately and we didn't observe any alleles, cannot determine frequency.
            if (!freq_flag2 && totall2 == 0) {
                killit = true;
            }
        }

        if (chr > nchrom) {
            killit = true;
        }

        // if reading freqs from file, read one (pop1)
        if (freq_flag1) {
            const std::string& input_line = freq1_data[iline-1];  // frequency file does not have a header
            parse_frequency_line(input_line, ffreq1, max_all, chr, pos[nsnp]);
        }

        // if reading freqs from file, read one (pop2)
        if (freq_flag2) {
            const std::string& input_line = freq2_data[iline-1];  // frequency file does not have a header
            parse_frequency_line(input_line, ffreq2, max_all, chr, pos[nsnp]);
        }

        nall[nsnp] = 0;
        majall = -1;
        maxfreq = 0;

        // process this variant -- calculate allele frequencies for each pop
        calculate_allele_frequencies(max_all, freq_flag1, freq1[nsnp], ffreq1, allcount1, totall1,
                                    iflag2, freq_flag2, freq2[nsnp], ffreq2, allcount2, totall2,
                                    fmean, maxfreq, majall, nall[nsnp]);

        if (killit) {continue;}

        prev_chrom = chr;

        if (nall[nsnp] > 2) {
            ++ntri;
        }

        // Tabulate differences by pair for calculating discordance
        tabulate_differences_by_pair_for_discordance(nsample1, use_sample1, iflag2, nsample2, use_sample2,
                                                    geno1, geno2, majall, same_min, diff, nmiss_bypair, nsnp);

        if (nsnp < start_chr[chr]) { start_chr[chr] = nsnp; }
        if (nsnp > end_chr[chr])   { end_chr[chr]   = nsnp; }

        ++nsnp;
    }

    Rprintf("Variants skipped for spacing: %d\n", nskipped);
    if (ex_all > 0) {
        Rprintf("Variants with too many alleles: %d\n", ex_all);
    }
    Rprintf("%d variants used,\t%d with >2 alleles\n", nsnp, ntri);

    filter_pairs_by_discordance_and_markers(
        nsample1, use_sample1, iflag2, nsample2, use_sample2, diff, same_min, discord, use_pair, max_discord, min_inform, min_discord, nuse_pair
    );

    maxlen = 0;
    for (chr = 1; chr <= nchrom; ++chr) {
        int length = end_chr[chr] - start_chr[chr] + 1;
        if (length > maxlen) {
            maxlen = length;
        }
    }

    fprintf(outf, "sample1\tsample2\tchr\tstart\tend\tdifferent\tNsnp\n");
    fprintf(pf, "sample1\tsample2\tN_informative_sites\tdiscordance\tlog_p\tN_fit_iteration\tN_generation");
    fprintf(pf, "\tN_state_transition\tseq_shared_best_traj\tfract_sites_IBD\tfract_vit_sites_IBD\n");

    output1_t output_fn1 = [outf](const std::string& sample1, const std::string& sample2, int chromosome, int position) {
        fprintf(outf, "%s\t%s\t%d\t%d", sample1.c_str(), sample2.c_str(), chromosome, position);
    };
    output2_t output_fn2 = [outf](int position, int different, int nsnp) {
        fprintf(outf, "\t%d\t%d\t%d\n", position, different, nsnp);
    };
    output3_t freq_fn = [pf](const std::string& sample1, const std::string& sample2,
                        int num_sites, double discord, double log_p, int fit_iteration,
                        double generation, int state_transition, double seq_shared_best_traj,
                        double fract_sites_ibd, double fract_vit_sites_ibd) {
        fprintf(pf, "%s\t%s\t%d\t%.4f\t%.5e\t%d\t%.3f\t%d\t%.5f\t%.5f\t%.5f\n",
                sample1.c_str(), sample2.c_str(),
                num_sites,
                discord,
                log_p,
                fit_iteration,
                generation,
                state_transition,
                seq_shared_best_traj,
                fract_sites_ibd,
                fract_vit_sites_ibd
        );
    };

    calculate_hmm_ibd(
        maxlen,
        nsample1,
        use_sample1,
        iflag2,
        nsample2,
        use_sample2,
        diff,
        same_min,
        use_pair,
        pi,
        pinit,
        k_rec_init,
        niter,
        nchrom,
        start_chr,
        end_chr,
        pos,
        geno1,
        geno2,
        eps,
        nall,
        freq1,
        freq2,
        rec_rate,
        sample1,
        sample2,
        output_fn1,
        output_fn2,
        freq_fn,
        nflag,
        k_rec_max,
        fit_thresh_dpi,
        fit_thresh_dk,
        fit_thresh_drelk,
        discord
    );

    fclose(pf);
    fclose(outf);

    return(0);
}

void parse_parameters(
    Rcpp::List& param_list,
    bool& freq_flag1, std::string& freq_file1,
    bool& freq_flag2, std::string& freq_file2,
    bool& bflag, std::string& bad_file,
    bool& gflag, std::string& good_file,
    bool& mflag, int& niter,
    bool& nflag, double& k_rec_max,
    bool& iflag1, std::string& data_file1,
    bool& iflag2, std::string& data_file2,
    bool& oflag, std::string& out_filebase)
{
    if (param_list["f"] != R_NilValue){
        freq_flag1 = true;
        freq_file1 = (const char*)param_list["f"];
        Rprintf("freq_file1 = '%s'\n", freq_file1.c_str());
    }
    if (param_list["F"] != R_NilValue){
        freq_flag2 = true;
        freq_file2 = (const char*)param_list["F"];
        Rprintf("freq_file2 = '%s'\n", freq_file2.c_str());
    }
    if (param_list["b"] != R_NilValue){
        bflag = true;
        bad_file = (const char*)param_list["b"];
        Rprintf("bad_file = '%s'\n", bad_file.c_str());
    }
    if (param_list["g"] != R_NilValue){
        gflag = true;
        good_file = (const char*)param_list["g"];
        Rprintf("good_file = '%s'\n", good_file.c_str());
    }
    if (param_list["m"] != R_NilValue){
        mflag = true;
        niter = Rcpp::as<int>(param_list["m"]);
    }
    if (param_list["n"] != R_NilValue){
        nflag = true;
        k_rec_max = Rcpp::as<double>(param_list["n"]);
    }
    if (param_list["i"] != R_NilValue){
        iflag1 = true;
        data_file1 = (const char*)param_list["i"];
        Rprintf("data_file1 = '%s'\n", data_file1.c_str());
    }
    if (param_list["I"] != R_NilValue){
        iflag2 = true;
        data_file2 = (const char*)param_list["I"];
        Rprintf("data_file2 = '%s'\n", data_file2.c_str());
    }
    if (param_list["o"] != R_NilValue){
        oflag = true;
        out_filebase = (const char*)param_list["o"];
        Rprintf("out_filebase = '%s'\n", out_filebase.c_str());
    }
}

using namespace Rcpp;

//' @title
//' IdmRcppTest
//' @description
//' IDM test for calling from/returning data to R
//'
//' @param dfin A dataframe of inputs \code{hmm_idm}
//'
//' @details
//' \code{IdmRcppTest} ... implements hidden Markov model for detecting segments
//'   of shared ancestry (identity by descent) in genetic sequence data.
//'
//' @export
// [[Rcpp::export]]
DataFrame IdmRcppTest(const DataFrame& dfin)
{
    IntegerVector chromosomes = dfin["chromosome"];
    IntegerVector positions = dfin["position"];

    // create a list with n slots
    int n = 13, nobs = 42;
    List res(n);
    CharacterVector list_names = CharacterVector(n);
  
    // populate the list elements
    for (int i = 0; i < n; i++) {
        list_names[i] = "col_" + std::to_string(i);
        res[i] = rnorm(nobs);
    }
  
    // set the names for the list elements
    res.names() = list_names;
    return res;
}
