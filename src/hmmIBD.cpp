#include <fstream>
#include <string>
#include <vector>

#include "hmmIBD.h"
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
//' hmmIBD
//' @description
//' hmmIBD implementation from https://github.com/glipsnort/hmmIBD
//'
//' @param param_list A list of parameters created with \code{hmm_ibd}
//'
//' @details
//' \code{hmmibd_c} implements hidden Markov model for detecting segments
//'   of shared ancestry (identity by descent) in genetic sequence data.
//'
//' @export
// [[Rcpp::export]]
int hmmibd_c(Rcpp::List param_list) {

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
    double k_rec, k_rec_max = 0;  // working and max value for same
    const int max_all = 8;
    int niter = 5;    // maximum number of iterations of fit; can be overriden by -m
    int max_snp = 30000;
    std::vector<std::string> sample1;
    std::vector<std::string> sample2;
    std::vector<std::string> bad_samp;
    std::vector<std::pair<std::string, std::string>> good_pairs;
    int nsample1=0, nsample2=0, chr, sum, all, js, snp_ind;
    int **geno1, **geno2, chr2, pos2, majall, npair_report;
    double **discord, pright, seq_ibd_fb=0, seq_dbd_fb=0, p_ibd, fmean, fmeani, fmeanj;
    double **freq1=nullptr, **freq2=nullptr, *ffreq1=nullptr, *ffreq2=nullptr, xisum, xi[2][2], trans_pred, trans_obs;
    double *phi[2], pinit[2], pi[2], *b[2], a[2][2], ptrans, *alpha[2], *beta[2], *scale;
    double maxval, max_phi=0, max_phiL, seq_ibd, seq_dbd, count_ibd_fb, count_dbd_fb;
    double gamma[2], last_pi=0, last_prob=0, last_krec=0, delpi, delprob, delk, maxfreq;
    FILE *outf=nullptr, *pf=nullptr;
    int *diff=nullptr, *same_min=nullptr, *allcount1, *allcount2;
    std::vector<bool> use_sample1;
    std::vector<bool> use_sample2;
    int *traj, add_seq, nsample_use2;
    int nsample_use1, nsnp, npair, isnp, chrlen, *pos, *psi[2], max;
    int *nmiss_bypair=nullptr, totall1, totall2, is, maxlen;
    std::vector<int> start_chr;
    std::vector<int> end_chr;
    bool **use_pair = nullptr;
    int *nall=nullptr, nuse_pair=0, gi, gj, delpos;
    bool killit;
    int ntri=0, start_snp, ex_all=0, last_snp;
    int fpos=0, fchr=0, iter, ntrans;
    bool finish_fit;
    int prev_chrom, nskipped=0, nsite;
    int count_ibd_vit, count_dbd_vit;

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

    start_chr.resize(nchrom+1, 10000000);
    start_chr[0] = 0;
    end_chr.resize(nchrom+1, -1);
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
    fprintf(outf, "sample1\tsample2\tchr\tstart\tend\tdifferent\tNsnp\n");

    pf = open_output_file(out_filebase + ".hmm_fract.txt");
    fprintf(pf, "sample1\tsample2\tN_informative_sites\tdiscordance\tlog_p\tN_fit_iteration\tN_generation");
    fprintf(pf, "\tN_state_transition\tseq_shared_best_traj\tfract_sites_IBD\tfract_vit_sites_IBD\n");

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

    nmiss_bypair = new int[npair];
    same_min = new int[npair];
    diff = new int[npair];

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

        if (!freq_flag1 && totall1 == 0) {
            killit = true; // no valid calls to calculate frequency
        }

        // Parse line, pop2
        if (iflag2) {

            const std::string& input_line = file2_data[iline];
            parse_input2_line(input_line, chr2, pos2, pos, nsnp, chr, all, max_all, killit, ex_all, geno2, use_sample2, allcount2, totall2);

            if (!freq_flag2 && totall2 == 0) {
                killit = true; // no valid calls to calculate frequency
            }
        }

        if (chr > nchrom) {
            killit = true;
        }

        // if reading freqs from file, read one (pop1)
        if (freq_flag1) {
            const std::string& input_line = freq1_data[iline-1];  // frequency file does not have a header
            parse_frequency_line(input_line, ffreq1, max_all, fpos, fchr, chr, pos, nsnp);
        }

        // if reading freqs from file, read one (pop2)
        if (freq_flag2) {
            const std::string& input_line = freq2_data[iline-1];  // frequency file does not have a header
            parse_frequency_line(input_line, ffreq2, max_all, fpos, fchr, chr, pos, nsnp);
        }

        nall[nsnp] = 0;
        majall = -1;
        maxfreq = 0;

        // process this variant -- calculate allele frequencies for each pop
        calculate_allele_frequencies(max_all, freq_flag1, freq1, nsnp, ffreq1, allcount1, totall1,
                                    iflag2, freq_flag2, freq2, ffreq2, allcount2, totall2,
                                    fmean, maxfreq, majall, nall);

        if (killit) {continue;}

        prev_chrom = chr;

        if (nall[nsnp] > 2) {
            ++ntri;
        }

        // Tabulate differences by pair for calculating discordance
        tabulate_differences_by_pair_for_discordance(nsample1, use_sample1, iflag2, nsample2, use_sample2,
                                                    geno1, geno2, majall, same_min, diff, nmiss_bypair, nsnp, chr, start_chr, end_chr);

        ++nsnp;
    }

    Rprintf("Variants skipped for spacing: %d\n", nskipped);
    if (ex_all > 0) {
        Rprintf("Variants with too many alleles: %d\n", ex_all);
    }
    Rprintf("%d variants used,\t%d with >2 alleles\n", nsnp, ntri);

    filter_pairs_by_discordance_and_markers(
        nsample1, use_sample1, iflag2, nsample2, use_sample2, sum, diff, same_min, discord, use_pair, max_discord, min_inform, min_discord, nuse_pair
    );

    maxlen = 0;
    for (chr = 1; chr <= nchrom; ++chr) {
        int length = end_chr[chr] - start_chr[chr] + 1;
        if (length > maxlen) {
            maxlen = length;
        }
    }

    for (is = 0; is < 2; ++is) {
        psi[is] = new int[maxlen];
        phi[is] = new double[maxlen];
        b[is] = new double[maxlen];
        alpha[is] = new double[maxlen];
        beta[is] = new double[maxlen];
    }

    scale = new double[maxlen];
    traj = new int[maxlen];

    {
        int ipair = 0;
        for (int isamp = 0; isamp < nsample1; ++isamp) {

            if (!use_sample1[isamp]) {continue;}

            int jstart = (iflag2 ? 0 : isamp+1);
            for (int jsamp = jstart; jsamp < nsample2; ++jsamp) {

                if (!use_sample2[jsamp]) {continue;}

                sum = diff[ipair] + same_min[ipair];

                if (use_pair[isamp][jsamp]) {

                    last_prob = 0;
                    last_pi = pi[0] = pinit[0];  // initialize with prior
                    pi[1] = pinit[1];
                    last_krec = k_rec = k_rec_init;
                    finish_fit = false;

                    for (iter = 0; iter < niter; ++iter) {

                        trans_obs = trans_pred = ntrans = 0;
                        seq_ibd = seq_dbd = count_ibd_fb = count_dbd_fb = seq_ibd_fb = seq_dbd_fb = 0;
                        count_ibd_vit = count_dbd_vit = 0;
                        max_phi = 0;

                        for (chr = 1; chr <= nchrom; ++chr) {

                            nsite = 0;

                            if (end_chr[chr] < 0) {continue;}

                            chrlen = pos[end_chr[chr]] - pos[start_chr[chr]];

                            for (isnp = start_chr[chr]; isnp <= end_chr[chr]; ++isnp) {

                                snp_ind = isnp - start_chr[chr];
                                gi = geno1[isamp][isnp];
                                gj = geno2[jsamp][isnp];
                                pright = 1 - eps * (nall[isnp] - 1);

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

                                    for (is = 0; is < 2; ++is) {
                                        phi[is][snp_ind] = log(pi[is]) + log(b[is][snp_ind]);
                                        alpha[is][snp_ind] = pi[is] * b[is][snp_ind];
                                        scale[is] = 1.;
                                    }
                                }
                                else {
                                    ptrans = k_rec * rec_rate * (pos[isnp] - pos[isnp-1]);
                                    a[0][1] = 1 - pi[0] - (1 - pi[0]) * exp(-ptrans);
                                    a[1][0] = 1 - pi[1] - (1 - pi[1]) * exp(-ptrans);
                                    a[0][0] = 1 - a[0][1];
                                    a[1][1] = 1 - a[1][0];

                                    for (js = 0; js < 2; ++js) {    // index over state of current snp

                                        maxval = -10000000;
                                        alpha[js][snp_ind] = scale[snp_ind] = 0;

                                        for (is = 0; is < 2; ++is) {    // index over state of previous snp

                                            if (phi[is][snp_ind-1] + log(a[is][js]) > maxval ) {
                                                maxval = phi[is][snp_ind-1] + log(a[is][js]);
                                                psi[js][snp_ind] = is;
                                            }

                                            phi[js][snp_ind] = maxval + log(b[js][snp_ind]);
                                            alpha[js][snp_ind] += alpha[is][snp_ind-1] * a[is][js] * b[js][snp_ind];
                                        }

                                        scale[snp_ind] += alpha[js][snp_ind];
                                    }

                                    for (js = 0; js < 2; ++js) {    // scale alpha to prevent underflow (Rabiner eqns 92)
                                        alpha[js][snp_ind] /= scale[snp_ind];
                                    }
                                }   // end if initializing/continuing
                            }  // end snp loop

                            last_snp = snp_ind;
                            max_phiL = phi[1][snp_ind];

                            if (phi[0][snp_ind] > phi[1][snp_ind]) {max_phiL = phi[0][snp_ind];}

                            max = (phi[1][snp_ind] > phi[0][snp_ind]) ? 1 : 0;
                            traj[snp_ind] = max;
                            max_phi += max_phiL;

                            // Loop backward to calculate betas (backward part of forward-backward)
                            isnp = end_chr[chr];
                            snp_ind = isnp - start_chr[chr];
                            beta[0][snp_ind] = beta[1][snp_ind] = 1;

                            for (isnp = end_chr[chr]-1; isnp >= start_chr[chr]; isnp--) {

                                snp_ind = isnp - start_chr[chr];
                                ptrans = k_rec * rec_rate * (pos[isnp+1] - pos[isnp]);
                                a[0][1] = 1 - pi[0] - (1 - pi[0]) * exp(-ptrans);
                                a[1][0] = 1 - pi[1] - (1 - pi[1]) * exp(-ptrans);
                                a[0][0] = 1 - a[0][1];
                                a[1][1] = 1 - a[1][0];

                                for (is = 0; is < 2; ++is) {    // index over state of current snp

                                    beta[is][snp_ind] = 0;

                                    for (js = 0; js < 2; ++js) {
                                        // index over state of previous snp (= snp+1, since looping backward)
                                        beta[is][snp_ind] += beta[js][snp_ind+1] * a[is][js] * b[js][snp_ind+1] / scale[snp_ind];
                                    }
                                }
                            }

                            // inverse loop for Viterbi
                            for (snp_ind = last_snp-1; snp_ind >= 0; snp_ind--) {
                                traj[snp_ind] = psi[max][snp_ind+1];
                                max = traj[snp_ind];
                            }

                            // tabulate FB results
                            for (snp_ind = 0; snp_ind <= last_snp; ++snp_ind) {

                                p_ibd = alpha[0][snp_ind] * beta[0][snp_ind] /
                                        (alpha[0][snp_ind] * beta[0][snp_ind] + alpha[1][snp_ind] * beta[1][snp_ind]);

                                count_ibd_fb += p_ibd;
                                count_dbd_fb += (1-p_ibd);

                                if (snp_ind < last_snp) {

                                    isnp = snp_ind + start_chr[chr];
                                    delpos = pos[isnp+1] - pos[isnp];
                                    ptrans = k_rec * rec_rate * delpos;
                                    a[0][1] = 1 - pi[0] - (1 - pi[0]) * exp(-ptrans);
                                    a[1][0] = 1 - pi[1] - (1 - pi[1]) * exp(-ptrans);
                                    a[0][0] = 1 - a[0][1];
                                    a[1][1] = 1 - a[1][0];
                                    xisum = 0;

                                    for (is = 0; is < 2; ++is) {
                                        for (js = 0; js < 2; ++js) {
                                            xi[is][js] = alpha[is][snp_ind] * a[is][js] * b[js][snp_ind+1] * beta[js][snp_ind+1];
                                            xisum += xi[is][js];
                                        }
                                    }

                                    for (is = 0; is < 2; ++is) {

                                        gamma[is] = 0;

                                        for (js = 0; js < 2; ++js) {
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

                                    fprintf(outf, "%s\t%s\t%d\t%d", sample1[isamp].c_str(), sample2[jsamp].c_str(), chr, pos[0+start_chr[chr]]);

                                    start_snp = 0;
                                    for (isnp = 1; isnp < end_chr[chr] - start_chr[chr] + 1; ++isnp) {

                                        if (traj[isnp] != traj[isnp-1]) {

                                            ++ntrans;
                                            add_seq = pos[isnp - 1 + start_chr[chr]] - pos[start_snp + start_chr[chr]] + 1;

                                            if (traj[isnp-1] == 0) {seq_ibd += add_seq;}
                                            else {seq_dbd += add_seq;}

                                            // end one and start another
                                            fprintf(outf, "\t%d\t%d\t%d\n",
                                                    pos[isnp - 1 + start_chr[chr]], traj[isnp-1], isnp - start_snp);

                                            fprintf(outf, "%s\t%s\t%d\t%d", sample1[isamp].c_str(),
                                                    sample2[jsamp].c_str(), chr, pos[isnp + start_chr[chr]]);

                                            start_snp = isnp;
                                        }
                                    }

                                    isnp = end_chr[chr] - start_chr[chr];
                                    add_seq = pos[isnp + start_chr[chr]] - pos[start_snp + start_chr[chr]] + 1;

                                    if (traj[isnp] == 0) {seq_ibd += add_seq;}
                                    else {seq_dbd += add_seq;}

                                    fprintf(outf, "\t%d\t%d\t%d\n", pos[end_chr[chr]], traj[isnp], isnp - start_snp + 1);

                                    // Tabulate sites by state for Viterbi trajectory
                                    for (isnp = 0; isnp < end_chr[chr] - start_chr[chr] + 1; ++isnp) {
                                        if (traj[isnp] == 0) {++count_ibd_vit;}
                                        else {++count_dbd_vit;}
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
                        delpi = pi[0] - last_pi;
                        delk = k_rec - last_krec;

                        if (nflag && k_rec > k_rec_max) {delk = k_rec_max - last_krec;}

                        delprob = max_phi - last_prob;
                        last_pi = pi[0];
                        last_krec = k_rec;
                        last_prob = max_phi;

                        // Evaluate fit
                        if ( (fabs(delpi) < fit_thresh_dpi) &&
                            (fabs(delk) < fit_thresh_dk || fabs(delk/k_rec) < fit_thresh_drelk) ) {

                            finish_fit = true;
                        }
                    }  // end parameter fitting loop

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

                }   // end if use pair

                ++ipair;

                if ((ipair%10000) == 0) {Rprintf("Starting pair %d\n", ipair);}
            }
        }
    }

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
