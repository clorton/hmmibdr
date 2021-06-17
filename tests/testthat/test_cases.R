# install.packages("devtools")
# devtools::load_all("/Users/christopher/Coding/hmmibdr")
# test_file("/Users/christopher/Coding/hmmibdr/tests/testthat/test_cases.R")

run_scenario <- function(name, params) {
    output = tempfile(name)
    params$output_file = output
    result = hmm_ibd(
        input_file           = params$input_file,
        allele_freqs         = params$allele_freqs,
        genotypes_sec_pop    = params$genotypes_sec_pop,
        allele_freqs_sec_pop = params$allele_freqs_sec_pop,
        analysis_ids         = params$analysis_ids,
        exclude_ids          = params$exclude_ids,
        output_file          = params$output_file)
    output
}

compare_files <- function(output, baseline) {
    # HMM
    hmm <- read.csv(paste(output, ".hmm.txt", sep=""), sep = "\t")
    # FRAC
    frac <- read.csv(paste(output, ".hmm_fract.txt", sep=""), sep = "\t")
}

test_that("unittests are run", {
    pop1 = system.file("extdata/pf3k_Cambodia_13.txt", package="hmmibdr")
    freq1 = system.file("extdata/freqs_pf3k_Cambodia_13.txt", package="hmmibdr")
    pop2 = system.file("extdata/pf3k_Ghana_13.txt", package="hmmibdr")
    freq2 = system.file("extdata/freqs_pf3k_Ghana_13.txt", package="hmmibdr")
    good_pairs = system.file("extdata/good_pairs.txt", package="hmmibdr")
    bad_samples = system.file("extdata/bad_samples.txt", package="hmmibdr")

	# single population
	# ./hmmIBD -i samp_data/pf3k_Cambodia_13.txt -o test_output/single_population
    expect_output(name <- run_scenario("single_population", list(input_file=pop1)))
    print(name)
    expect_equal(read.csv(paste(name, ".hmm.txt", sep=""), sep="\t"), read.csv(system.file("baseline/single_population.hmm.txt", package="hmmibdr"), sep="\t"))
    expect_equal(read.csv(paste(name, ".hmm_fract.txt", sep=""), sep="\t"), read.csv(system.file("baseline/single_population.hmm_fract.txt", package="hmmibdr"), sep="\t"))

	# single population with frequencies ~ samp_data/output_Cambodia.*
	# ./hmmIBD -i samp_data/pf3k_Cambodia_13.txt -f samp_data/freqs_pf3k_Cambodia_13.txt -o test_output/single_population_with_frequencies
    expect_output(name <- run_scenario("single_population_with_frequencies", list(input_file=pop1, allele_freqs=freq1)))
    print(name)
    expect_equal(read.csv(paste(name, ".hmm.txt", sep=""), sep="\t"), read.csv(system.file("baseline/single_population_with_frequencies.hmm.txt", package="hmmibdr"), sep="\t"))
    expect_equal(read.csv(paste(name, ".hmm_fract.txt", sep=""), sep="\t"), read.csv(system.file("baseline/single_population_with_frequencies.hmm_fract.txt", package="hmmibdr"), sep="\t"))

	# two populations
	# ./hmmIBD -i samp_data/pf3k_Cambodia_13.txt -I samp_data/pf3k_Ghana_13.txt -o test_output/two_populations
    expect_output(name <- run_scenario("two_populations", list(input_file=pop1, genotypes_sec_pop=pop2)))
    print(name)
    expect_equal(read.csv(paste(name, ".hmm.txt", sep=""), sep="\t"), read.csv(system.file("baseline/two_populations.hmm.txt", package="hmmibdr"), sep="\t"))
    expect_equal(read.csv(paste(name, ".hmm_fract.txt", sep=""), sep="\t"), read.csv(system.file("baseline/two_populations.hmm_fract.txt", package="hmmibdr"), sep="\t"))

	# two populations, allele frequencies for population 1
	# ./hmmIBD -i samp_data/pf3k_Cambodia_13.txt -f samp_data/freqs_pf3k_Cambodia_13.txt -I samp_data/pf3k_Ghana_13.txt -o test_output/two_populations_one_frequency
    expect_output(name <- run_scenario("two_populations_one_frequency", list(input_file=pop1, allele_freqs=freq1, genotypes_sec_pop=pop2)))
    print(name)
    expect_equal(read.csv(paste(name, ".hmm.txt", sep=""), sep="\t"), read.csv(system.file("baseline/two_populations_one_frequency.hmm.txt", package="hmmibdr"), sep="\t"))
    expect_equal(read.csv(paste(name, ".hmm_fract.txt", sep=""), sep="\t"), read.csv(system.file("baseline/two_populations_one_frequency.hmm_fract.txt", package="hmmibdr"), sep="\t"))

	# two populations, allele frequencies for both populations ~ samp_data/output_Cambodia_Ghana
	# ./hmmIBD -i samp_data/pf3k_Cambodia_13.txt -f samp_data/freqs_pf3k_Cambodia_13.txt -I samp_data/pf3k_Ghana_13.txt -F samp_data/freqs_pf3k_Ghana_13.txt -o test_output/two_populations_with_frequencies
    expect_output(name <- run_scenario("two_populations_with_frequencies", list(input_file=pop1, allele_freqs=freq1, genotypes_sec_pop=pop2, allele_freqs_sec_pop=freq2)))
    print(name)
    expect_equal(read.csv(paste(name, ".hmm.txt", sep=""), sep="\t"), read.csv(system.file("baseline/two_populations_with_frequencies.hmm.txt", package="hmmibdr"), sep="\t"))
    expect_equal(read.csv(paste(name, ".hmm_fract.txt", sep=""), sep="\t"), read.csv(system.file("baseline/two_populations_with_frequencies.hmm_fract.txt", package="hmmibdr"), sep="\t"))

	# single population with good pairs file
	# ./hmmIBD -i samp_data/pf3k_Cambodia_13.txt -o test_output/single_population_good_pairs
    expect_output(name <- run_scenario("single_population_good_pairs", list(input_file=pop1, analysis_ids=good_pairs)))
    print(name)
    expect_equal(read.csv(paste(name, ".hmm.txt", sep=""), sep="\t"), read.csv(system.file("baseline/single_population_good_pairs.hmm.txt", package="hmmibdr"), sep="\t"))
    expect_equal(read.csv(paste(name, ".hmm_fract.txt", sep=""), sep="\t"), read.csv(system.file("baseline/single_population_good_pairs.hmm_fract.txt", package="hmmibdr"), sep="\t"))

	# single population with bad samples file
	# ./hmmIBD -i samp_data/pf3k_Cambodia_13.txt -o test_output/single_population_bad_samples
    expect_output(name <- run_scenario("single_population_bad_samples", list(input_file=pop1, exclude_ids=bad_samples)))
    print(name)
    expect_equal(read.csv(paste(name, ".hmm.txt", sep=""), sep="\t"), read.csv(system.file("baseline/single_population_bad_samples.hmm.txt", package="hmmibdr"), sep="\t"))
    expect_equal(read.csv(paste(name, ".hmm_fract.txt", sep=""), sep="\t"), read.csv(system.file("baseline/single_population_bad_samples.hmm_fract.txt", package="hmmibdr"), sep="\t"))
})
