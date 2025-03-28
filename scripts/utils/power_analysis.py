import numpy as np
import scipy.stats as st

# Helper function to compute the sample size
# needed to power a comparison test for marker gene
# detection.

# The default test is the same as the default
# test used by sc.get.rank_genes_groups_df in the function
# 'get_marker_genes', i.e. 'wilcoxon'.
# To do: add power analyses for other tests.
# The significance level should be the same as the one
# used for the function 'get_marker_genes', after
# applying any multiple-comparisons corrections that are
# also applied in that function.
def get_sample_size(test_type='wilcoxon',
        min_l2fc=1,
        adj_p_val=0.05,
        target_power=0.80,
        rand_seed=0):
    """
    Perform a simulation for same-sized samples
    (if one is larger, power increases) from two
    distributions that differ only by their means
    (a conservative case of difference), to determine
    the sample size required to detect the specified
    mean difference at the specified power.

    Here, we assume that the scaled expression is
    roughly lognormal (i.e. a right-tailed bell-like
    shape with cutoff at 0), since this matches a conservative
    sketch of gene expression after scaling.
    We then work with log1p expression values, which would
    give 2 normal underlying distributions.

    A zero-centered, scaled lognormal distribution would have mean 1 (mu=0 plus 1
    for log1p) and standard deviation 1.
    The log would then result in a normal distribution
    with mu=ln(1/sqrt(2))=-0.347 and std.=ln(2)=0.693.

    We set the standard deviation for both underlying expression
    distributions to std=0.693; we set the mean difference
    to match the specified log2 fold change.
    """
    # Initialize random number generator
    numpy_random_state_experiment = np.random.RandomState(seed=rand_seed)
    # Set sample size parameters to test
    ssize_to_test = 10
    ssize_increment = 2
    n_iters_per_test = 10000
    power_curr = 0
    # Perform power analysis for the specified test
    if test_type == 'wilcoxon':
        # Perform a simulation for same-sized
        # distributions (if one is larger power increases)
        # Here, we assume that the scaled expression
        # distribution roughly lognormal (i.e. a right-tailed
        # bell-like shape with cutoff at 0), since this
        # matches a conservative sketch of gene expression
        # after scaling. This dist. would have mean 1 (mu=0
        # plus 1 for log1p) and std. 1.
        # The log would then result in a normal distribution
        # with mu=ln(1/sqrt(2))=-0.347, std=ln(2)=0.693.
        # Distributions with equal std. should be hardest
        # to distinguish so we set std=0.693 here.
        # We set the mean difference to match the specified
        # log2 fold change. Here, mu2=(2^(l2fc))*mu1.
        mu1 = 1
        mu2 = 2**min_l2fc
        std12 = 0.693
        power_achieved = False
        while not power_achieved:
            n_sub_ps = 0
            # Perform repeated Wilcoxon tests
            for iter in range(n_iters_per_test):
                # Pull samples
                s1 = numpy_random_state_experiment.normal(
                        loc=mu1,
                        scale=std12,
                        size=ssize_to_test)
                s2 = numpy_random_state_experiment.normal(
                        loc=mu2,
                        scale=std12,
                        size=ssize_to_test)
                # Perform Wilcoxon test
                result = st.wilcoxon(s1,s2)
                if result.pvalue < adj_p_val:
                    n_sub_ps += 1
            # Compute the fraction of true positives
            power_curr = n_sub_ps*1.0/n_iters_per_test
            # Compare power with target
            if power_curr >= target_power:
                # Sample size is sufficient; break
                power_achieved = True
            else:
                # Sample size is insufficient;
                # increment sample size and re-test
                ssize_to_test += ssize_increment
    return ssize_to_test, power_curr

