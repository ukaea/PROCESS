'''
Master script to generate input, run calculations and collect results of the bt scan
06.2025 Walkowiak
'''
from stellarator_test.autostart import generate_input, run_cases, collect_results

case_name = 'squid_v1_backward_test'
prefix = 'squid'

# case_name = 'rebuild'
# prefix = 'rebuild'

# case_name = 'updated_beta5'
# prefix = 'updated'

generate_input.main(case_name, prefix = prefix, B_min=6.0, B_max=6.5, clean_start=True)
run_cases.main(case_name, prefix = prefix)
collect_results.main(case_name, prefix = prefix)
