'''
Master script to generate input, run calculations and collect results of the bt scan
06.2025 Walkowiak
'''
from stellarator_test.autostart import generate_input, run_cases, collect_results

case_name = 'helias5_7T'
prefix = 'helias5_7T'

# case_name = 'rebuild'
# prefix = 'rebuild'

# case_name = 'updated_beta5'
# prefix = 'updated'

generate_input.main(case_name, prefix = prefix, B_min=7.0, B_max=7.2)
run_cases.main(case_name, prefix = prefix, skip_calculated=True)
collect_results.main(case_name, prefix = prefix)
