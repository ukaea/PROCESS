'''
Master script to generate input, run calculations and collect results of the bt scan
06.2025 Walkowiak
'''
from stellarator_test.autostart import generate_input, run_cases, collect_results

case_name = 'updated_beta5'
prefix = 'updated'

generate_input.main(case_name, prefix = prefix)
run_cases.main(case_name, prefix = prefix, skip_calculated=True)
collect_results.main(case_name, prefix = prefix)
